/* This file is a part of the CompLaB program.  AGPL-3.0-or-later.
 * Meile Lab, University of Georgia.  shahram.asgari@uga.edu
*/

/* ================================================================================================
 * complab3d_diagnostics.hh  --  WRITE THE NUMBERS, NOT JUST THE PICTURES
 * ================================================================================================
 *
 *  WHAT THIS REPLACES
 *
 *  Until now CompLB3D wrote .vti snapshots and nothing else. No totals, no porosity history, no
 *  mass balance. Every quantitative statement about a run had to be reconstructed afterwards from
 *  the VTK files by a script the user wrote. This writes the numbers as the run produces them.
 *
 *      <diagnostics>
 *          <enabled>true</enabled>
 *          <summary_csv>summary.csv</summary_csv>
 *          <interval>500</interval>
 *          <conserve>Fe2+FeS</conserve>       <!-- must stay constant -->
 *          <conserve>HS+FeS</conserve>
 *          <tolerance>1e-6</tolerance>
 *      </diagnostics>
 *
 *  One row per interval: iteration, porosity, and the total, mean, minimum and maximum of every
 *  substrate and every microbe -- over OPEN voxels only, so a run that seals pore space is not
 *  compared against a moving denominator.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHY THE CONSERVATION CHECK LIVES HERE AND NOT IN A SCRIPT
 *
 *  Because it can stop the run. A conservation drift means the chemistry is wrong or the time step
 *  is too long, and finding that out at step 200 is worth more than finding it out after a
 *  fortnight of wall-clock. The check prints a warning as soon as the drift exceeds the tolerance,
 *  and the final row records the verdict.
 *
 *  A sum like "Fe2+FeS" names substrates by their <name_of_substrates>. Anything the reaction
 *  network cannot create or destroy belongs here: one mole of Fe2+ removed must appear as one mole
 *  of FeS.
 *
 *  ------------------------------------------------------------------------------------------------
 *  MPI. Every total is a global reduction. The caller supplies the already-reduced values, because
 *  this header deliberately knows nothing about Palabos -- which is what lets the test suite drive
 *  it with synthetic fields and check the arithmetic.
 *
 *  Only rank 0 writes. A CSV opened by 36 ranks at once is a corrupted CSV.
 * ================================================================================================
 */
#ifndef COMPLAB3D_DIAGNOSTICS_HH
#define COMPLAB3D_DIAGNOSTICS_HH

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace complab_diag {

struct FieldStat {
    std::string name;
    double total, mean, minv, maxv;
    FieldStat() : total(0), mean(0), minv(0), maxv(0) {}
};

struct Row {
    long iteration;
    double porosity;
    long openVoxels;
    std::vector<FieldStat> fields;
    Row() : iteration(0), porosity(0), openVoxels(0) {}
};

/* A conserved sum, named by substrate. */
struct ConserveCheck {
    std::string expr;
    std::vector<std::string> terms;
    double first;
    bool haveFirst;
    double worstDrift;
    bool warned;
    bool skipped;
    ConserveCheck() : first(0), haveFirst(false), worstDrift(0), warned(false), skipped(false) {}
};

inline std::vector<std::string> splitPlus(const std::string &s)
{
    std::vector<std::string> out;
    std::string cur;
    for (size_t i = 0; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == '+') {
            std::string t;
            for (size_t k = 0; k < cur.size(); ++k)
                if (cur[k] != ' ' && cur[k] != '\t') t += cur[k];
            if (!t.empty()) out.push_back(t);
            cur.clear();
        } else cur += s[i];
    }
    return out;
}

/* ------------------------------------------------------------------------------------------------ */
class Diagnostics {
public:
    Diagnostics() : enabled_(false), tol_(1e-6), rows_(0), anyFail_(false), master_(true) {}

    void configure(bool enabled, const std::string &path, double tol, bool isMaster)
    {
        enabled_ = enabled;
        path_ = path;
        tol_ = tol;
        master_ = isMaster;
    }

    void addConserve(const std::string &expr)
    {
        ConserveCheck c;
        c.expr = expr;
        c.terms = splitPlus(expr);
        if (!c.terms.empty()) checks_.push_back(c);
    }

    bool active() const { return enabled_; }
    bool failed() const { return anyFail_; }

    /* Returns any message the caller should print. Empty means all is well. */
    std::string record(const Row &r)
    {
        if (!enabled_) return "";
        std::string msg;

        if (master_) {
            if (rows_ == 0) writeHeader(r);
            writeRow(r);
        }
        ++rows_;

        for (size_t c = 0; c < checks_.size(); ++c) {
            ConserveCheck &k = checks_[c];
            double sum = 0.0;
            bool complete = true;
            for (size_t t = 0; t < k.terms.size(); ++t) {
                bool found = false;
                for (size_t f = 0; f < r.fields.size(); ++f)
                    if (r.fields[f].name == k.terms[t]) { sum += r.fields[f].total; found = true; break; }
                if (!found) { complete = false; break; }
            }
            if (!complete) {
                if (!k.warned) {
                    k.warned = true;
                    k.skipped = true;
                    msg += "  [DIAG] <conserve>" + k.expr + "</conserve> names a substrate that "
                           "does not exist; the check is skipped.\n";
                }
                continue;
            }
            if (!k.haveFirst) { k.first = sum; k.haveFirst = true; continue; }

            const double scale = (std::fabs(k.first) > 1e-300) ? std::fabs(k.first) : 1e-300;
            const double drift = std::fabs(sum - k.first) / scale;
            if (drift > k.worstDrift) k.worstDrift = drift;
            if (drift > tol_ && !k.warned) {
                k.warned = true;
                anyFail_ = true;
                char b[512];
                std::snprintf(b, sizeof(b),
                    "  [DIAG] MASS BALANCE: %s drifted %.3e (tolerance %.1e) by step %ld.\n"
                    "  [DIAG]   started %.10g, now %.10g.\n"
                    "  [DIAG]   That is chemistry that does not close, or a time step too long\n"
                    "  [DIAG]   for the rate you specified. It will not fix itself.\n",
                    k.expr.c_str(), drift, tol_, r.iteration, k.first, sum);
                msg += b;
            }
        }

        /* A negative concentration is not physical and is worth saying once, loudly. */
        for (size_t f = 0; f < r.fields.size(); ++f) {
            if (r.fields[f].minv < -1e-9 && !negWarned_.count(r.fields[f].name)) {
                negWarned_.insert(r.fields[f].name);
                anyFail_ = true;
                char b[256];
                std::snprintf(b, sizeof(b),
                    "  [DIAG] %s went negative (%.3e) at step %ld.\n",
                    r.fields[f].name.c_str(), r.fields[f].minv, r.iteration);
                msg += b;
            }
        }
        return msg;
    }

    std::string finalReport() const
    {
        if (!enabled_) return "";
        std::string s = "\n  [DIAG] summary written to " + path_ + "\n";
        char b[256];
        for (size_t c = 0; c < checks_.size(); ++c) {
            /* A check that never ran is reported as SKIPPED, not PASS. A skipped check that
             * looks like a pass is how a mass-balance guarantee quietly becomes worthless. */
            std::snprintf(b, sizeof(b), "  [DIAG] %-24s worst relative drift %.3e  %s\n",
                          checks_[c].expr.c_str(), checks_[c].worstDrift,
                          checks_[c].skipped ? "SKIPPED (unknown substrate)"
                                             : (checks_[c].worstDrift <= tol_ ? "PASS" : "FAIL"));
            s += b;
        }
        if (!checks_.empty())
            s += std::string("  [DIAG] verdict: ") + (anyFail_ ? "FAIL" : "PASS") + "\n";
        return s;
    }

private:
    void writeHeader(const Row &r)
    {
        std::FILE *f = std::fopen(path_.c_str(), "w");
        if (!f) return;
        std::fprintf(f, "# CompLB3D run summary, one row per diagnostic interval\n");
        std::fprintf(f, "# totals and means are over OPEN voxels only\n");
        std::fprintf(f, "iteration,porosity,open_voxels");
        for (size_t i = 0; i < r.fields.size(); ++i) {
            const char *n = r.fields[i].name.c_str();
            std::fprintf(f, ",%s_total,%s_mean,%s_min,%s_max", n, n, n, n);
        }
        std::fprintf(f, "\n");
        std::fclose(f);
    }

    void writeRow(const Row &r)
    {
        std::FILE *f = std::fopen(path_.c_str(), "a");
        if (!f) return;
        std::fprintf(f, "%ld,%.17g,%ld", r.iteration, r.porosity, r.openVoxels);
        for (size_t i = 0; i < r.fields.size(); ++i)
            std::fprintf(f, ",%.17g,%.17g,%.17g,%.17g",
                         r.fields[i].total, r.fields[i].mean, r.fields[i].minv, r.fields[i].maxv);
        std::fprintf(f, "\n");
        std::fclose(f);          // reopened per row on purpose: a run killed by the queue still
                                 // leaves a complete, readable CSV behind
    }

    bool enabled_;
    std::string path_;
    double tol_;
    long rows_;
    bool anyFail_;
    bool master_;
    std::vector<ConserveCheck> checks_;

    struct StrSet {
        std::vector<std::string> v;
        bool count(const std::string &s) const {
            for (size_t i = 0; i < v.size(); ++i) if (v[i] == s) return true;
            return false;
        }
        void insert(const std::string &s) { if (!count(s)) v.push_back(s); }
    };
    StrSet negWarned_;
};

}  // namespace complab_diag

#endif  // COMPLAB3D_DIAGNOSTICS_HH
