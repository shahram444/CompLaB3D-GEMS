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
    /* What computeDensity() cannot see. Palabos's BounceBack and NoDynamics both answer
     * computeDensity() from a stored number and ignore the populations they are holding, so mass
     * momentarily in flight at a wall or inside a grain is absent from `total`. `held` is that
     * mass, measured from the populations themselves. It is reported as its own column, and the
     * conservation check adds it back, so a run is judged on what the lattice actually holds
     * rather than on what the dynamics are willing to say. */
    double held;
    FieldStat() : total(0), mean(0), minv(0), maxv(0), held(0) {}
    /* The conserved quantity: everything the lattice is holding, wherever it is sitting. */
    double conserved() const { return total + held; }
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
    /* [v1.3] How many times the sum was actually COMPARED against its baseline. record()
     * establishes the baseline and returns, so a run producing fewer than two diagnostic rows
     * left worstDrift at 0 and reported PASS having compared nothing -- which is precisely the
     * "a skipped check that looks like a pass" failure the note in finalReport() warns about. */
    long comparisons;
    double worstDrift;
    bool warned;
    bool skipped;
    ConserveCheck() : first(0), haveFirst(false), comparisons(0), worstDrift(0), warned(false), skipped(false) {}
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
                    if (r.fields[f].name == k.terms[t]) { sum += r.fields[f].conserved(); found = true; break; }
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

            /* Normalise by the LARGER of the two values, not by the first one.
             *
             * Dividing by the first value alone breaks whenever a run starts from an empty
             * domain: the baseline is 0, and the first non-zero total produces a "drift" of
             * 1e+302, which tells the user nothing except that something is wrong somewhere.
             * With max(|first|, |now|) the number is always between 0 and 1 and reads as what
             * it is -- the fraction of the sum that appeared or disappeared. */
            double scale = std::fabs(k.first);
            if (std::fabs(sum) > scale) scale = std::fabs(sum);
            if (scale < 1e-300) scale = 1e-300;
            const double drift = std::fabs(sum - k.first) / scale;
            ++k.comparisons;
            if (drift > k.worstDrift) k.worstDrift = drift;
            if (drift > tol_ && !k.warned) {
                k.warned = true;
                anyFail_ = true;
                char b[1024];
                std::snprintf(b, sizeof(b),
                    "  [DIAG] MASS BALANCE: %s drifted %.3e (tolerance %.1e) by step %ld.\n"
                    "  [DIAG]   started %.10g, now %.10g.\n"
                    "  [DIAG]   Three things cause this, in order of how often they do:\n"
                    "  [DIAG]   1. An OPEN BOUNDARY. A substrate held at a fixed concentration on\n"
                    "  [DIAG]      an inlet is being supplied from outside, so its total is not\n"
                    "  [DIAG]      conserved and never will be. Do not <conserve> it -- conserve a\n"
                    "  [DIAG]      sum that is closed, such as every species sharing one element.\n"
                    "  [DIAG]   2. Reaction stoichiometry that does not balance. Check the sum you\n"
                    "  [DIAG]      named really is a conserved moiety of your reaction network.\n"
                    "  [DIAG]   3. A time step too long for the rate you specified. This one shows\n"
                    "  [DIAG]      up together with negative concentrations, below.\n",
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
            const char *verdict =
                checks_[c].skipped        ? "SKIPPED (unknown substrate)"
              : checks_[c].comparisons==0 ? "SKIPPED (never compared: fewer than two "
                                            "diagnostic rows, so there was no baseline to "
                                            "compare against)"
              : (checks_[c].worstDrift <= tol_ ? "PASS" : "FAIL");
            std::snprintf(b, sizeof(b), "  [DIAG] %-24s worst relative drift %.3e  %s\n",
                          checks_[c].expr.c_str(), checks_[c].worstDrift, verdict);
            s += b;
        }
        /* A verdict is only meaningful if at least one check actually ran. */
        size_t ran = 0;
        for (size_t c = 0; c < checks_.size(); ++c)
            if (!checks_[c].skipped && checks_[c].comparisons > 0) ++ran;
        if (!checks_.empty()) {
            if (ran == 0)
                s += "  [DIAG] verdict: NOT CHECKED -- no conservation check ever ran. Lower "
                     "<interval>, or\n         raise <ade_max_iT>, so the run records at least "
                     "two diagnostic rows.\n";
            else
                s += std::string("  [DIAG] verdict: ") + (anyFail_ ? "FAIL" : "PASS") + "\n";
        }
        return s;
    }

private:
    void writeHeader(const Row &r)
    {
        std::FILE *f = std::fopen(path_.c_str(), "w");
        if (!f) return;
        std::fprintf(f, "# CompLB3D run summary, one row per diagnostic interval\n");
        /* This line used to claim the totals were over open voxels only. They are not, and they
         * must not be: an immobile species keeps its inventory inside SOLID voxels, so excluding
         * those would drop a dissolving mineral out of its own mass balance and the check would
         * fail for a reason that has nothing to do with the chemistry. The totals are over the
         * physical domain; the open-voxel count is reported separately, as porosity, and is what
         * the means are divided by. */
        std::fprintf(f, "# totals are over the physical domain x=1..nx-2; porosity and the means\n");
        std::fprintf(f, "# use the OPEN voxel count, taken from the live mask so it tracks a\n");
        std::fprintf(f, "# pore space that precipitation or dissolution is changing.\n");
        std::fprintf(f, "# <name>_held is the mass sitting in wall and grain voxels, which\n");
        std::fprintf(f, "# computeDensity cannot report because BounceBack and NoDynamics answer\n");
        std::fprintf(f, "# from a stored number. The conserved quantity is _total + _held, and\n");
        std::fprintf(f, "# that is what <conserve> is checked against.\n");
        std::fprintf(f, "iteration,porosity,open_voxels");
        for (size_t i = 0; i < r.fields.size(); ++i) {
            const char *n = r.fields[i].name.c_str();
            std::fprintf(f, ",%s_total,%s_mean,%s_min,%s_max,%s_held", n, n, n, n, n);
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
            std::fprintf(f, ",%.17g,%.17g,%.17g,%.17g,%.17g",
                         r.fields[i].total, r.fields[i].mean, r.fields[i].minv, r.fields[i].maxv,
                         r.fields[i].held);
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
