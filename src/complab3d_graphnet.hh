/* ================================================================================================
 *  complab3d_graphnet.hh  --  A BIPARTITE GRAPH NETWORK OVER A REACTION NETWORK
 *
 *  Part of the CompLaB program.  GNU Affero General Public License v3 or later.
 *  Meile Lab, University of Georgia.  shahram.asgari@uga.edu
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT THIS IS AND WHY IT IS NOT THE .srg SURROGATE
 *
 *  complab3d_surrogate.hh fits one number: the growth rate, from a handful of uptake bounds.  That
 *  is the right shape when the linear program's objective is all you need, and it is why the
 *  surrogate leaves Fout untouched.
 *
 *  This is the other case.  When you need the rate of every species, not just growth, a dense
 *  network has to learn from data which species affect which -- and that is information you already
 *  have.  It is the stoichiometric matrix.  Supplying it as STRUCTURE rather than making the network
 *  infer it means fewer weights, less training data, and a model that cannot invent a coupling the
 *  chemistry does not contain.
 *
 *  ------------------------------------------------------------------------------------------------
 *  THE GRAPH
 *
 *  Bipartite: one node per species, one node per reaction, an edge wherever s[i][r] is non-zero,
 *  carrying that coefficient as its weight.  No approximation -- this is the stoichiometric matrix
 *  read as a graph.
 *
 *  ONE ROUND OF MESSAGE PASSING
 *
 *      reactions   hR_r = tanh( Wsr . SUM_i s[i][r] hS_i  +  Wda da_r  +  bR )
 *      species     hS_i = tanh( Wss . hS_i  +  Wrs . SUM_r s[i][r] hR_r  +  bS )
 *
 *  The sums are stoichiometry-weighted, so a species consuming two moles per reaction speaks twice
 *  as loudly into that reaction as one consuming a single mole.  da_r is the per-reaction Damkohler
 *  number: an edge-level parameter, so one trained network can cover a range of transport regimes
 *  rather than a single one.
 *
 *  ------------------------------------------------------------------------------------------------
 *  THE FILE
 *
 *  Plain text, for the same reason the .srg file is plain text: a trained network is a scientific
 *  result and should be readable, diffable and quotable without a tool.  tools/train_graphnet.py
 *  writes it and this reads it; the two were checked against each other to machine precision, so
 *  the format is not a convention anyone has to maintain by hand.
 *
 *  RANGES ARE ENFORCED, NOT ADVISORY.  Every evaluation is clamped to trainmin..trainmax and the
 *  clamping is counted.  A neural network evaluated outside its training box does not fail; it
 *  returns a confident number.
 *
 *  UNITS.  "per_hour" divides every returned rate by 3600 on load, matching complab3d_symbolic.hh.
 *  A network fitted to flux balance output is per hour; defineKinetics.hh is per second.  Getting
 *  that wrong is a factor of 3600, so the loaded setting is printed in the log.
 * ================================================================================================ */

#ifndef COMPLAB3D_GRAPHNET_HH
#define COMPLAB3D_GRAPHNET_HH

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace complab_gnn {

/* ------------------------------------------------------------------------------------------------
 *  A dense matrix stored row major, so a row is contiguous and the inner loop of every matrix
 *  product walks memory forwards.
 * ---------------------------------------------------------------------------------------------- */
struct Mat {
    int rows, cols;
    std::vector<double> a;
    Mat() : rows(0), cols(0) {}
    void resize(int r, int c) { rows = r; cols = c; a.assign((size_t) r * (size_t) c, 0.0); }
    double  operator()(int r, int c) const { return a[(size_t) r * (size_t) cols + (size_t) c]; }
    double &operator()(int r, int c)       { return a[(size_t) r * (size_t) cols + (size_t) c]; }
    size_t  size() const { return a.size(); }
};

struct Layer { Mat Wsr, Wda, Wss, Wrs; std::vector<double> bR, bS; };

struct Network {
    int nS, nR, width, rounds;
    bool hasGrowth;

    /* ---- HOW THE RATES ARE READ OFF THE NETWORK ------------------------------------------------
     *
     * SPECIES (the original, and what a file with no `readout` line means).  The readout is applied
     * to each SPECIES node and gives that species its rate directly.  The stoichiometric matrix is
     * used as the edge weights of the message passing, so it shapes the answer -- but nothing in
     * the arithmetic forces the returned rates into stoichiometric ratio, and measured on the
     * shipped AOM example they are not: HS produced over CH4 consumed ran between 0.76 and 1.11
     * where the structure implies exactly 1.
     *
     * EXTENT.  The readout is applied to each REACTION node instead and gives one extent per
     * reaction, and the species rates are then formed as
     *
     *     r_i  =  SUM_r  S(i,r) * xi_r
     *
     * which is stoichiometrically exact by construction, at every point of the input space
     * including the ones nobody checked.  It also uses fewer numbers: nR extents rather than nS
     * independently scaled outputs.
     *
     * What it costs is generality.  A rate vector that does NOT lie in the column space of S cannot
     * be represented at all -- which is the point, since such a vector is not a set of reaction
     * rates.  The trainer projects the training data onto that space and reports how much of it did
     * not fit, so a data set that is inconsistent with the declared stoichiometry says so rather
     * than being approximated. */
    enum Readout { SPECIES = 0, EXTENT = 1 };
    int readout;

    /* EXTENT mode only: one scale per reaction, turning the network's O(1) output into an extent in
     * the file's units.  A pure scale with no offset, because an offset does not commute with the
     * sum above and would break the exactness this mode exists for. */
    std::vector<double> xiScale;

    std::vector<std::string> species, reactions;
    Mat  S;                              /* nS x nR stoichiometry, also the edge weights */
    std::vector<double> da;              /* per reaction */

    std::vector<double> xOffset, xGain;  /* per species */
    double xYmin;
    std::vector<double> yOffset, yGain;  /* per output: nS species rates, then growth */
    double yYmin;

    std::vector<double> trainMin, trainMax;

    Mat Wenc; std::vector<double> bEnc;
    std::vector<Layer> layer;
    Mat Wout; std::vector<double> bOut;

    double unitScale;
    std::string unitName, provenance;

    Network() : nS(0), nR(0), width(0), rounds(0), hasGrowth(false), readout(SPECIES),
                xYmin(-1.0), yYmin(-1.0), unitScale(1.0), unitName("per_second") {}

    /* How many numbers come back from eval(): one per species, plus growth if the file says so. */
    int nOut() const { return hasGrowth ? nS + 1 : nS; }

    /* How many ROWS the readout matrix has, which is NOT the same thing.  Row 0 is shared by every
     * species node -- the species are told apart by their own hidden state and by the per-species
     * yoffset/ygain, not by separate readout weights -- and row 1, when present, reads growth off
     * the mean node.  Giving each species its own row would be nS-1 extra rows the trainer never
     * touches, and would invite a reader to assume row i belongs to species i.  It does not. */
    int nHead() const { return hasGrowth ? 2 : 1; }

    /* How many yoffset/ygain entries the file must carry.  In SPECIES mode the readout produces one
     * independently scaled number per species, plus growth.  In EXTENT mode the species rates are
     * not scaled at all -- xiScale does that job, per reaction -- so the only entry left is growth,
     * and a file with no growth output needs none. */
    int nYScale() const { return (readout == EXTENT) ? (hasGrowth ? 1 : 0) : nOut(); }

    bool valid() const {
        if (nS <= 0 || nR <= 0 || width <= 0 || rounds <= 0) return false;
        if ((int) layer.size() != rounds) return false;
        if (S.rows != nS || S.cols != nR) return false;
        if ((int) da.size() != nR) return false;
        if (Wenc.rows != width || Wenc.cols != 1) return false;
        if ((int) bEnc.size() != width) return false;
        if (Wout.rows != nHead() || Wout.cols != width) return false;
        if ((int) bOut.size() != nHead()) return false;
        if ((int) xOffset.size() != nS || (int) xGain.size() != nS) return false;
        if ((int) yOffset.size() != nYScale() || (int) yGain.size() != nYScale()) return false;
        if (readout == EXTENT && (int) xiScale.size() != nR) return false;
        if ((int) trainMin.size() != nS || (int) trainMax.size() != nS) return false;
        /* Every layer array must be fully present.  A keyword that never appeared leaves its
         * vector empty, and an empty bS[] would be read past the end in the first round. */
        for (size_t l = 0; l < layer.size(); ++l) {
            const Layer &L = layer[l];
            if (L.Wsr.rows != width || L.Wsr.cols != width) return false;
            if (L.Wss.rows != width || L.Wss.cols != width) return false;
            if (L.Wrs.rows != width || L.Wrs.cols != width) return false;
            if (L.Wda.rows != width || L.Wda.cols != 1)     return false;
            if ((int) L.bR.size() != width) return false;
            if ((int) L.bS.size() != width) return false;
        }
        return true;
    }

    bool inTrainingRange(const std::vector<double> &c) const {
        for (int i = 0; i < nS && i < (int) c.size(); ++i)
            if (c[(size_t) i] < trainMin[(size_t) i] || c[(size_t) i] > trainMax[(size_t) i])
                return false;
        return true;
    }

    /* Working room for one forward pass.  The caller owns it and hands the same one back every
     * time, so the five buffers are allocated once for a sweep rather than once per voxel.
     *
     * Be honest about what this is worth: on the small AOM example it is about 3% (2.55 to 2.47
     * microseconds per voxel), because that network is dominated by its 84 tanh calls, which are
     * roughly 40% of the total on their own.  It matters more as the network grows, since the
     * buffers grow with it while the allocation count does not.  The reason to keep it is that it
     * costs nothing and removes a per-voxel allocation from an inner loop; it is not a speed-up
     * anyone should quote. */
    struct Scratch { std::vector<double> hS, mR, hR, mS, hNew; };

    /* --------------------------------------------------------------------------------------------
     *  Forward pass.  conc has one entry per species; out comes back with one rate per species and,
     *  when the file says so, a growth rate in the last slot.
     *
     *  The short form allocates its own scratch and is meant for tests and one-off calls.  Anything
     *  in an inner loop should hold a Scratch and use the other one.
     * ------------------------------------------------------------------------------------------ */
    void eval(const std::vector<double> &conc, std::vector<double> &out,
              bool clamp = true, long *clampCount = 0) const
    {
        Scratch tmp;
        eval(conc, out, tmp, clamp, clampCount);
    }

    void eval(const std::vector<double> &conc, std::vector<double> &out, Scratch &sc,
              bool clamp = true, long *clampCount = 0) const
    {
        out.assign((size_t) nOut(), 0.0);
        if (!valid()) return;

        const size_t W = (size_t) width;
        std::vector<double> &hS = sc.hS;
        hS.assign((size_t) nS * W, 0.0);

        /* encode: each species node starts from its own scaled concentration.
         * [FIX] One count per EVALUATION that clamped anything, not one per species:
         * runtimeReport() divides this by the evaluation count, so counting species made a
         * four-species network able to report 400% of its evaluations as out of range. */
        bool anyClamped = false;
        for (int i = 0; i < nS; ++i) {
            double c = (i < (int) conc.size()) ? conc[(size_t) i] : 0.0;
            if (clamp) {
                if (c < trainMin[(size_t) i])      { c = trainMin[(size_t) i]; anyClamped = true; }
                else if (c > trainMax[(size_t) i]) { c = trainMax[(size_t) i]; anyClamped = true; }
            }
            const double x = (c - xOffset[(size_t) i]) * xGain[(size_t) i] + xYmin;
            for (size_t h = 0; h < W; ++h)
                hS[(size_t) i * W + h] = std::tanh(x * Wenc((int) h, 0) + bEnc[h]);
        }
        if (anyClamped && clampCount) ++*clampCount;

        std::vector<double> &mR = sc.mR; std::vector<double> &hR = sc.hR;
        std::vector<double> &mS = sc.mS; std::vector<double> &hNew = sc.hNew;
        mR.resize((size_t) nR * W); hR.resize((size_t) nR * W);
        mS.resize((size_t) nS * W); hNew.resize((size_t) nS * W);

        for (int l = 0; l < rounds; ++l) {
            const Layer &L = layer[(size_t) l];

            /* species -> reactions, weighted by the stoichiometry */
            std::fill(mR.begin(), mR.end(), 0.0);
            for (int i = 0; i < nS; ++i)
                for (int r = 0; r < nR; ++r) {
                    const double s = S(i, r);
                    if (s == 0.0) continue;                 /* no edge */
                    for (size_t h = 0; h < W; ++h)
                        mR[(size_t) r * W + h] += s * hS[(size_t) i * W + h];
                }
            for (int r = 0; r < nR; ++r)
                for (size_t h = 0; h < W; ++h) {
                    double acc = L.bR[h] + da[(size_t) r] * L.Wda((int) h, 0);
                    for (size_t k = 0; k < W; ++k)
                        acc += L.Wsr((int) h, (int) k) * mR[(size_t) r * W + k];
                    hR[(size_t) r * W + h] = std::tanh(acc);
                }

            /* reactions -> species, same weights */
            std::fill(mS.begin(), mS.end(), 0.0);
            for (int i = 0; i < nS; ++i)
                for (int r = 0; r < nR; ++r) {
                    const double s = S(i, r);
                    if (s == 0.0) continue;
                    for (size_t h = 0; h < W; ++h)
                        mS[(size_t) i * W + h] += s * hR[(size_t) r * W + h];
                }
            for (int i = 0; i < nS; ++i)
                for (size_t h = 0; h < W; ++h) {
                    double acc = L.bS[h];
                    for (size_t k = 0; k < W; ++k)
                        acc += L.Wss((int) h, (int) k) * hS[(size_t) i * W + k]
                             + L.Wrs((int) h, (int) k) * mS[(size_t) i * W + k];
                    hNew[(size_t) i * W + h] = std::tanh(acc);
                }
            hS.swap(hNew);
        }

        if (readout == EXTENT) {
            /* One extent per reaction, off that reaction's node, then the species rates as
             * r_i = SUM_r S(i,r) xi_r.  hR still holds the last round's reaction states.
             *
             * The ratio between two species is therefore whatever S says it is, at every input,
             * and not something the fit has to be trusted to have got approximately right. */
            std::vector<double> &xi = sc.mR;         /* mR is finished with; reuse it */
            for (int r = 0; r < nR; ++r) {
                double y = bOut[0];
                for (size_t h = 0; h < W; ++h) y += Wout(0, (int) h) * hR[(size_t) r * W + h];
                if (y != y) y = 0.0;
                xi[(size_t) r] = y * xiScale[(size_t) r];
            }
            for (int i = 0; i < nS; ++i) {
                double v = 0.0;
                for (int r = 0; r < nR; ++r) v += S(i, r) * xi[(size_t) r];
                out[(size_t) i] = v * unitScale;
            }
        } else {
            /* read out: row 0 of Wout gives every species its rate from its own node */
            for (int i = 0; i < nS; ++i) {
                double y = bOut[0];
                for (size_t h = 0; h < W; ++h) y += Wout(0, (int) h) * hS[(size_t) i * W + h];
                double v = (y - yYmin) / yGain[(size_t) i] + yOffset[(size_t) i];
                if (v != v) v = 0.0;
                out[(size_t) i] = v * unitScale;
            }
        }
        if (hasGrowth) {
            /* growth reads the mean species node, so it sees the whole network rather than one
             * species; this is what the trainer fits.  It is not a species rate and no
             * stoichiometry constrains it, so it keeps its own offset and gain in both modes --
             * the last entry in SPECIES mode, the only one in EXTENT mode. */
            const size_t gs = (readout == EXTENT) ? 0 : (size_t) nS;
            double y = bOut[1];
            for (size_t h = 0; h < W; ++h) {
                double m = 0.0;
                for (int i = 0; i < nS; ++i) m += hS[(size_t) i * W + h];
                y += Wout(1, (int) h) * (m / (double) nS);
            }
            double v = (y - yYmin) / yGain[gs] + yOffset[gs];
            if (v != v) v = 0.0;
            out[(size_t) nS] = v * unitScale;
        }
    }
};

/* ================================================================================================
 *  Runtime state, held the way complab_srg::runtime() holds the surrogate's
 * ================================================================================================ */
struct Binding {
    const Network *net;
    std::vector<int> subsOfSpecies;   /* species slot -> substrate index */
    int growthSlot;                   /* index into eval() output, or -1 */
    Binding() : net(0), growthSlot(-1) {}
};
struct Runtime {
    std::vector<Binding> byMicrobe;

    /* The abiotic network, if there is one.  Held apart from byMicrobe because it belongs to no
     * organism and is swept over every fluid voxel.  net == 0 means none was loaded. */
    Binding abiotic;

    long evaluations, clamped;
    Runtime() : evaluations(0), clamped(0) {}
};
inline Runtime &runtime() { static Runtime R; return R; }

/* ================================================================================================
 *  Reading the file
 * ================================================================================================ */
namespace detail {
    inline bool readN(std::FILE *f, std::vector<double> &v, size_t n) {
        v.assign(n, 0.0);
        for (size_t i = 0; i < n; ++i) if (std::fscanf(f, "%lf", &v[i]) != 1) return false;
        return true;
    }
    inline bool readMat(std::FILE *f, Mat &m, int r, int c) {
        m.resize(r, c);
        for (size_t i = 0; i < m.size(); ++i) if (std::fscanf(f, "%lf", &m.a[i]) != 1) return false;
        return true;
    }
    inline void words(const std::string &s, std::vector<std::string> &out) {
        out.clear(); size_t i = 0;
        while (i < s.size()) {
            while (i < s.size() && std::isspace((unsigned char) s[i])) ++i;
            size_t st = i;
            while (i < s.size() && !std::isspace((unsigned char) s[i])) ++i;
            if (i > st) out.push_back(s.substr(st, i - st));
        }
    }
}

inline bool load(Network &N, const std::string &path, std::string *err = 0)
{
    std::FILE *f = std::fopen(path.c_str(), "r");
    if (!f) { if (err) *err = "cannot open '" + path + "'"; return false; }
    N = Network();

    char key[128];
    int curLayer = -1;
    bool sawShape = false;

    while (std::fscanf(f, "%127s", key) == 1) {
        std::string k(key);

        if (k.size() && k[0] == '#') {                       /* comment: skipped to end of line */
            char line[2048];
            if (!std::fgets(line, sizeof line, f)) break;
            continue;
        }
        if (k == "provenance") {                             /* where this network came from */
            char line[2048];
            if (!std::fgets(line, sizeof line, f)) break;
            std::string c(line);
            while (!c.empty() && (c[c.size()-1]=='\n' || c[c.size()-1]=='\r')) c.erase(c.size()-1);
            size_t b = c.find_first_not_of(" \t");
            if (b != std::string::npos) {
                c = c.substr(b);
                if (!N.provenance.empty()) N.provenance += " | ";
                N.provenance += c;
            }
            continue;
        }
        if (k == "version") { int v; if (std::fscanf(f, "%d", &v) != 1) break; continue; }
        if (k == "units") {
            if (std::fscanf(f, "%127s", key) != 1) break;
            std::string u(key);
            for (size_t i = 0; i < u.size(); ++i) u[i] = (char) std::tolower((unsigned char) u[i]);
            if      (u == "per_second") { N.unitScale = 1.0;          N.unitName = "per_second"; }
            else if (u == "per_hour")   { N.unitScale = 1.0 / 3600.0; N.unitName = "per_hour";   }
            else { std::fclose(f); if (err) *err = "'units' must be per_second or per_hour"; return false; }
            continue;
        }
        if (k == "species" || k == "reactions") {
            char line[4096];
            if (!std::fgets(line, sizeof line, f)) break;
            std::vector<std::string> w; detail::words(line, w);
            if (k == "species") N.species = w; else N.reactions = w;
            continue;
        }
        if (k == "rounds") { if (std::fscanf(f, "%d", &N.rounds) != 1) break; continue; }
        if (k == "width")  { if (std::fscanf(f, "%d", &N.width)  != 1) break; continue; }
        if (k == "growth") { int g; if (std::fscanf(f, "%d", &g) != 1) break; N.hasGrowth = (g != 0); continue; }
        if (k == "readout") {
            if (std::fscanf(f, "%127s", key) != 1) break;
            std::string m(key);
            for (size_t i = 0; i < m.size(); ++i) m[i] = (char) std::tolower((unsigned char) m[i]);
            if      (m == "species") N.readout = Network::SPECIES;
            else if (m == "extent")  N.readout = Network::EXTENT;
            else {
                std::fclose(f);
                if (err) *err = "'readout' must be species or extent, not '" + m + "'";
                return false;
            }
            continue;
        }
        if (k == "stoich") {
            if (std::fscanf(f, "%d %d", &N.nS, &N.nR) != 2) break;
            if (!detail::readMat(f, N.S, N.nS, N.nR)) break;
            sawShape = true;
            continue;
        }
        if (!sawShape) {
            std::fclose(f);
            if (err) *err = "'" + path + "': '" + k + "' appears before the stoich block";
            return false;
        }

        if      (k == "da")       { if (!detail::readN(f, N.da,       (size_t) N.nR)) break; }
        else if (k == "xoffset")  { if (!detail::readN(f, N.xOffset,  (size_t) N.nS)) break; }
        else if (k == "xgain")    { if (!detail::readN(f, N.xGain,    (size_t) N.nS)) break; }
        else if (k == "xymin")    { if (std::fscanf(f, "%lf", &N.xYmin) != 1) break; }
        /* nYScale() is nOut() in SPECIES mode and just the growth entry in EXTENT mode, so a file
         * written for one mode and labelled the other fails here rather than reading the wrong
         * count of numbers and shifting every key after it. */
        else if (k == "yoffset")  { if (!detail::readN(f, N.yOffset,  (size_t) N.nYScale())) break; }
        else if (k == "ygain")    { if (!detail::readN(f, N.yGain,    (size_t) N.nYScale())) break; }
        else if (k == "xiscale")  { if (!detail::readN(f, N.xiScale,  (size_t) N.nR)) break; }
        else if (k == "yymin")    { if (std::fscanf(f, "%lf", &N.yYmin) != 1) break; }
        else if (k == "trainmin") { if (!detail::readN(f, N.trainMin, (size_t) N.nS)) break; }
        else if (k == "trainmax") { if (!detail::readN(f, N.trainMax, (size_t) N.nS)) break; }
        else if (k == "Wenc")     { if (!detail::readMat(f, N.Wenc, N.width, 1)) break; }
        else if (k == "benc")     { if (!detail::readN(f, N.bEnc, (size_t) N.width)) break; }
        else if (k == "layer") {
            if (std::fscanf(f, "%d", &curLayer) != 1) break;
            if ((int) N.layer.size() <= curLayer) N.layer.resize((size_t) curLayer + 1);
        }
        else if (k == "Wsr" || k == "Wss" || k == "Wrs" || k == "Wda" || k == "bR" || k == "bS") {
            if (curLayer < 0 || curLayer >= (int) N.layer.size()) {
                std::fclose(f);
                if (err) *err = "'" + path + "': '" + k + "' outside any layer block";
                return false;
            }
            Layer &L = N.layer[(size_t) curLayer];
            bool ok = true;
            if      (k == "Wsr") ok = detail::readMat(f, L.Wsr, N.width, N.width);
            else if (k == "Wss") ok = detail::readMat(f, L.Wss, N.width, N.width);
            else if (k == "Wrs") ok = detail::readMat(f, L.Wrs, N.width, N.width);
            else if (k == "Wda") ok = detail::readMat(f, L.Wda, N.width, 1);
            else if (k == "bR")  ok = detail::readN(f, L.bR, (size_t) N.width);
            else if (k == "bS")  ok = detail::readN(f, L.bS, (size_t) N.width);
            if (!ok) break;
        }
        else if (k == "Wout") { if (!detail::readMat(f, N.Wout, N.nHead(), N.width)) break; }
        else if (k == "bout") { if (!detail::readN(f, N.bOut, (size_t) N.nHead())) break; }
        else { std::fclose(f); if (err) *err = "'" + path + "': unknown keyword '" + k + "'"; return false; }
    }
    std::fclose(f);

    if (!N.valid()) {
        if (err) *err = "'" + path + "' is incomplete or inconsistent: it should declare stoich, "
                        "rounds, width, one layer block per round, Wenc/benc and Wout/bout";
        return false;
    }
    if (N.trainMin.size() != (size_t) N.nS || N.trainMax.size() != (size_t) N.nS) {
        if (err) *err = "'" + path + "' has no trainmin/trainmax, so its valid range is unknown";
        return false;
    }
    return true;
}

inline std::string describe(const Network &N, const std::string &path)
{
    char b[512];
    std::string s = "  [GNN] graph network from '" + path + "'\n";
    if (!N.provenance.empty()) s += "  [GNN]   provenance: " + N.provenance + "\n";
    std::sprintf(b, "  [GNN]   %d species, %d reactions, width %d, %d message-passing round%s\n",
                 N.nS, N.nR, N.width, N.rounds, N.rounds == 1 ? "" : "s");
    s += b;
    s += "  [GNN]   units: " + N.unitName;
    if (N.unitScale != 1.0) s += "  (every rate divided by 3600 on load)";
    s += "\n";
    s += "  [GNN]   outputs: a rate per species";
    s += N.hasGrowth ? ", plus growth\n" : "\n";
    if (N.readout == Network::EXTENT)
        s += "  [GNN]   readout: EXTENT -- one extent per reaction, species rates formed as\n"
             "  [GNN]            r_i = SUM_r S(i,r) xi_r, so the stoichiometric ratio is exact\n"
             "  [GNN]            at every input rather than approximately fitted.\n";
    else
        s += "  [GNN]   readout: SPECIES -- each species reads its own rate off its own node.\n"
             "  [GNN]            The stoichiometry shapes the message passing but does NOT\n"
             "  [GNN]            constrain the answer, so the returned rates are only\n"
             "  [GNN]            approximately in ratio. Retrain with --readout extent to make\n"
             "  [GNN]            them exact.\n";
    s += "  [GNN]   valid over:";
    for (int i = 0; i < N.nS; ++i) {
        std::sprintf(b, " %s %.6g..%.6g",
                     (i < (int) N.species.size() ? N.species[(size_t) i].c_str() : "?"),
                     N.trainMin[(size_t) i], N.trainMax[(size_t) i]);
        s += b;
    }
    s += "\n";
    return s;
}

/* ================================================================================================
 *  BINDING A NETWORK TO THIS RUN'S LATTICES
 *
 *  The .gnn file names the species it was fitted to.  This matches those names against
 *  <name_of_substrates> and records the mapping, so a network fitted to CH4 SO4 HS HCO3 cannot be
 *  applied to a simulation that calls them something else.  Returns "" on success, or the message
 *  to print before terminating.
 * ================================================================================================ */
inline std::string bindToSubstrates(const Network &N,
                                    const std::vector<std::string> &subsNames,
                                    Binding &out,
                                    bool abiotic = false)
{
    out = Binding();
    out.net = &N;

    if (!N.valid()) return "the network file did not load into a usable network.";
    if ((int) N.species.size() != N.nS) {
        char b[192];
        std::sprintf(b, "the file names %d species but declares %d.",
                     (int) N.species.size(), N.nS);
        return std::string(b);
    }

    out.subsOfSpecies.assign((size_t) N.nS, -1);
    for (int i = 0; i < N.nS; ++i) {
        int hit = -1;
        for (size_t s = 0; s < subsNames.size(); ++s)
            if (subsNames[s] == N.species[(size_t) i]) { hit = (int) s; break; }
        if (hit < 0) {
            std::string m = "the network was trained on species '" + N.species[(size_t) i]
                          + "', which is not in <name_of_substrates>. This run has:";
            for (size_t s = 0; s < subsNames.size(); ++s) m += " " + subsNames[s];
            m += ".";
            return m;
        }
        out.subsOfSpecies[(size_t) i] = hit;
    }

    /* eval() returns one value per species, then growth if the file declares it. */
    if (N.hasGrowth) {
        if (abiotic)
            return "the abiotic network declares a growth output. An abiotic network belongs to no "
                   "organism, so there is no biomass for it to apply to. Retrain it without the "
                   "growth column, or give it to a microbe through <network_file>.";
        out.growthSlot = N.nS;
    } else {
        out.growthSlot = -1;
    }
    return std::string();
}

inline void registerNetwork(int microbe, const Binding &b, int nMicrobes)
{
    Runtime &R = runtime();
    if ((int) R.byMicrobe.size() < nMicrobes) R.byMicrobe.resize((size_t) nMicrobes);
    if (microbe < 0) {
        for (size_t i = 0; i < R.byMicrobe.size(); ++i) R.byMicrobe[i] = b;
    } else if (microbe < (int) R.byMicrobe.size()) {
        R.byMicrobe[(size_t) microbe] = b;
    }
}

inline void registerAbiotic(const Binding &b) { runtime().abiotic = b; }

inline bool haveAbiotic()
{
    const Binding &b = runtime().abiotic;
    return b.net != 0 && b.net->valid();
}

inline std::string runtimeReport()
{
    Runtime &R = runtime();
    if (R.evaluations == 0) return std::string();
    char b[512];
    const double pct = 100.0 * (double) R.clamped / (double) R.evaluations;
    std::sprintf(b, "  [GNN] %ld network evaluations, %ld inputs clamped to the training box (%.2f%%)\n",
                 R.evaluations, R.clamped, pct);
    std::string s(b);
    if (pct > 5.0)
        s += "  [GNN] WARNING: more than 5% of inputs were outside the box this network was\n"
             "  [GNN] trained on. The results should not be used until the ranges are widened\n"
             "  [GNN] and the network refitted.\n";
    return s;
}

}  // namespace complab_gnn

#endif  // COMPLAB3D_GRAPHNET_HH
