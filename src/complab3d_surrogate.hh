/* This file is a part of the CompLaB program.
 *
 * The CompLaB3D software is developed since 2024 by the University of Georgia
 * (Meile Lab, Department of Marine Sciences).  Contact: Shahram Asgari,
 * Christof Meile -- shahram.asgari@uga.edu
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* ================================================================================================
 * complab3d_surrogate.hh  --  TRAIN AND USE A SURROGATE NETWORK, ENTIRELY FROM THE XML
 * ================================================================================================
 *
 *  WHAT THIS REPLACES
 *
 *  Fitting a surrogate used to take six manual steps and a rebuild:
 *      generateTrainingData.py -> trainSurrogate.py or .m -> exportSurrogateHeader.m
 *      -> verifyExport.py -> paste an #include and a branch into surrogateModel.hh -> rebuild
 *
 *  All of it is here, in C++, driven from one XML block:
 *
 *      <surrogate>
 *          <enabled>true</enabled>
 *          <weights_file>geobacter_net.srg</weights_file>
 *          <train_if_missing>true</train_if_missing>
 *          <train>
 *              <samples>4000</samples>
 *              <input0>acetate 1e-3 10  log</input0>
 *              <input1>Fe3     1e-5 0.5 log</input1>
 *              <layers>10 10 10 10</layers>
 *              <restarts>5</restarts>
 *              <verify_points>512</verify_points>
 *          </train>
 *      </surrogate>
 *
 *  THE CHANGE THAT MAKES THIS POSSIBLE is that the weights now live in a FILE the solver reads at
 *  run time, instead of being compiled into surrogateModel.hh. That single change removes the
 *  rebuild, and it lets the training range travel with the weights -- so the solver can refuse to
 *  extrapolate instead of returning confident nonsense, which is the failure mode that makes
 *  surrogates dangerous.
 *
 *  ------------------------------------------------------------------------------------------------
 *  HOW THE TRAINING SWEEP GETS ITS ANSWERS
 *
 *  It calls back into whatever flux balance solver the caller supplies:
 *
 *      typedef double (*FbaFn)(const std::vector<double>& uptake, void* ctx);
 *
 *  complab.cpp wires that to the GLPK path that is already linked. Passing it in rather than
 *  calling GLPK directly here keeps this header testable on its own -- the test suite fits it to
 *  an analytic function with no LP solver in sight -- and means the same trainer would work
 *  against COBRApy or any future solver without change.
 *
 *  ------------------------------------------------------------------------------------------------
 *  THE NETWORK
 *
 *  Identical in shape and arithmetic to the one CompLaB has shipped since v1.0: mapminmax scaling
 *  on both ends, tansig hidden layers, one linear output. A .srg file and the old compiled-in
 *  weights therefore describe the same function, and inspectSurrogate.py reads both.
 *
 *  Training is Adam with restarts and early stopping on a held-out split. That is not what MATLAB's
 *  trainlm does, and the weights will differ; what must agree is the FIT, and verifyAgainstFba()
 *  below measures exactly that by re-solving the LP at points the network never saw.
 *
 *  ------------------------------------------------------------------------------------------------
 *  KNOWN LIMITS
 *
 *    * Training needs a flux balance solver, so ENABLE_GLPK (or COBRApy) must be on. Using an
 *      already-trained .srg file needs neither.
 *    * 1 to 3 inputs. Past that the sample count needed grows faster than the fit improves.
 *    * Training happens on rank 0 and the result is broadcast; it is not itself parallel. On a
 *      genome-scale model a 4000-point sweep is a few minutes, once.
 * ================================================================================================
 */
#ifndef COMPLAB3D_SURROGATE_HH
#define COMPLAB3D_SURROGATE_HH

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace complab_srg {

/* ================================================================================================
 *  A deterministic RNG.
 *
 *  Written out rather than using <random> so that a given seed gives the same network on every
 *  compiler and every machine. std::mt19937 is portable but the DISTRIBUTIONS are not: the same
 *  seed through std::uniform_real_distribution gives different numbers on libstdc++ and libc++,
 *  which would make a trained network irreproducible across platforms.
 * ================================================================================================ */
struct Rng {
    unsigned long long s;
    explicit Rng(unsigned long long seed = 1) : s(seed ? seed : 1) {}
    unsigned long long next() {                 // xorshift64*
        s ^= s >> 12; s ^= s << 25; s ^= s >> 27;
        return s * 2685821657736338717ULL;
    }
    double uniform() { return (double) (next() >> 11) / 9007199254740992.0; }
    double uniform(double a, double b) { return a + (b - a) * uniform(); }
    double normal() {                           // Box-Muller
        double u1 = uniform(), u2 = uniform();
        if (u1 < 1e-300) u1 = 1e-300;
        return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    }
};

/* ================================================================================================
 *  The network
 * ================================================================================================ */
struct Network {
    int nIn;
    std::vector<int> sizes;                 // nIn, h1..hn, 1
    std::vector<std::vector<double> > W;    // per layer, row major (rows x cols)
    std::vector<std::vector<double> > B;

    std::vector<double> xOffset, xGain;     // mapminmax on the inputs
    double xYmin;
    double yOffset, yGain, yYmin;
    bool logOutput;

    std::vector<double> trainMin, trainMax; // the box the data covered
    double outMin, outMax;

    std::string provenance;                 // model, sample count, date: written into the file

    Network() : nIn(0), xYmin(-1.0), yOffset(0.0), yGain(1.0), yYmin(-1.0),
                logOutput(false), outMin(0.0), outMax(0.0) {}

    bool valid() const { return nIn > 0 && !W.empty() && W.size() == B.size(); }

    static double tansig(double n) { return 2.0 / (1.0 + std::exp(-2.0 * n)) - 1.0; }

    /* True when every input is inside the box the network was trained on. Outside it the return
     * value of eval() is an extrapolation and means nothing. */
    bool inTrainingRange(const std::vector<double> &x) const {
        if ((int) x.size() < nIn) return false;
        for (int i = 0; i < nIn; ++i)
            if (x[i] < trainMin[i] || x[i] > trainMax[i]) return false;
        return true;
    }

    double eval(const std::vector<double> &x) const {
        if (!valid() || (int) x.size() < nIn) return 0.0;
        std::vector<double> a((size_t) nIn);
        for (int i = 0; i < nIn; ++i) a[(size_t) i] = (x[(size_t) i] - xOffset[(size_t) i]) * xGain[(size_t) i] + xYmin;

        for (size_t L = 0; L < W.size(); ++L) {
            const int rows = sizes[L + 1], cols = sizes[L];
            std::vector<double> out((size_t) rows);
            for (int r = 0; r < rows; ++r) {
                double s = B[L][(size_t) r];
                const double *w = &W[L][(size_t) r * (size_t) cols];
                for (int c = 0; c < cols; ++c) s += w[c] * a[(size_t) c];
                out[(size_t) r] = (L + 1 == W.size()) ? s : tansig(s);
            }
            a.swap(out);
        }
        double g = (a[0] - yYmin) / yGain + yOffset;
        if (logOutput) g = std::pow(10.0, g);
        if (!(g > 1e-8)) g = 0.0;               // also catches NaN
        return g;
    }
};

/* ================================================================================================
 *  The weights file
 *
 *  Plain text on purpose. A trained network is a scientific result: it should be readable, diffable
 *  and quotable in a paper without a tool. It is small -- a few hundred numbers -- so nothing is
 *  gained by making it binary. Written with %.17g, which round-trips an IEEE double exactly.
 * ================================================================================================ */
inline bool save(const Network &N, const std::string &path, std::string *err = 0)
{
    std::FILE *f = std::fopen(path.c_str(), "w");
    if (!f) { if (err) *err = "cannot write '" + path + "'"; return false; }

    std::fprintf(f, "# CompLB3D surrogate network\n");
    std::fprintf(f, "# %s\n", N.provenance.c_str());
    std::fprintf(f, "# THE NETWORK IS ONLY VALID INSIDE trainMin..trainMax. Outside that box a\n");
    std::fprintf(f, "# neural network does not fail, it returns confident nonsense.\n");
    std::fprintf(f, "version 1\n");
    std::fprintf(f, "inputs %d\n", N.nIn);
    std::fprintf(f, "layers");
    for (size_t i = 0; i < N.sizes.size(); ++i) std::fprintf(f, " %d", N.sizes[i]);
    std::fprintf(f, "\n");
    std::fprintf(f, "logoutput %d\n", N.logOutput ? 1 : 0);

    std::fprintf(f, "xoffset"); for (int i = 0; i < N.nIn; ++i) std::fprintf(f, " %.17g", N.xOffset[(size_t) i]); std::fprintf(f, "\n");
    std::fprintf(f, "xgain");   for (int i = 0; i < N.nIn; ++i) std::fprintf(f, " %.17g", N.xGain[(size_t) i]);   std::fprintf(f, "\n");
    std::fprintf(f, "xymin %.17g\n", N.xYmin);
    std::fprintf(f, "yoffset %.17g\nygain %.17g\nyymin %.17g\n", N.yOffset, N.yGain, N.yYmin);
    std::fprintf(f, "trainmin"); for (int i = 0; i < N.nIn; ++i) std::fprintf(f, " %.17g", N.trainMin[(size_t) i]); std::fprintf(f, "\n");
    std::fprintf(f, "trainmax"); for (int i = 0; i < N.nIn; ++i) std::fprintf(f, " %.17g", N.trainMax[(size_t) i]); std::fprintf(f, "\n");
    std::fprintf(f, "outrange %.17g %.17g\n", N.outMin, N.outMax);

    for (size_t L = 0; L < N.W.size(); ++L) {
        std::fprintf(f, "W %d\n", (int) L);
        const int rows = N.sizes[L + 1], cols = N.sizes[L];
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) std::fprintf(f, "%.17g ", N.W[L][(size_t) r * (size_t) cols + (size_t) c]);
            std::fprintf(f, "\n");
        }
        std::fprintf(f, "B %d\n", (int) L);
        for (int r = 0; r < rows; ++r) std::fprintf(f, "%.17g ", N.B[L][(size_t) r]);
        std::fprintf(f, "\n");
    }
    std::fclose(f);
    return true;
}

inline bool load(Network &N, const std::string &path, std::string *err = 0)
{
    std::FILE *f = std::fopen(path.c_str(), "r");
    if (!f) { if (err) *err = "cannot open '" + path + "'"; return false; }

    char key[64];
    N = Network();
    bool haveLayers = false;
    while (std::fscanf(f, "%63s", key) == 1) {
        if (key[0] == '#') { int ch; while ((ch = std::fgetc(f)) != EOF && ch != '\n') {} continue; }
        std::string k(key);
        if (k == "version") { int v; if (std::fscanf(f, "%d", &v) != 1) break; }
        else if (k == "inputs") { if (std::fscanf(f, "%d", &N.nIn) != 1) break; }
        else if (k == "layers") {
            /* read integers until the line ends */
            N.sizes.clear();
            int ch, v;
            while (true) {
                ch = std::fgetc(f);
                if (ch == '\n' || ch == EOF) break;
                if (ch == ' ' || ch == '\t' || ch == '\r') continue;
                std::ungetc(ch, f);
                if (std::fscanf(f, "%d", &v) != 1) break;
                N.sizes.push_back(v);
            }
            haveLayers = true;
        }
        else if (k == "logoutput") { int v = 0; if (std::fscanf(f, "%d", &v) != 1) break; N.logOutput = (v != 0); }
        else if (k == "xoffset") { N.xOffset.assign((size_t) N.nIn, 0.0); for (int i = 0; i < N.nIn; ++i) if (std::fscanf(f, "%lf", &N.xOffset[(size_t) i]) != 1) break; }
        else if (k == "xgain")   { N.xGain.assign((size_t) N.nIn, 1.0);   for (int i = 0; i < N.nIn; ++i) if (std::fscanf(f, "%lf", &N.xGain[(size_t) i]) != 1) break; }
        else if (k == "xymin")   { if (std::fscanf(f, "%lf", &N.xYmin) != 1) break; }
        else if (k == "yoffset") { if (std::fscanf(f, "%lf", &N.yOffset) != 1) break; }
        else if (k == "ygain")   { if (std::fscanf(f, "%lf", &N.yGain) != 1) break; }
        else if (k == "yymin")   { if (std::fscanf(f, "%lf", &N.yYmin) != 1) break; }
        else if (k == "trainmin"){ N.trainMin.assign((size_t) N.nIn, 0.0); for (int i = 0; i < N.nIn; ++i) if (std::fscanf(f, "%lf", &N.trainMin[(size_t) i]) != 1) break; }
        else if (k == "trainmax"){ N.trainMax.assign((size_t) N.nIn, 0.0); for (int i = 0; i < N.nIn; ++i) if (std::fscanf(f, "%lf", &N.trainMax[(size_t) i]) != 1) break; }
        else if (k == "outrange"){ if (std::fscanf(f, "%lf %lf", &N.outMin, &N.outMax) != 2) break; }
        else if (k == "W" || k == "B") {
            int L = 0;
            if (std::fscanf(f, "%d", &L) != 1) break;
            if (!haveLayers || L + 1 >= (int) N.sizes.size()) break;
            const int rows = N.sizes[(size_t) L + 1], cols = N.sizes[(size_t) L];
            if ((int) N.W.size() < L + 1) { N.W.resize((size_t) L + 1); N.B.resize((size_t) L + 1); }
            const int n = (k == "W") ? rows * cols : rows;
            std::vector<double> v((size_t) n, 0.0);
            bool truncated = false;
            for (int i = 0; i < n && !truncated; ++i)
                if (std::fscanf(f, "%lf", &v[(size_t) i]) != 1) truncated = true;
            if (truncated) {
                std::fclose(f);
                if (err) *err = "'" + path + "' is truncated: not enough numbers for layer "
                                + (char) ('0' + L);
                return false;
            }
            if (k == "W") N.W[(size_t) L].swap(v); else N.B[(size_t) L].swap(v);
        }
    }
    std::fclose(f);
    if (!N.valid()) { if (err) *err = "'" + path + "' is not a usable surrogate file"; return false; }
    return true;
}

/* ================================================================================================
 *  Training
 * ================================================================================================ */
struct InputSpec {
    std::string name;
    double lo, hi;
    bool logScale;
    InputSpec() : lo(0), hi(1), logScale(false) {}
};

struct TrainOptions {
    std::vector<InputSpec> inputs;
    std::vector<int> hidden;
    int samples;
    int restarts;
    int epochs;
    int verifyPoints;
    bool logOutput;
    unsigned long long seed;

    TrainOptions() : samples(2000), restarts(5), epochs(4000), verifyPoints(512),
                     logOutput(false), seed(1) {
        hidden.push_back(10); hidden.push_back(10); hidden.push_back(10); hidden.push_back(10);
    }
};

struct TrainReport {
    int nSamples, nGrowing;
    double rmseTrain, rmseTest, r2Test, maxErrTest;
    double verifyRmse, verifyMaxRel;
    std::string note;
    TrainReport() : nSamples(0), nGrowing(0), rmseTrain(0), rmseTest(0), r2Test(0),
                    maxErrTest(0), verifyRmse(0), verifyMaxRel(0) {}
};

typedef double (*FbaFn)(const std::vector<double> &uptake, void *ctx);

/* ---- Latin hypercube, so every 1-D projection is covered even with few samples ---- */
inline void sampleInputs(const TrainOptions &opt, int n, Rng &rng,
                         std::vector<std::vector<double> > &X)
{
    const int d = (int) opt.inputs.size();
    X.assign((size_t) n, std::vector<double>((size_t) d, 0.0));
    std::vector<int> perm((size_t) n);
    for (int j = 0; j < d; ++j) {
        for (int i = 0; i < n; ++i) perm[(size_t) i] = i;
        for (int i = n - 1; i > 0; --i) {                    // Fisher-Yates
            const int k = (int) (rng.uniform() * (double) (i + 1));
            const int t = perm[(size_t) i]; perm[(size_t) i] = perm[(size_t) k]; perm[(size_t) k] = t;
        }
        const InputSpec &s = opt.inputs[(size_t) j];
        for (int i = 0; i < n; ++i) {
            const double u = ((double) perm[(size_t) i] + rng.uniform()) / (double) n;
            X[(size_t) i][(size_t) j] = s.logScale
                ? std::pow(10.0, std::log10(s.lo) + u * (std::log10(s.hi) - std::log10(s.lo)))
                : s.lo + u * (s.hi - s.lo);
        }
    }
}

namespace detail {

struct Fit {
    std::vector<int> sizes;
    std::vector<double> theta;                  // all weights then biases, layer by layer

    int nParam() const {
        int n = 0;
        for (size_t L = 0; L + 1 < sizes.size(); ++L) n += sizes[L + 1] * sizes[L] + sizes[L + 1];
        return n;
    }
    void unpack(std::vector<std::vector<double> > &W, std::vector<std::vector<double> > &B) const {
        W.clear(); B.clear();
        size_t k = 0;
        for (size_t L = 0; L + 1 < sizes.size(); ++L) {
            const int r = sizes[L + 1], c = sizes[L];
            W.push_back(std::vector<double>(theta.begin() + (long) k, theta.begin() + (long) k + r * c)); k += (size_t) (r * c);
            B.push_back(std::vector<double>(theta.begin() + (long) k, theta.begin() + (long) k + r));     k += (size_t) r;
        }
    }
};

/* Forward pass, keeping the activations backprop needs. */
inline double forward(const Fit &F, const std::vector<double> &x,
                      std::vector<std::vector<double> > &act)
{
    act.clear();
    act.push_back(x);
    size_t k = 0;
    for (size_t L = 0; L + 1 < F.sizes.size(); ++L) {
        const int r = F.sizes[L + 1], c = F.sizes[L];
        const double *w = &F.theta[k]; k += (size_t) (r * c);
        const double *b = &F.theta[k]; k += (size_t) r;
        std::vector<double> out((size_t) r);
        const std::vector<double> &in = act.back();
        for (int i = 0; i < r; ++i) {
            double s = b[i];
            for (int j = 0; j < c; ++j) s += w[i * c + j] * in[(size_t) j];
            out[(size_t) i] = (L + 2 == F.sizes.size()) ? s : Network::tansig(s);
        }
        act.push_back(out);
    }
    return act.back()[0];
}

/* Mean squared error and its gradient. */
inline double lossGrad(const Fit &F, const std::vector<std::vector<double> > &X,
                       const std::vector<double> &y, const std::vector<int> &idx,
                       std::vector<double> &grad)
{
    grad.assign(F.theta.size(), 0.0);
    double loss = 0.0;
    const double invN = 1.0 / (double) idx.size();
    std::vector<std::vector<double> > act;

    for (size_t t = 0; t < idx.size(); ++t) {
        const size_t s = (size_t) idx[t];
        const double out = forward(F, X[s], act);
        const double e = out - y[s];
        loss += e * e * invN;

        std::vector<double> delta(1, 2.0 * e * invN);       // d loss / d pre-activation, output
        long k = (long) F.theta.size();
        for (long L = (long) F.sizes.size() - 2; L >= 0; --L) {
            const int r = F.sizes[(size_t) L + 1], c = F.sizes[(size_t) L];
            k -= r;               const long bOff = k;
            k -= (long) r * c;    const long wOff = k;
            const std::vector<double> &in = act[(size_t) L];

            for (int i = 0; i < r; ++i) {
                grad[(size_t) (bOff + i)] += delta[(size_t) i];
                for (int j = 0; j < c; ++j)
                    grad[(size_t) (wOff + i * c + j)] += delta[(size_t) i] * in[(size_t) j];
            }
            if (L > 0) {
                std::vector<double> nd((size_t) c, 0.0);
                for (int j = 0; j < c; ++j) {
                    double s = 0.0;
                    for (int i = 0; i < r; ++i) s += F.theta[(size_t) (wOff + i * c + j)] * delta[(size_t) i];
                    const double a = in[(size_t) j];
                    nd[(size_t) j] = s * (1.0 - a * a);      // tansig'
                }
                delta.swap(nd);
            }
        }
    }
    return loss;
}

}  // namespace detail

/* ------------------------------------------------------------------------------------------------
 *  Fit a network to (X, y). Adam with restarts and early stopping on a held-out split.
 * ------------------------------------------------------------------------------------------------ */
inline void fitNetwork(Network &N, const std::vector<std::vector<double> > &Xraw,
                       const std::vector<double> &yraw, const TrainOptions &opt,
                       TrainReport &rep)
{
    const int n = (int) Xraw.size();
    const int d = (int) opt.inputs.size();

    /* ---- mapminmax, exactly as MATLAB defines it ---- */
    N.nIn = d;
    N.logOutput = opt.logOutput;
    N.xOffset.assign((size_t) d, 0.0); N.xGain.assign((size_t) d, 1.0); N.xYmin = -1.0;
    N.trainMin.assign((size_t) d, 0.0); N.trainMax.assign((size_t) d, 0.0);
    for (int j = 0; j < d; ++j) {
        double lo = Xraw[0][(size_t) j], hi = lo;
        for (int i = 1; i < n; ++i) {
            const double v = Xraw[(size_t) i][(size_t) j];
            if (v < lo) lo = v;
            if (v > hi) hi = v;
        }
        N.trainMin[(size_t) j] = lo; N.trainMax[(size_t) j] = hi;
        N.xOffset[(size_t) j] = lo;
        N.xGain[(size_t) j] = (hi > lo) ? 2.0 / (hi - lo) : 1.0;
    }

    std::vector<double> yfit(yraw);
    N.outMin = N.outMax = yraw.empty() ? 0.0 : yraw[0];
    for (int i = 0; i < n; ++i) { if (yraw[(size_t) i] < N.outMin) N.outMin = yraw[(size_t) i];
                                  if (yraw[(size_t) i] > N.outMax) N.outMax = yraw[(size_t) i]; }
    if (opt.logOutput) {
        const double floorv = (N.outMax > 0 ? N.outMax * 1e-9 : 1e-12);
        for (int i = 0; i < n; ++i) yfit[(size_t) i] = std::log10(yraw[(size_t) i] > floorv ? yraw[(size_t) i] : floorv);
    }
    double ylo = yfit[0], yhi = yfit[0];
    for (int i = 1; i < n; ++i) { if (yfit[(size_t) i] < ylo) ylo = yfit[(size_t) i];
                                  if (yfit[(size_t) i] > yhi) yhi = yfit[(size_t) i]; }
    N.yOffset = ylo; N.yGain = (yhi > ylo) ? 2.0 / (yhi - ylo) : 1.0; N.yYmin = -1.0;

    std::vector<std::vector<double> > Xs((size_t) n, std::vector<double>((size_t) d));
    std::vector<double> ys((size_t) n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < d; ++j)
            Xs[(size_t) i][(size_t) j] = (Xraw[(size_t) i][(size_t) j] - N.xOffset[(size_t) j]) * N.xGain[(size_t) j] + N.xYmin;
        ys[(size_t) i] = (yfit[(size_t) i] - N.yOffset) * N.yGain + N.yYmin;
    }

    /* ---- 70 / 15 / 15, the same split MATLAB's dividerand uses by default ---- */
    Rng rng(opt.seed);
    std::vector<int> perm((size_t) n);
    for (int i = 0; i < n; ++i) perm[(size_t) i] = i;
    for (int i = n - 1; i > 0; --i) {
        const int k = (int) (rng.uniform() * (double) (i + 1));
        const int t = perm[(size_t) i]; perm[(size_t) i] = perm[(size_t) k]; perm[(size_t) k] = t;
    }
    const int n1 = (int) (0.70 * n), n2 = (int) (0.85 * n);
    std::vector<int> itr(perm.begin(), perm.begin() + n1);
    std::vector<int> iva(perm.begin() + n1, perm.begin() + n2);
    std::vector<int> ite(perm.begin() + n2, perm.end());
    if (iva.empty()) iva = itr;
    if (ite.empty()) ite = itr;

    detail::Fit best;
    double bestVal = 1e300;

    for (int attempt = 0; attempt < opt.restarts; ++attempt) {
        detail::Fit F;
        F.sizes.push_back(d);
        for (size_t i = 0; i < opt.hidden.size(); ++i) F.sizes.push_back(opt.hidden[i]);
        F.sizes.push_back(1);

        Rng r2(opt.seed + 1000ULL * (unsigned long long) (attempt + 1));
        F.theta.assign((size_t) F.nParam(), 0.0);
        {   /* fan-in scaling, so tansig units start responsive rather than saturated */
            size_t k = 0;
            for (size_t L = 0; L + 1 < F.sizes.size(); ++L) {
                const int rr = F.sizes[L + 1], cc = F.sizes[L];
                const double sd = std::sqrt(2.0 / (double) (rr + cc));
                for (int i = 0; i < rr * cc; ++i) F.theta[k++] = r2.normal() * sd;
                for (int i = 0; i < rr; ++i) F.theta[k++] = 0.0;
            }
        }

        /* Adam. Chosen over L-BFGS because it needs no line search and no history, which keeps
         * this header short and its behaviour identical on every platform. */
        std::vector<double> m(F.theta.size(), 0.0), v(F.theta.size(), 0.0), g;
        const double lr = 0.01, b1 = 0.9, b2 = 0.999, eps = 1e-8;
        double localBestVal = 1e300;
        std::vector<double> localBestTheta = F.theta;
        int stale = 0;

        for (int ep = 1; ep <= opt.epochs; ++ep) {
            detail::lossGrad(F, Xs, ys, itr, g);
            const double c1 = 1.0 - std::pow(b1, (double) ep), c2 = 1.0 - std::pow(b2, (double) ep);
            for (size_t k = 0; k < F.theta.size(); ++k) {
                m[k] = b1 * m[k] + (1 - b1) * g[k];
                v[k] = b2 * v[k] + (1 - b2) * g[k] * g[k];
                F.theta[k] -= lr * (m[k] / c1) / (std::sqrt(v[k] / c2) + eps);
            }
            if (ep % 25 == 0 || ep == opt.epochs) {
                std::vector<double> dummy;
                const double vl = detail::lossGrad(F, Xs, ys, iva, dummy);
                if (vl < localBestVal - 1e-12) { localBestVal = vl; localBestTheta = F.theta; stale = 0; }
                else if (++stale > 20) break;        // early stopping
            }
        }
        F.theta = localBestTheta;
        if (localBestVal < bestVal) { bestVal = localBestVal; best = F; }
    }

    best.unpack(N.W, N.B);
    N.sizes = best.sizes;

    /* ---- report, on data the network never trained on ---- */
    double se = 0.0, mx = 0.0, mean = 0.0, ss = 0.0;
    for (size_t t = 0; t < ite.size(); ++t) mean += yraw[(size_t) ite[t]];
    mean /= (double) ite.size();
    for (size_t t = 0; t < ite.size(); ++t) {
        const size_t s = (size_t) ite[t];
        const double e = N.eval(Xraw[s]) - yraw[s];
        se += e * e;
        if (std::fabs(e) > mx) mx = std::fabs(e);
        ss += (yraw[s] - mean) * (yraw[s] - mean);
    }
    rep.nSamples = n;
    rep.rmseTest = std::sqrt(se / (double) ite.size());
    rep.maxErrTest = mx;
    rep.r2Test = (ss > 0) ? 1.0 - se / ss : 0.0;
    double se2 = 0.0;
    for (size_t t = 0; t < itr.size(); ++t) {
        const size_t s = (size_t) itr[t];
        const double e = N.eval(Xraw[s]) - yraw[s];
        se2 += e * e;
    }
    rep.rmseTrain = std::sqrt(se2 / (double) itr.size());
}

/* ------------------------------------------------------------------------------------------------
 *  The whole thing: sweep the LP, fit, then check the fit by re-solving the LP at fresh points.
 * ------------------------------------------------------------------------------------------------ */
inline bool train(Network &N, const TrainOptions &opt, FbaFn fba, void *ctx, TrainReport &rep,
                  std::string *err = 0)
{
    if (opt.inputs.empty()) { if (err) *err = "no <inputN> lines: nothing to train on"; return false; }
    for (size_t j = 0; j < opt.inputs.size(); ++j) {
        const InputSpec &s = opt.inputs[j];
        if (s.hi <= s.lo) { if (err) *err = "input '" + s.name + "' has hi <= lo"; return false; }
        if (s.logScale && s.lo <= 0) { if (err) *err = "input '" + s.name + "' is logarithmic but starts at or below zero"; return false; }
    }

    Rng rng(opt.seed);
    std::vector<std::vector<double> > X;
    sampleInputs(opt, opt.samples, rng, X);

    std::vector<double> y((size_t) opt.samples, 0.0);
    int growing = 0;
    for (int i = 0; i < opt.samples; ++i) {
        double g = fba(X[(size_t) i], ctx);
        if (!(g > 0.0) || g != g) g = 0.0;       // infeasible or NaN is zero growth, not an error
        y[(size_t) i] = g;
        if (g > 1e-12) ++growing;
    }
    rep.nGrowing = growing;
    if (growing == 0) {
        if (err) *err = "nothing grew anywhere in the sweep. Widen the ranges, check the exchange "
                        "indices, or check that the model's objective is the biomass reaction.";
        return false;
    }

    fitNetwork(N, X, y, opt, rep);

    /* ---- verification: fresh points, solved by the LP, compared against the network ----
     * This is the step that makes the result trustworthy. The fit statistics above are computed on
     * a held-out split of the SAME sweep; this re-solves the linear program at points drawn
     * independently, which is the only way to catch a network that has learned the sampling rather
     * than the function. */
    if (opt.verifyPoints > 0) {
        Rng vr(opt.seed + 99991ULL);
        std::vector<std::vector<double> > Xv;
        sampleInputs(opt, opt.verifyPoints, vr, Xv);
        double se = 0.0, worstRel = 0.0;
        for (int i = 0; i < opt.verifyPoints; ++i) {
            const double truth = fba(Xv[(size_t) i], ctx);
            const double pred = N.eval(Xv[(size_t) i]);
            const double e = pred - truth;
            se += e * e;
            const double denom = (std::fabs(truth) > 1e-12) ? std::fabs(truth) : 1e-12;
            const double rel = std::fabs(e) / denom;
            if (truth > 1e-9 && rel > worstRel) worstRel = rel;
        }
        rep.verifyRmse = std::sqrt(se / (double) opt.verifyPoints);
        rep.verifyMaxRel = worstRel;
    }
    return true;
}

inline std::string reportText(const Network &N, const TrainReport &r)
{
    char b[1024];
    std::string s;
    std::snprintf(b, sizeof(b), "    samples %d (%d with growth), architecture", r.nSamples, r.nGrowing);
    s += b;
    for (size_t i = 0; i < N.sizes.size(); ++i) { std::snprintf(b, sizeof(b), " %d", N.sizes[i]); s += b; }
    s += "\n";
    std::snprintf(b, sizeof(b), "    held-out RMSE %.6g, R^2 %.6f, worst error %.6g\n",
                  r.rmseTest, r.r2Test, r.maxErrTest);
    s += b;
    std::snprintf(b, sizeof(b), "    verification against fresh LP solves: RMSE %.6g, worst relative %.3g\n",
                  r.verifyRmse, r.verifyMaxRel);
    s += b;
    s += "    valid input range (outside it the network extrapolates and means nothing):\n";
    for (int i = 0; i < N.nIn; ++i) {
        std::snprintf(b, sizeof(b), "      input %d : %.10g .. %.10g\n", i, N.trainMin[(size_t) i], N.trainMax[(size_t) i]);
        s += b;
    }
    return s;
}

}  // namespace complab_srg

#endif  // COMPLAB3D_SURROGATE_HH
