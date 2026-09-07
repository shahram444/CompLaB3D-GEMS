/* ================================================================================================
 *  test_surrogate_parity.cpp -- the two surrogate evaluators must agree
 *
 *  THE SAME ARITHMETIC IS WRITTEN TWICE IN THIS REPOSITORY, and that is worth being honest about:
 *
 *    1. src/complab3d_surrogate.hh  Network::evalAll() -- generic, loops over the layer sizes it
 *       read from the file, so one function evaluates any architecture.  This is the run-time path:
 *       it is what <weights_file> loads.
 *
 *    2. the generated header from trainSurrogate.py -- the same forward pass with the loop bounds
 *       written out as literals and the weights as fixed-size arrays.  This is the compiled path,
 *       and it exists because it is the format MATLAB's genFunction emits, which is what the 2D
 *       CompLaB used and what every network exported before this repository looks like.
 *
 *  Two implementations of one calculation is a liability unless something checks they agree.  This
 *  is that check: the SAME fit is written to both formats in one pass by trainSurrogate.py, then
 *  loaded here through the run-time reader and compared against the reference table the trainer
 *  produced from its own arithmetic.  A divergence means one of the two has drifted, and it would
 *  otherwise show up as a run whose results depend on whether the network was compiled in or read
 *  from a file -- which nobody would think to suspect.
 *
 *  The test needs a .srg and its matching reference table.  Both are produced by
 *
 *      python3 tools/surrogate/trainSurrogate.py sweep.csv -o net.hh --srg net.srg
 *
 *  and the test SKIPS, loudly, when they are absent, because the sweep needs cobra and a metabolic
 *  model and neither belongs in a unit test.
 * ================================================================================================ */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include "complab3d_surrogate.hh"

int main(int argc, char **argv)
{
    const std::string net = (argc > 1) ? argv[1] : "parity_net.srg";
    const std::string ref = (argc > 2) ? argv[2] : "parity_net_reference.csv";

    std::FILE *probe = std::fopen(net.c_str(), "r");
    if (!probe) {
        std::printf("  SKIPPED: %s not present.\n", net.c_str());
        std::printf("  Produce it with:  python3 tools/surrogate/trainSurrogate.py \\\n");
        std::printf("                        sweep.csv -o net.hh --srg %s\n", net.c_str());
        std::printf("  A skipped check is not a pass, and this one is skipped because the\n");
        std::printf("  sweep needs cobra and a metabolic model, not because it is optional.\n");
        return 0;
    }
    std::fclose(probe);

    complab_srg::Network N;
    std::string err;
    if (!complab_srg::load(N, net, &err)) {
        std::printf("  FAIL: %s\n", err.c_str());
        return 1;
    }
    const int nIn = N.nIn, nOut = N.nOut();
    std::printf("  loaded %s: %d input(s), %d output(s)\n", net.c_str(), nIn, nOut);

    std::FILE *f = std::fopen(ref.c_str(), "r");
    if (!f) { std::printf("  FAIL: %s is missing\n", ref.c_str()); return 1; }

    char line[65536];
    int rows = 0;
    std::vector<double> worst((size_t) nOut, 0.0);
    while (std::fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        std::vector<double> v;
        char *p = line;
        while (*p) {
            char *end = 0;
            const double x = std::strtod(p, &end);
            if (end == p) break;
            v.push_back(x);
            p = end;
            while (*p == ',' || *p == ' ') ++p;
        }
        if ((int) v.size() < nIn + nOut) continue;
        std::vector<double> x(v.begin(), v.begin() + nIn);
        std::vector<double> got;
        N.evalAll(x, got);
        for (int k = 0; k < nOut; ++k) {
            const double want = v[(size_t) (nIn + k)];
            const double d = std::fabs(got[(size_t) k] - want) / (std::fabs(want) > 1e-30 ? std::fabs(want) : 1e-30);
            if (d > worst[(size_t) k]) worst[(size_t) k] = d;
        }
        ++rows;
    }
    std::fclose(f);

    if (rows == 0) { std::printf("  FAIL: no usable rows in %s\n", ref.c_str()); return 1; }

    std::printf("  compared %d point(s) against the trainer's own arithmetic\n", rows);
    double overall = 0.0;
    for (int k = 0; k < nOut; ++k) {
        const char *nm = (k < (int) N.outputNames.size()) ? N.outputNames[(size_t) k].c_str() : "output";
        std::printf("    %-20s worst relative difference %.3e\n", nm, worst[(size_t) k]);
        if (worst[(size_t) k] > overall) overall = worst[(size_t) k];
    }

    if (overall > 1e-10) {
        std::printf("\n  FAIL: the run-time evaluator does not reproduce the trainer.\n");
        std::printf("  One of the two forward passes has drifted from the other.\n");
        return 1;
    }
    std::printf("\n  the run-time and compiled paths agree to %.1e relative\n", overall);
    return 0;
}
