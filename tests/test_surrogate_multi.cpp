/* ================================================================================================
 *  test_surrogate_multi.cpp -- the surrogate network returns more than growth
 *
 *  A single-output surrogate predicts the growth rate and nothing else, so the solver falls back
 *  to a Monod term for how much of each substrate was consumed.  That fallback is exact only where
 *  the swept uptake bound was the binding constraint.  Everywhere else it is wrong, and a
 *  growth-only network cannot release a product at all.
 *
 *  These checks cover the three ways that can go wrong silently:
 *    1. a multi-output file that does not survive a save/load round trip;
 *    2. a single-output file, written before any of this, that stops loading;
 *    3. an inconsistent file being accepted rather than refused.
 * ================================================================================================ */
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "complab3d_surrogate.hh"

static int failures = 0;

static void check(const char *what, bool ok)
{
    std::printf("  %-58s %s\n", what, ok ? "ok" : "FAIL");
    if (!ok) ++failures;
}

static void checkNear(const char *what, double got, double want, double tol)
{
    const bool ok = std::fabs(got - want) <= tol;
    std::printf("  %-58s got %-16.10g want %-16.10g %s\n", what, got, want, ok ? "ok" : "FAIL");
    if (!ok) ++failures;
}

/* A small network with a known answer: identity input scaling, one hidden layer, three outputs. */
static complab_srg::Network makeNet(int nOut)
{
    complab_srg::Network N;
    N.nIn = 2;
    N.sizes.push_back(2); N.sizes.push_back(4); N.sizes.push_back(nOut);
    N.W.resize(2); N.B.resize(2);
    N.W[0].assign(8, 0.3); N.B[0].assign(4, 0.05);
    N.W[1].assign((size_t)(4 * nOut), 0.2); N.B[1].assign((size_t) nOut, -0.1);
    N.xOffset.assign(2, 0.0); N.xGain.assign(2, 1.0); N.xYmin = -1.0;
    for (int k = 0; k < nOut; ++k) { N.yOffset.push_back(0.1 * k); N.yGain.push_back(1.0 + k); }
    N.yYmin = -1.0;
    N.trainMin.assign(2, 0.0); N.trainMax.assign(2, 10.0);
    N.outMin.assign((size_t) nOut, 0.0); N.outMax.assign((size_t) nOut, 1.0);
    N.outputNames.push_back("growth");
    if (nOut > 1) N.outputNames.push_back("r_acetate");
    if (nOut > 2) N.outputNames.push_back("r_Fe3");
    N.provenance = "test_surrogate_multi";
    return N;
}

int main()
{
    std::string err;
    const char *path = "multi_test.srg";
    std::vector<double> x(2); x[0] = 3.0; x[1] = 4.0;

    std::printf("--- a three-output network survives the file unchanged\n");
    complab_srg::Network N = makeNet(3);
    check("it reports three outputs", N.nOut() == 3);
    check("save", complab_srg::save(N, path, &err));

    complab_srg::Network M;
    check("load", complab_srg::load(M, path, &err));
    check("the loaded net still reports three outputs", M.nOut() == 3);
    check("the output names came back", M.outputNames.size() == 3
                                        && M.outputNames[1] == "r_acetate");

    std::vector<double> a, b;
    N.evalAll(x, a);
    M.evalAll(x, b);
    check("three values came back", a.size() == 3 && b.size() == 3);
    double worst = 0.0;
    for (size_t k = 0; k < a.size() && k < b.size(); ++k)
        worst = std::max(worst, std::fabs(a[k] - b[k]));
    checkNear("round trip is exact, not merely close", worst, 0.0, 0.0);

    /* The whole point of the per-output scaling: two outputs with different
     * gains must NOT come back equal.  If they do, the scaling is being applied
     * from a single scalar and every flux is wrong by a constant factor. */
    check("outputs are scaled independently", std::fabs(b[0] - b[1]) > 1e-12);

    std::printf("\n--- growth still means growth\n");
    checkNear("eval() returns output 0 and nothing else", M.eval(x), b[0], 0.0);

    std::printf("\n--- a single-output network is unaffected\n");
    complab_srg::Network S = makeNet(1);
    check("save", complab_srg::save(S, "single_test.srg", &err));
    complab_srg::Network T;
    check("load", complab_srg::load(T, "single_test.srg", &err));
    check("it reports one output", T.nOut() == 1);
    std::vector<double> t;
    T.evalAll(x, t);
    check("evalAll returns exactly one value", t.size() == 1);
    checkNear("and it agrees with eval()", T.eval(x), t[0], 0.0);

    std::printf("\n--- an inconsistent file is refused, not guessed at\n");
    /* Claim four outputs in a file whose last layer has three.  Accepting this
     * would write fluxes into the wrong substrates, silently. */
    {
        std::FILE *f = std::fopen(path, "r");
        std::string body;
        int c; while ((c = std::fgetc(f)) != EOF) body.push_back((char) c);
        std::fclose(f);
        const std::string from = "outputs 3", to = "outputs 4";
        const size_t at = body.find(from);
        if (at != std::string::npos) body.replace(at, from.size(), to);
        f = std::fopen("bad_test.srg", "w");
        std::fwrite(body.data(), 1, body.size(), f);
        std::fclose(f);
    }
    complab_srg::Network Bad;
    const bool loaded = complab_srg::load(Bad, "bad_test.srg", &err);
    check("a file whose declared output count is a lie is rejected", !loaded);
    if (!loaded) std::printf("     reported: %s\n", err.c_str());

    std::remove(path); std::remove("single_test.srg"); std::remove("bad_test.srg");

    std::printf("\n%s\n", failures ? "SOME CHECKS FAILED" : "all multi-output checks passed");
    return failures ? 1 : 0;
}
