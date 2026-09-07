/* ================================================================================================
 *  test_gnn.cpp  --  complab3d_graphnet.hh checked without Palabos and without a trainer
 *
 *  tests/xval_gnn.py checks the C++ forward pass against the Python one.  That catches a
 *  disagreement between the two, but if both were wrong in the same way it would say nothing.  So
 *  this file does the other half: it builds a network small enough to work out on paper, writes it
 *  as a .gnn file, loads it, and compares eval() against the message-passing rule written out by
 *  hand as a single expression.  Different code, same arithmetic.
 *
 *      g++ -O2 -Wall -Wextra -std=c++11 -I.. -o t test_gnn.cpp && ./t
 * ================================================================================================ */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "complab3d_graphnet.hh"

static int fails = 0;

static void ck(const char *what, double got, double want, double tol = 1e-12) {
    const bool ok = std::fabs(got - want) <= tol * (1.0 + std::fabs(want));
    if (!ok) ++fails;
    std::printf("  %-52s got %-16.10g want %-16.10g %s\n",
                what, got, want, ok ? "ok" : "** FAIL **");
}

static void ckTrue(const char *what, bool ok) {
    if (!ok) ++fails;
    std::printf("  %-52s %s\n", what, ok ? "ok" : "** FAIL **");
}

static void write(const char *path, const char *body) {
    std::FILE *f = std::fopen(path, "w");
    std::fputs(body, f);
    std::fclose(f);
}

/* ------------------------------------------------------------------------------------------------
 *  Two species, one reaction, width one, one round.  Everything is identity-scaled so the numbers
 *  that come out of eval() are the raw network output and nothing is hidden in the scaling.
 *
 *      S = [ -1 ]      da = [0.5]      Wenc = [1]   benc = [0]
 *          [  2 ]
 *      Wsr = [1]  Wda = [2]  bR = [0]  Wss = [0.5]  Wrs = [1]  bS = [0]
 *      Wout = [1]  bout = [0]
 * ---------------------------------------------------------------------------------------------- */
static const char *TINY =
    "# a network small enough to check by hand\n"
    "provenance written by test_gnn.cpp\n"
    "version 1\n"
    "units per_second\n"
    "species A B\n"
    "reactions R1\n"
    "rounds 1\n"
    "width 1\n"
    "growth 0\n"
    "stoich 2 1\n"
    "-1\n"
    "2\n"
    "da 0.5\n"
    "xoffset 0 0\n"
    "xgain 1 1\n"
    "xymin 0\n"
    "yoffset 0 0\n"
    "ygain 1 1\n"
    "yymin 0\n"
    "trainmin 0 0\n"
    "trainmax 10 10\n"
    "Wenc 1\n"
    "benc 0\n"
    "layer 0\n"
    "Wsr 1\n"
    "Wda 2\n"
    "bR 0\n"
    "Wss 0.5\n"
    "Wrs 1\n"
    "bS 0\n"
    "Wout 1\n"
    "bout 0\n";

/* the same rule written out directly, with no loops and no matrices */
static void byHand(double cA, double cB, double &oA, double &oB) {
    const double hA = std::tanh(cA);
    const double hB = std::tanh(cB);
    const double mR = (-1.0) * hA + (2.0) * hB;            /* stoichiometry-weighted */
    const double hR = std::tanh(1.0 * mR + 2.0 * 0.5);     /* Wsr*mR + Wda*da */
    oA = std::tanh(0.5 * hA + 1.0 * ((-1.0) * hR));        /* Wss*hA + Wrs*(s*hR) */
    oB = std::tanh(0.5 * hB + 1.0 * ((2.0) * hR));
}

int main() {
    std::printf("--- a network small enough to check by hand\n");
    write("tiny.gnn", TINY);

    complab_gnn::Network N;
    std::string err;
    ckTrue("loads", complab_gnn::load(N, "tiny.gnn", &err));
    if (!err.empty()) std::printf("     %s\n", err.c_str());
    ckTrue("2 species, 1 reaction, width 1, 1 round",
           N.nS == 2 && N.nR == 1 && N.width == 1 && N.rounds == 1);
    ckTrue("no growth output", !N.hasGrowth);
    ckTrue("readout has one row, not one per species", N.Wout.rows == 1);
    ckTrue("provenance kept, '#' comment dropped",
           N.provenance == "written by test_gnn.cpp");

    std::vector<double> c(2), out;
    double wA, wB;

    const double probe[4][2] = { {1.0, 2.0}, {0.0, 0.0}, {3.5, 0.25}, {10.0, 10.0} };
    for (int t = 0; t < 4; ++t) {
        c[0] = probe[t][0]; c[1] = probe[t][1];
        N.eval(c, out, true, 0);
        byHand(c[0], c[1], wA, wB);
        char lbl[96];
        std::sprintf(lbl, "rate A at (%.4g, %.4g)", c[0], c[1]); ck(lbl, out[0], wA);
        std::sprintf(lbl, "rate B at (%.4g, %.4g)", c[0], c[1]); ck(lbl, out[1], wB);
    }

    std::printf("--- the training box is enforced, not advised\n");
    long clamped = 0;
    c[0] = 50.0; c[1] = 2.0;                 /* A is way outside trainmax = 10 */
    N.eval(c, out, true, &clamped);
    byHand(10.0, 2.0, wA, wB);
    ck("out-of-box input evaluates at the boundary", out[0], wA);
    ck("and the other species is unaffected",        out[1], wB);
    ckTrue("the clamp was counted", clamped == 1);

    clamped = 0;
    c[0] = -1.0; c[1] = 99.0;                /* both outside, one low one high */
    N.eval(c, out, true, &clamped);
    /* ONE count, not two.  The counter is compared against the number of evaluations in
     * runtimeReport(), so it counts evaluations that clamped anything rather than clamped
     * values -- otherwise a four-species network can report 400% of its evaluations as
     * out of range, and the percentage a reader uses to judge a result becomes unreadable. */
    ckTrue("an evaluation that clamps two species still counts once", clamped == 1);

    clamped = 0;
    c[0] = 1.0; c[1] = 2.0;
    N.eval(c, out, true, &clamped);
    ckTrue("an in-box input is not counted", clamped == 0);
    ckTrue("inTrainingRange agrees",  N.inTrainingRange(c));
    c[0] = 50.0;
    ckTrue("inTrainingRange rejects", !N.inTrainingRange(c));

    std::printf("--- units are a factor of 3600, so they are checked\n");
    {
        std::string body(TINY);
        const size_t p = body.find("per_second");
        body.replace(p, 10, "per_hour");
        write("tiny_hr.gnn", body.c_str());

        complab_gnn::Network H;
        ckTrue("per_hour file loads", complab_gnn::load(H, "tiny_hr.gnn", &err));
        std::vector<double> o2;
        c[0] = 1.0; c[1] = 2.0;
        H.eval(c, o2, true, 0);
        N.eval(c, out, true, 0);
        ck("per_hour rate is the per_second one over 3600", o2[0] * 3600.0, out[0]);
        ckTrue("and the log says which was loaded",
               complab_gnn::describe(H, "tiny_hr.gnn").find("per_hour") != std::string::npos);
    }

    std::printf("--- a file that is wrong is refused, not half-loaded\n");
    {
        struct Case { const char *what; std::string body; };
        std::vector<Case> bad;

        std::string b;

        b = TINY;                                        /* bS never appears */
        b.erase(b.find("bS 0\n"), 5);
        { Case cse = { "a missing bS line", b }; bad.push_back(cse); }

        b = TINY;                                        /* layer weights with no layer block */
        b.erase(b.find("layer 0\n"), 8);
        { Case cse = { "layer weights outside any layer", b }; bad.push_back(cse); }

        b = TINY;                                        /* Wenc before the shape is known */
        b.insert(b.find("stoich"), "Wenc 1\n");
        { Case cse = { "a weight before the stoich block", b }; bad.push_back(cse); }

        b = TINY;
        b.replace(b.find("per_second"), 10, "per_fortnight");
        { Case cse = { "an unrecognised units word", b }; bad.push_back(cse); }

        b = TINY;
        b += "wibble 3\n";
        { Case cse = { "an unknown keyword", b }; bad.push_back(cse); }

        b = TINY;                                        /* no valid range declared */
        b.erase(b.find("trainmin 0 0\n"), 13);
        { Case cse = { "no trainmin, so no known valid range", b }; bad.push_back(cse); }

        b = TINY;                                        /* truncated stoichiometry */
        b.replace(b.find("stoich 2 1"), 10, "stoich 2 3");
        { Case cse = { "a stoich block that does not match its shape", b }; bad.push_back(cse); }

        for (size_t i = 0; i < bad.size(); ++i) {
            write("bad.gnn", bad[i].body.c_str());
            complab_gnn::Network B;
            std::string e;
            const bool loaded = complab_gnn::load(B, "bad.gnn", &e);
            if (loaded) ++fails;
            std::printf("  %-52s %s\n", bad[i].what,
                        loaded ? "** FAIL: accepted **" : "refused ok");
        }

        complab_gnn::Network M;
        ckTrue("a file that is not there is refused",
               !complab_gnn::load(M, "no_such_file.gnn", &err));
    }

    std::printf("--- the runtime report\n");
    {
        complab_gnn::Runtime &R = complab_gnn::runtime();
        R.evaluations = 1000; R.clamped = 0;
        ckTrue("a clean run reports no warning",
               complab_gnn::runtimeReport().find("WARNING") == std::string::npos);
        R.clamped = 200;                                  /* 20% */
        ckTrue("20% out of box warns",
               complab_gnn::runtimeReport().find("WARNING") != std::string::npos);
        R.evaluations = 0; R.clamped = 0;
        ckTrue("a run that never evaluated says nothing",
               complab_gnn::runtimeReport().empty());
    }

    std::printf("--- the EXTENT readout, where the stoichiometry is exact by construction\n");
    {
        /* The same tiny shape, read out the other way.  With S = [-1; 2] and one reaction, every
         * output pair must satisfy rate_B / rate_A = 2 / -1 exactly, at every input, whatever the
         * weights happen to be.  That is the whole claim of this mode, so it is checked at inputs
         * chosen to be awkward rather than at one convenient point.
         *
         *   xi   = tanh-free linear readout on the REACTION node:  bout + Wout . hR
         *   r_A  = -1 * xi * xiscale ,   r_B = +2 * xi * xiscale
         */
        std::string body(TINY);
        body.insert(body.find("stoich"), "readout extent\n");
        /* EXTENT mode carries no per-species scaling: xiscale replaces it. */
        body.replace(body.find("yoffset 0 0\n"), 12, "xiscale 3\n");
        body.replace(body.find("ygain 1 1\n"), 10, "");
        write("tiny_ext.gnn", body.c_str());

        complab_gnn::Network E;
        ckTrue("the extent file loads", complab_gnn::load(E, "tiny_ext.gnn", &err));
        if (!err.empty()) std::printf("     %s\n", err.c_str());
        ckTrue("the mode was read", E.readout == complab_gnn::Network::EXTENT);
        ckTrue("no per-species scaling is required", E.yOffset.empty() && E.yGain.empty());
        ckTrue("one extent scale, for the one reaction", E.xiScale.size() == 1);
        ckTrue("and the log says which readout it is",
               complab_gnn::describe(E, "tiny_ext.gnn").find("EXTENT") != std::string::npos);

        const double probe2[5][2] = { {1.0, 2.0}, {0.0, 0.0}, {3.5, 0.25},
                                      {10.0, 10.0}, {0.125, 9.75} };
        std::vector<double> c2(2), o2;
        for (int t = 0; t < 5; ++t) {
            c2[0] = probe2[t][0]; c2[1] = probe2[t][1];
            E.eval(c2, o2, true, 0);

            /* the same message passing written out by hand, then the extent readout */
            const double hA = std::tanh(c2[0]);
            const double hB = std::tanh(c2[1]);
            const double mR = (-1.0) * hA + (2.0) * hB;
            const double hR = std::tanh(1.0 * mR + 2.0 * 0.5);
            const double xi = (0.0 + 1.0 * hR) * 3.0;          /* bout + Wout.hR, times xiscale */

            char lbl[96];
            std::sprintf(lbl, "rate A at (%.4g, %.4g)", c2[0], c2[1]);
            ck(lbl, o2[0], -1.0 * xi);
            std::sprintf(lbl, "rate B at (%.4g, %.4g)", c2[0], c2[1]);
            ck(lbl, o2[1], 2.0 * xi);

            /* the claim, stated as the test: the ratio is the stoichiometry, exactly */
            if (std::fabs(o2[0]) > 1e-14) {
                std::sprintf(lbl, "B over A is exactly -2 at (%.4g, %.4g)", c2[0], c2[1]);
                ck(lbl, o2[1] / o2[0], -2.0, 1e-14);
            }
        }

        std::printf("--- a file labelled one mode and written for the other is refused\n");
        {
            /* Extent mode without xiscale: the species rates would have no scale at all. */
            std::string b2(TINY);
            b2.insert(b2.find("stoich"), "readout extent\n");
            b2.replace(b2.find("yoffset 0 0\n"), 12, "");
            b2.replace(b2.find("ygain 1 1\n"), 10, "");
            write("bad.gnn", b2.c_str());
            complab_gnn::Network B;
            std::string e;
            ckTrue("extent mode with no xiscale is refused",
                   !complab_gnn::load(B, "bad.gnn", &e));

            /* Extent mode still carrying the species-mode scaling: two entries where the mode
             * wants none, which would silently consume the next key's numbers. */
            std::string b3(TINY);
            b3.insert(b3.find("stoich"), "readout extent\n");
            b3 += "xiscale 3\n";
            write("bad.gnn", b3.c_str());
            complab_gnn::Network B2;
            ckTrue("extent mode with species-mode yoffset/ygain is refused",
                   !complab_gnn::load(B2, "bad.gnn", &e));
        }

        std::remove("tiny_ext.gnn");
    }

    std::remove("tiny.gnn"); std::remove("tiny_hr.gnn"); std::remove("bad.gnn");
    std::printf("\n%s\n", fails ? "SOME CHECKS FAILED" : "all graph network checks passed");
    return fails;
}
