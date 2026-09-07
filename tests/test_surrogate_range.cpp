/* ================================================================================================
 *  test_surrogate_range.cpp  --  the COMPILED surrogate must not extrapolate silently
 *
 *  Every other fitted rate path in this repository holds its inputs to the box it was fitted over
 *  and counts how often it had to: a .sym law clamps against its `range` lines, a .gnn network
 *  against `trainmin`/`trainmax`, and a .srg network against the same, through evalBound().
 *
 *  The COMPILED surrogate -- weights pasted into surrogateModel.hh and built into the executable --
 *  did none of that. Outside its training box it returned a confident number, nothing counted it,
 *  and nothing in the output distinguished it from a real one. That is the single failure a fitted
 *  model has that a hand-written rate law does not, and it was the only path where it was invisible.
 *
 *  This file compiles the shipped surrogateModel.hh against the Palabos stub and checks the three
 *  things that make the fix real rather than nominal:
 *
 *      the answer outside the box equals the answer AT the box, so the clamp is the one that ran;
 *      an evaluation that clamped is counted once, not once per input;
 *      the end-of-run report names the compiled path and warns when the fraction is large.
 *
 *      g++ -O2 -Wall -Wextra -std=c++11 -I. -I.. -I../src -o t test_surrogate_range.cpp && ./t
 * ================================================================================================ */

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "palabos_stub.hh"
#include "../src/surrogateModel.hh"

static int fails = 0;

static void ck(const char *what, double got, double want, double tol = 1e-12) {
    const bool ok = std::fabs(got - want) <= tol * (1.0 + std::fabs(want));
    if (!ok) ++fails;
    std::printf("  %-56s got %-15.9g want %-15.9g %s\n",
                what, got, want, ok ? "ok" : "** FAIL **");
}

static void ckTrue(const char *what, bool ok) {
    if (!ok) ++fails;
    std::printf("  %-56s %s\n", what, ok ? "ok" : "** FAIL **");
}

/* The training box of the shipped network, as the comment in surrogateModel.hh states it and as
 * the clamp is now built from: lo = offset, hi = offset + (1 - ymin)/gain. Written out here a
 * second time on purpose -- if the header's numbers drift from these, that is a real change to
 * what the network claims to be valid for, and it should fail here rather than pass quietly. */
static const double LO0 = 0.000890159834356918;
static const double LO1 = 3.51303254819135e-05;
static const double HI0 = 0.000890159834356918 + 2.0 / 0.200059307889652;   /* ~10.0 */
static const double HI1 = 3.51303254819135e-05 + 2.0 / 4.00231590228802;    /* ~0.4998 */

/* One call to the shipped model, for microbe 0, with two uptake bounds. */
static double growthAt(double donor, double acceptor)
{
    std::vector< std::vector<double> > Fin(1, std::vector<double>(2, 0.0));
    std::vector< std::vector<double> > Fout(1, std::vector<double>(2, 0.0));
    std::vector<double> bioR(1, 0.0);
    Fin[0][0] = donor;
    Fin[0][1] = acceptor;
    Fout = Fin;
    defineSurrogateModel(0, Fin, Fout, bioR, 2);
    return bioR[0];
}

static long evals()   { return complab_srg::runtime().compiledEvaluations; }
static long clamped() { return complab_srg::runtime().compiledClamped; }

int main()
{
    complab_srg::runtime() = complab_srg::Runtime();

    std::printf("--- inside the box, the network answers and nothing is clamped\n");
    {
        const double mid0 = 0.5 * (LO0 + HI0), mid1 = 0.5 * (LO1 + HI1);
        const double g = growthAt(mid0, mid1);
        ckTrue("a mid-box point gives a finite, non-negative growth rate",
               g == g && g >= 0.0 && g < 1e3);
        ckTrue("the evaluation was counted", evals() == 1);
        ckTrue("and it was not counted as clamped", clamped() == 0);
        std::printf("       growth at the centre of the box: %.6g 1/h\n", g);
    }

    std::printf("--- outside the box, the answer is the answer AT the box\n");
    {
        complab_srg::runtime() = complab_srg::Runtime();

        /* Far above the donor ceiling. If the clamp is working, this is the same number as the
         * ceiling itself -- and if it is NOT working, it is an extrapolation, which for a tanh
         * network several boundary-widths out is a different number entirely. */
        const double atEdge  = growthAt(HI0, HI1);
        const double wayOut  = growthAt(HI0 * 50.0, HI1 * 50.0);
        ck("50x past both ceilings equals the value at the ceiling", wayOut, atEdge);

        const double atFloor = growthAt(LO0, LO1);
        const double below   = growthAt(-1.0, -1.0);
        ck("below both floors equals the value at the floor", below, atFloor);

        ckTrue("all four evaluations were counted", evals() == 4);
        ckTrue("exactly the two out-of-box ones were clamped", clamped() == 2);
    }

    std::printf("--- an evaluation that clamps two inputs still counts once\n");
    {
        complab_srg::runtime() = complab_srg::Runtime();
        growthAt(HI0 * 10.0, -5.0);              /* one above, one below */
        ckTrue("one evaluation", evals() == 1);
        /* Counting values rather than evaluations is how a two-input network comes to report
         * 200% of its evaluations as out of range, which makes the one number a reader uses to
         * judge the result unreadable. The .sym and .gnn loaders had exactly this bug. */
        ckTrue("one clamp, not two", clamped() == 1);
    }

    std::printf("--- only one input outside is still one clamp\n");
    {
        complab_srg::runtime() = complab_srg::Runtime();
        const double mid1 = 0.5 * (LO1 + HI1);
        growthAt(HI0 * 3.0, mid1);
        ckTrue("counted as clamped", clamped() == 1);
        const double g1 = growthAt(HI0, mid1);
        const double g2 = growthAt(HI0 * 3.0, mid1);
        ck("and the clamped input is the only one that moved", g2, g1);
    }

    std::printf("--- the closing report names the compiled path\n");
    {
        complab_srg::runtime() = complab_srg::Runtime();
        ckTrue("a run that evaluated nothing says nothing",
               complab_srg::runtimeReport().empty());

        const double mid0 = 0.5 * (LO0 + HI0), mid1 = 0.5 * (LO1 + HI1);
        for (int i = 0; i < 100; ++i) growthAt(mid0, mid1);
        std::string r = complab_srg::runtimeReport();
        ckTrue("100 clean evaluations are reported",
               r.find("compiled network") != std::string::npos);
        ckTrue("and no warning is raised", r.find("OUTSIDE THE TRAINING RANGE") == std::string::npos);

        for (int i = 0; i < 20; ++i) growthAt(HI0 * 5.0, mid1);      /* 20/120 = 16.7% */
        r = complab_srg::runtimeReport();
        ckTrue("a large out-of-box fraction is warned about",
               r.find("OUTSIDE THE TRAINING RANGE") != std::string::npos);
        std::printf("%s", r.c_str());
    }

    std::printf("--- the runtime and compiled counters do not share a total\n");
    {
        complab_srg::runtime() = complab_srg::Runtime();
        complab_srg::runtime().evaluations = 7;      /* as if a .srg network had run */
        complab_srg::runtime().clamped = 1;
        const double mid0 = 0.5 * (LO0 + HI0), mid1 = 0.5 * (LO1 + HI1);
        growthAt(mid0, mid1);
        ckTrue("the compiled call did not touch the runtime counters",
               complab_srg::runtime().evaluations == 7 && complab_srg::runtime().clamped == 1);
        const std::string r = complab_srg::runtimeReport();
        ckTrue("and the report carries both lines",
               r.find("runtime network") != std::string::npos
               && r.find("compiled network") != std::string::npos);
    }

    std::printf("\n%s\n", fails ? "SOME CHECKS FAILED" : "all compiled-surrogate range checks passed");
    return fails;
}
