/* ================================================================================================
 *  test_upscale.cpp  --  complab3d_upscale.hh checked without Palabos
 *
 *  This file exists because of one bug.  classicalEta() was written for the GENERALIZED Thiele
 *  modulus, (R/3) sqrt(k/D), while complab.cpp computed and every README printed the modulus built
 *  on the sphere RADIUS, R sqrt(k/D).  The two differ by a factor of three in the argument, so the
 *  classical effectiveness factor came out at eta(3 phi): about 0.15 too low across the range the
 *  shipped upscaling case covers.  Nothing crashed.  The only symptom was a column in a CSV, and a
 *  conclusion in a README that the measured factor sat ABOVE the classical curve when in fact it
 *  sits below it.
 *
 *  So the checks below do not merely evaluate the function.  They pin the CONVENTION, three ways
 *  that a single shared mistake cannot satisfy at once:
 *
 *    1. against the closed form written out longhand here, independently of the header;
 *    2. against the series expansion at small phi, which is 1 - phi^2/15 for the radius
 *       convention and 1 - 0.6 phi^2 for the generalized one, so the two are told apart by a
 *       factor of nine in the first correction term;
 *    3. against the diffusion-limited asymptote, eta -> 3/phi, which the generalized convention
 *       would give as 1/phi.
 *
 *  A regression to the old formula fails 2 and 3 outright and fails 1 everywhere.
 *
 *      g++ -O2 -Wall -Wextra -std=c++11 -I../src -o t test_upscale.cpp && ./t
 * ================================================================================================ */

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "complab3d_thiele.hh"

static int fails = 0;

static void ck(const char *what, double got, double want, double tol = 1e-10) {
    const bool ok = std::fabs(got - want) <= tol * (1.0 + std::fabs(want));
    if (!ok) ++fails;
    std::printf("  %-58s got %-15.9g want %-15.9g %s\n",
                what, got, want, ok ? "ok" : "** FAIL **");
}

static void ckTrue(const char *what, bool ok) {
    if (!ok) ++fails;
    std::printf("  %-58s %s\n", what, ok ? "ok" : "** FAIL **");
}

/* The sphere with first-order kinetics, for the modulus on the RADIUS.  Written out here
 * independently of the header so a mistake has to be made twice, differently. */
static double etaRadius(double phi) {
    return (3.0 / (phi * phi)) * (phi / std::tanh(phi) - 1.0);
}

/* The same solution for the GENERALIZED modulus, which is what the header used to contain.
 * Present only so the tests below can assert that the header is NOT this. */
static double etaGeneralized(double phi) {
    const double x = 3.0 * phi;
    return (1.0 / phi) * (1.0 / std::tanh(x) - 1.0 / x);
}

int main() {
    std::printf("### the effectiveness factor, and the Thiele convention it belongs to\n");

    std::printf("--- against the closed form, written out separately\n");
    const double phis[] = { 0.05, 0.2, 0.5, 0.924008, 1.5, 3.0, 8.0 };
    for (int i = 0; i < 7; ++i) {
        char buf[80];
        std::snprintf(buf, sizeof buf, "classicalEta(%.6g)", phis[i]);
        ck(buf, complab_upscale::classicalEta(phis[i]), etaRadius(phis[i]), 1e-9);
    }

    std::printf("--- the small-phi series tells the two conventions apart\n");
    /* At phi = 0.01 the radius convention gives 1 - 1e-4/15 = 0.99999333 and the generalized one
     * gives 1 - 6e-5 = 0.99994.  They differ in the fifth decimal, which is far outside the
     * tolerance below, so this check alone catches the regression. */
    ck("classicalEta(0.01) is 1 - phi^2/15",
       complab_upscale::classicalEta(0.01), 1.0 - 1e-4 / 15.0, 1e-9);
    ckTrue("and is NOT 1 - 0.6 phi^2 (the generalized convention)",
           std::fabs(complab_upscale::classicalEta(0.01) - (1.0 - 0.6 * 1e-4)) > 1e-6);

    std::printf("--- the closed form and the series agree across the branch point\n");
    /* The header switches to the series below phi = 1e-3.  If the two branches disagree there,
     * every sweep that crosses it has a step in its classical column. */
    ck("just below the branch point", complab_upscale::classicalEta(9.99e-4),
       1.0 - 9.99e-4 * 9.99e-4 / 15.0, 1e-9);
    ck("just above it, from the closed form", complab_upscale::classicalEta(1.01e-3),
       etaRadius(1.01e-3), 1e-7);
    ckTrue("and the two branches meet",
           std::fabs(complab_upscale::classicalEta(9.99e-4)
                   - complab_upscale::classicalEta(1.01e-3)) < 1e-8);

    std::printf("--- the diffusion-limited asymptote is 3/phi, not 1/phi\n");
    /* At large phi, coth phi -> 1 and eta -> 3/phi.  The generalized convention would give 1/phi,
     * a factor of three, which is the same factor the argument was wrong by. */
    const double big = 60.0;
    ck("classicalEta(60) approaches 3/phi", complab_upscale::classicalEta(big),
       3.0 / big, 2e-3);
    ckTrue("and is about three times the generalized value there",
           std::fabs(complab_upscale::classicalEta(big) / etaGeneralized(big) - 3.0) < 0.05);

    std::printf("--- the shipped example-20 point, which is the one in the README\n");
    /* phi = 0.924008 is what the shipped case prints.  0.947345 is the value the README, the
     * [UPSCALE] line and offline/expected_sweep.csv all now carry.  The old formula gave
     * 0.700324 here, which is what made the measured 0.823 look like an over-performance. */
    ck("classicalEta(0.924008)", complab_upscale::classicalEta(0.924008), 0.947345, 1e-5);
    ckTrue("the measured 0.823 sits BELOW it, as the gate implies",
           0.823348 < complab_upscale::classicalEta(0.924008));

    std::printf("--- the guards\n");
    ck("phi = 0 gives 1",   complab_upscale::classicalEta(0.0),  1.0);
    ck("phi < 0 gives 1",   complab_upscale::classicalEta(-1.0), 1.0);
    ckTrue("eta is between 0 and 1 over four decades of phi", [] {
        for (double p = 1e-4; p < 1e2; p *= 1.2) {
            const double e = complab_upscale::classicalEta(p);
            if (!(e > 0.0 && e <= 1.0 + 1e-12)) return false;
        }
        return true;
    }());
    ckTrue("and falls monotonically with phi", [] {
        double prev = 2.0;
        for (double p = 1e-3; p < 5e1; p *= 1.1) {
            const double e = complab_upscale::classicalEta(p);
            if (e > prev + 1e-12) return false;
            prev = e;
        }
        return true;
    }());

    if (fails) { std::printf("\n%d upscaling check(s) FAILED\n", fails); return 1; }
    std::printf("\nall upscaling checks passed\n");
    return 0;
}
