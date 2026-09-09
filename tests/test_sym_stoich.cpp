/* ================================================================================================
 *  test_sym_stoich.cpp  --  the reaction/yield/biomass block of a .sym file
 *
 *  Part of the CompLaB program.  GNU Affero General Public License v3 or later.
 *
 *  A hand-written substrate line carries a ratio that nothing checks. Written as
 *
 *      rate acetate = -2.5 * growth * Bug
 *      rate o2      = -5.0 * growth * Bug
 *
 *  the 5.0 over 2.5 IS the reaction, two oxygen per acetate, and one mistyped digit makes a run
 *  that finishes, looks smooth, and creates or destroys oxygen every step. Declaring the reaction
 *  instead makes that ratio arithmetic rather than typing.
 *
 *  What this file proves: the declared form gives bit-identical rates to the hand-written one, the
 *  ratio comes out exact, a product and a coefficient other than one are handled, every way of
 *  writing the block wrongly stops the run rather than being guessed at, and a file with no block
 *  at all still loads exactly as before.
 *
 *      g++ -std=c++11 -I../src -o t test_sym_stoich.cpp && ./t
 * ================================================================================================ */

#include "complab3d_symbolic.hh"
#include <cstdio>
#include <cmath>
using namespace complab_sym;

static int fails = 0;
static void ok(bool c, const char *what) {
    std::printf("  %-58s %s\n", what, c ? "ok" : "FAIL");
    if (!c) ++fails;
}

static void write(const char *p, const char *s) { FILE*f=fopen(p,"w"); fputs(s,f); fclose(f); }

int main() {
    std::printf("the reaction/yield/biomass block\n");

    /* the hand-written form, exactly as example 17 ships it */
    write("./sym_hand.sym",
      "version 1\nunits per_hour\nvars acetate o2 Bug\n"
      "rate growth  = 0.35 * acetate / (0.05 + acetate) * o2 / (0.01 + o2)\n"
      "rate acetate = -2.5 * growth * Bug\n"
      "rate o2      = -5.0 * growth * Bug\n");

    /* the same chemistry, declared */
    write("./sym_decl.sym",
      "version 1\nunits per_hour\nvars acetate o2 Bug\n"
      "reaction acetate -1  o2 -2\nyield acetate 0.4\nbiomass Bug\n"
      "rate growth = 0.35 * acetate / (0.05 + acetate) * o2 / (0.01 + o2)\n");

    Program A, B; std::string e;
    ok(load(A, "./sym_hand.sym", &e), "the hand-written file loads");
    ok(load(B, "./sym_decl.sym", &e), (std::string("the declared file loads: ")+e).c_str());
    ok(B.nDerived == 2, "two substrate lines were derived");
    ok(A.rates.size() == B.rates.size(), "both forms end with the same number of rate lines");

    /* they must agree everywhere, not just at one point */
    double worst = 0.0;
    for (int i = 1; i <= 40; ++i) {
        std::vector<double> v(3);
        v[0] = 5e-3 * i / 40.0; v[1] = 2e-3 * i / 40.0; v[2] = 1e-5 * i;
        std::vector<double> ra, rb;
        evaluate(A, v, ra, false, 0);
        evaluate(B, v, rb, false, 0);
        for (size_t k = 0; k < ra.size(); ++k) {
            const double d = std::fabs(ra[k] - rb[k]);
            const double s = std::fabs(ra[k]) > 1e-300 ? std::fabs(ra[k]) : 1.0;
            if (d / s > worst) worst = d / s;
        }
    }
    ok(worst == 0.0, "the two forms give bit-identical rates at 40 compositions");

    /* the ratio the whole feature exists to protect */
    std::vector<double> v(3); v[0]=2e-3; v[1]=1e-3; v[2]=1e-5;
    std::vector<double> r; evaluate(B, v, r, false, 0);
    const double ratio = r[2] / r[1];
    ok(std::fabs(ratio - 2.0) < 1e-15, "derived o2 over acetate is exactly 2");

    /* every way of getting it wrong must stop the run */
    Program C;
    write("./sym_e1.sym", "version 1\nvars a b Bug\nreaction a -1 b -2\nbiomass Bug\nrate growth = a\n");
    ok(!load(C, "./sym_e1.sym", &e), "reaction without yield is refused");
    write("./sym_e2.sym", "version 1\nvars a b Bug\nreaction a -1 b -2\nyield a 0.4\nrate growth = a\n");
    ok(!load(C, "./sym_e2.sym", &e), "reaction without biomass is refused");
    write("./sym_e3.sym", "version 1\nvars a b Bug\nreaction a -1 b -2\nyield a 0.4\nbiomass Bug\nrate a = 1\n");
    ok(!load(C, "./sym_e3.sym", &e), "reaction without a growth line is refused");
    write("./sym_e4.sym", "version 1\nvars a b Bug\nreaction a -1 b -2\nyield c 0.4\nbiomass Bug\nrate growth = a\n");
    ok(!load(C, "./sym_e4.sym", &e), "a yield quoted against a species not in the reaction is refused");
    write("./sym_e5.sym", "version 1\nvars a b Bug\nreaction a -1 b -2\nyield a 0.4\nbiomass Bug\n"
                         "rate growth = a\nrate a = -1 * growth * Bug\n");
    ok(!load(C, "./sym_e5.sym", &e), "a species with both a rate line and a reaction place is refused");
    write("./sym_e6.sym", "version 1\nvars a b Bug\nreaction a -1 b -2\nyield a -0.4\nbiomass Bug\nrate growth = a\n");
    ok(!load(C, "./sym_e6.sym", &e), "a negative yield is refused");
    write("./sym_e7.sym", "version 1\nvars a b Bug\nreaction a -1 c -2\nyield a 0.4\nbiomass Bug\nrate growth = a\n");
    ok(!load(C, "./sym_e7.sym", &e), "a reaction species not in vars is refused");
    write("./sym_e8.sym", "version 1\nvars a b Bug\nreaction a -1 b\nyield a 0.4\nbiomass Bug\nrate growth = a\n");
    ok(!load(C, "./sym_e8.sym", &e), "a species with no coefficient is refused");
    write("./sym_e9.sym", "version 1\nvars a b Bug\nyield a 0.4\nrate growth = a\n");
    ok(!load(C, "./sym_e9.sym", &e), "a yield with no reaction is refused");

    /* a product, and a yield quoted against a species with coefficient other than one */
    write("./sym_p.sym",
      "version 1\nunits per_hour\nvars lac so4 ac Bug\n"
      "reaction lac -2  so4 -1  ac 2\nyield lac 0.1\nbiomass Bug\nrate growth = lac\n");
    Program D; e.clear();
    ok(load(D, "./sym_p.sym", &e), (std::string("a three-species reaction loads ")+e).c_str());
    std::vector<double> v2(4); v2[0]=1.0; v2[1]=1.0; v2[2]=0.0; v2[3]=1.0;
    std::vector<double> r2; evaluate(D, v2, r2, false, 0);
    /* per unit biomass: lac 2/2/0.1 = 10 consumed, so4 1/2/0.1 = 5 consumed, ac 2/2/0.1 = 10 made */
    /* the file says per_hour, so every rate leaves divided by 3600 */
    const double H = 1.0 / 3600.0;
    ok(std::fabs(r2[1] + 10.0 * H) < 1e-15, "lactate coefficient is -10");
    ok(std::fabs(r2[2] +  5.0 * H) < 1e-15, "sulfate coefficient is -5");
    ok(std::fabs(r2[3] - 10.0 * H) < 1e-15, "acetate coefficient is +10, a product");
    ok(std::fabs(r2[1] / r2[2] - 2.0) < 1e-15, "lactate over sulfate is exactly 2, as written");

    /* the old form must still load untouched */
    Program E;
    ok(load(E, "./sym_hand.sym", &e) && E.nDerived == 0, "a file with no reaction block is unchanged");

    std::printf("\n%s\n", fails ? "FAILED" : "all checks passed");
    return fails;
}
