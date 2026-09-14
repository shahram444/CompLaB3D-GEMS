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

/* the substrate list for the two-organism checks: everything either file names */
static std::vector<std::string> subsAB() {
    std::vector<std::string> v;
    v.push_back("CH4"); v.push_back("SO4"); v.push_back("H2");
    return v;
}

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


    /* ============================================================================================
     *  AN ABIOTIC LAW: one extent expression, and the stoichiometry does the rest
     * ========================================================================================== */
    std::printf("\nthe abiotic form: one extent line\n");

    /* the shipped abiotic file's shape: the same product typed once per species */
    write("./sym_ahand.sym",
      "units per_second\nvars Fe HS\n"
      "rate Fe = -1.6e2 * Fe * HS\n"
      "rate HS = -1.6e2 * Fe * HS\n");

    /* the same chemistry, declared, and with the product written too */
    write("./sym_adecl.sym",
      "units per_second\nvars Fe HS FeS\n"
      "reaction Fe -1  HS -1  FeS 1\n"
      "rate extent = 1.6e2 * Fe * HS\n");

    Program AH, AD; e.clear();
    ok(load(AH, "./sym_ahand.sym", &e), "the hand-written abiotic file loads");
    e.clear();
    ok(load(AD, "./sym_adecl.sym", &e), (std::string("the declared abiotic file loads ")+e).c_str());
    ok(AD.nDerived == 3, "three species lines were derived from the reaction");

    /* the two consumption rates must match the hand-written ones everywhere */
    double aworst = 0.0;
    for (int i = 1; i <= 40; ++i) {
        std::vector<double> vh(2), vd(3);
        vh[0] = 1e-3 * i / 40.0; vh[1] = 1e-3 * (41 - i) / 40.0;
        vd[0] = vh[0]; vd[1] = vh[1]; vd[2] = 0.0;
        std::vector<double> rh, rd;
        evaluate(AH, vh, rh, false, 0);
        evaluate(AD, vd, rd, false, 0);
        /* AH: [Fe, HS].  AD: [extent, Fe, HS, FeS] */
        for (int k = 0; k < 2; ++k) {
            const double d = std::fabs(rh[(size_t) k] - rd[(size_t) k + 1]);
            const double sc = std::fabs(rh[(size_t) k]) > 1e-300 ? std::fabs(rh[(size_t) k]) : 1.0;
            if (d / sc > aworst) aworst = d / sc;
        }
    }
    ok(aworst == 0.0, "abiotic: derived Fe and HS match the hand-written pair exactly");

    {
        std::vector<double> vd(3); vd[0] = 4e-4; vd[1] = 7e-4; vd[2] = 0.0;
        std::vector<double> rd; evaluate(AD, vd, rd, false, 0);
        ok(std::fabs(rd[1] - rd[2]) == 0.0, "Fe and HS consumed at exactly the same rate");
        ok(std::fabs(rd[3] + rd[1]) == 0.0, "FeS produced at exactly the rate Fe is consumed");
        ok(rd[3] > 0.0 && rd[1] < 0.0,      "the product rises and the reactant falls");
    }

    /* no biomass anywhere in an abiotic law: an extent line must not carry one */
    write("./sym_a1.sym", "vars Fe HS Bug\nreaction Fe -1 HS -1\nbiomass Bug\nrate extent = Fe\n");
    ok(!load(C, "./sym_a1.sym", &e), "an extent line with a biomass line is refused");
    write("./sym_a2.sym", "vars Fe HS\nreaction Fe -1 HS -1\nyield Fe 0.4\nrate extent = Fe\n");
    ok(!load(C, "./sym_a2.sym", &e), "an extent line with a yield line is refused");
    write("./sym_a3.sym", "vars Fe HS Bug\nreaction Fe -1 HS -1\nyield Fe 0.4\nbiomass Bug\n"
                          "rate growth = Fe\nrate extent = Fe\n");
    ok(!load(C, "./sym_a3.sym", &e), "a file with BOTH a growth and an extent line is refused");
    write("./sym_a4.sym", "vars Fe HS\nrate extent = Fe * HS\n");
    ok(!load(C, "./sym_a4.sym", &e), "an extent line with no reaction is refused");
    write("./sym_a5.sym", "vars Fe HS\nreaction Fe -1 HS -1\nrate Fe = -1 * Fe\n");
    ok(!load(C, "./sym_a5.sym", &e), "a reaction with neither growth nor extent is refused");
    write("./sym_a6.sym", "vars Fe HS\nreaction Fe -1 HS -1\nrate extent = Fe\nrate Fe = -1\n");
    ok(!load(C, "./sym_a6.sym", &e), "abiotic: a species with a rate line and a reaction place is refused");

    /* the wrong tag: which sweep a file belongs to is decided by which line it carries */
    std::vector<std::string> subs;
    subs.push_back("Fe"); subs.push_back("HS"); subs.push_back("FeS");
    Binding bd;
    std::string berr = bindToSubstrates(AD, subs, "Bug", bd, false);
    ok(!berr.empty(), "an abiotic file loaded through <expressions_file> is refused");
    berr = bindToSubstrates(AD, subs, "", bd, true);
    ok(berr.empty(), (std::string("the same file through <abiotic_file> binds ")+berr).c_str());

    Program BI; e.clear();
    ok(load(BI, "./sym_decl.sym", &e), "reloading the biotic declared file");
    std::vector<std::string> subs2; subs2.push_back("acetate"); subs2.push_back("o2");
    berr = bindToSubstrates(BI, subs2, "", bd, true);
    ok(!berr.empty(), "a biotic file loaded through <abiotic_file> is refused");
    berr = bindToSubstrates(BI, subs2, "Bug", bd, false);
    ok(berr.empty(), (std::string("the same file through <expressions_file> binds ")+berr).c_str());

    /* 'extent' must not be hunted for as a substrate lattice */
    ok(bd.subsOfRate.size() == BI.rates.size(), "every rate line got a binding slot");


    /* ============================================================================================
     *  TWO ORGANISMS, EACH WITH ITS OWN RATE LAW
     *
     *  <symbolic><expressions_file> names one file, and it used to go to every organism on the
     *  symbolic path.  A biotic .sym names its organism on its vars line, so a file written for
     *  ANME could not bind to SRB and the run stopped.  Two symbolic organisms with different
     *  names could not run together at all.
     * ========================================================================================== */
    std::printf("\nan organism may name its own rate law file\n");

    write("./sym_anme.sym",
      "units per_hour\nvars CH4 SO4 ANME\n"
      "reaction CH4 -1  SO4 -1\nyield CH4 0.02\nbiomass ANME\n"
      "rate growth = 0.1 * CH4 / (0.001 + CH4)\n");
    write("./sym_srb.sym",
      "units per_hour\nvars H2 SO4 SRB\n"
      "reaction H2 -4  SO4 -1\nyield H2 0.01\nbiomass SRB\n"
      "rate growth = 0.2 * H2 / (0.0005 + H2)\n");

    std::vector<std::string> names;
    names.push_back("ANME"); names.push_back("SRB");
    std::vector<int> users;
    users.push_back(0); users.push_back(1);

    /* the old behaviour, to show what was actually broken */
    {
        std::vector<std::string> none(2), want, paths;
        std::string cerr = chooseFiles(none, "./sym_anme.sym", users, names, 2, want, paths);
        ok(cerr.empty(), "one shared file is still chosen for both organisms");
        ok(paths.size() == 1, "and it is loaded once, not twice");
        Program P; e.clear();
        ok(load(P, "./sym_anme.sym", &e), "that file loads");
        Binding b;
        ok(bindToSubstrates(P, subsAB(), "ANME", b, false).empty(),
           "it binds to the organism it was written for");
        ok(!bindToSubstrates(P, subsAB(), "SRB", b, false).empty(),
           "and REFUSES the other one, which is the problem being fixed");
    }

    /* each organism naming its own */
    {
        std::vector<std::string> per(2), want, paths;
        per[0] = "./sym_anme.sym";
        per[1] = "./sym_srb.sym";
        std::string cerr = chooseFiles(per, "", users, names, 2, want, paths);
        ok(cerr.empty(), (std::string("two organisms, two files, accepted ") + cerr).c_str());
        ok(paths.size() == 2, "two distinct files to load");
        ok(want[0] == "./sym_anme.sym" && want[1] == "./sym_srb.sym",
           "each organism gets the file it named");

        Program A2, B2; e.clear();
        ok(load(A2, want[0], &e) && load(B2, want[1], &e), "both files load");
        Binding ba, bb2;
        ok(bindToSubstrates(A2, subsAB(), "ANME", ba, false).empty(), "the first binds to ANME");
        ok(bindToSubstrates(B2, subsAB(), "SRB", bb2, false).empty(), "the second binds to SRB");
        ok(A2.nDerived == 2 && B2.nDerived == 2,
           "each one derived its own substrate lines from its own reaction");
    }

    /* one names its own, the other falls back to the shared file */
    {
        std::vector<std::string> per(2), want, paths;
        per[1] = "./sym_srb.sym";
        std::string cerr = chooseFiles(per, "./sym_anme.sym", users, names, 2, want, paths);
        ok(cerr.empty(), "a mix of one override and one fallback is accepted");
        ok(want[0] == "./sym_anme.sym" && want[1] == "./sym_srb.sym",
           "the organism with no file of its own gets the shared one");
    }

    /* two organisms naming the SAME file: loaded once, bound twice */
    {
        std::vector<std::string> per(2), want, paths;
        per[0] = per[1] = "./sym_anme.sym";
        std::string cerr = chooseFiles(per, "", users, names, 2, want, paths);
        ok(cerr.empty() && paths.size() == 1,
           "two organisms naming one file load it once");
    }

    /* every way of getting it wrong */
    {
        std::vector<std::string> per(2), want, paths;
        std::vector<int> only0; only0.push_back(0);
        per[1] = "./sym_srb.sym";
        ok(!chooseFiles(per, "./sym_anme.sym", only0, names, 2, want, paths).empty(),
           "a file named against an organism that is not symbolic is refused");

        std::vector<std::string> per2(2);
        ok(!chooseFiles(per2, "", users, names, 2, want, paths).empty(),
           "a symbolic organism with no file anywhere is refused");

        std::vector<std::string> per3(5);
        per3[4] = "./sym_srb.sym";
        ok(!chooseFiles(per3, "./sym_anme.sym", users, names, 2, want, paths).empty(),
           "a file named against an organism that does not exist is refused");
    }

    /* the old form must still load untouched */
    Program E;
    ok(load(E, "./sym_hand.sym", &e) && E.nDerived == 0, "a file with no reaction block is unchanged");

    std::printf("\n%s\n", fails ? "FAILED" : "all checks passed");
    return fails;
}
