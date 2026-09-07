/* ================================================================================================
 *  test_thermo.cpp  --  complab3d_thermo.hh checked without Palabos
 *
 *  The gate multiplies every rate in the solver, so a sign error in it does not crash anything: it
 *  produces a run that looks finished and is wrong everywhere.  That is the failure this file is
 *  written against.  Every number below is either worked out by hand in the comment beside it or
 *  computed a second time from the equation written out longhand, so a shared mistake between the
 *  header and the test has to be made twice, differently.
 *
 *      g++ -O2 -Wall -Wextra -std=c++11 -I../src -o t test_thermo.cpp && ./t
 * ================================================================================================ */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "complab3d_thermo.hh"

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

static void write(const char *path, const std::string &body) {
    std::FILE *f = std::fopen(path, "w");
    std::fputs(body.c_str(), f);
    std::fclose(f);
}

/* ------------------------------------------------------------------------------------------------
 *  Anaerobic oxidation of methane, the reaction the feature was built for.
 *
 *      CH4 + SO4(2-)  ->  HS(-) + HCO3(-) + H2O
 *
 *  dG0 = -33.12 kJ/mol at 4 C (277.15 K), which is the value used by He et al. (2021).  Water is
 *  not a substrate of the run and its activity is taken as one, so it does not appear in stoich.
 *  m dG_ATP = 0.25 x 50 = 12.5 kJ/mol, inside the 10-20 kJ/mol band for anaerobic syntrophs.
 * ---------------------------------------------------------------------------------------------- */
static const char *AOM =
    "# anaerobic oxidation of methane, coupled to sulfate reduction\n"
    "provenance dG0 from He et al. 2021; threshold from Hoehler et al. 2001\n"
    "version 1\n"
    "temperature 277.15\n"
    "energy_units kJ\n"
    "\n"
    "reaction AOM\n"
    "  microbe   ANME\n"
    "  dG0       -33.12\n"
    "  stoich    CH4 -1  SO4 -1  HS 1  HCO3 1\n"
    "  atp       0.25\n"
    "  dGatp     50.0\n"
    "  chi       1.0\n"
    "  cmin      1e-9\n"
    "  yield     on\n"
    "  dGana     -30.0\n"
    "  carbon    1\n"
    "  reduction 8.0\n"
    "  yieldmax  0.2\n";

static const double RGAS = 8.314462618e-3;      /* kJ / mol / K */

/* the same four equations written out with no loops and no structs */
static double byHand(double cCH4, double cSO4, double cHS, double cHCO3, double *dGout)
{
    const double T   = 277.15;
    const double RT  = RGAS * T;
    const double lnQ = -1.0 * std::log(cCH4) - 1.0 * std::log(cSO4)
                     + 1.0 * std::log(cHS)  + 1.0 * std::log(cHCO3);
    const double dG  = -33.12 + RT * lnQ;
    if (dGout) *dGout = dG;
    const double f   = dG + 0.25 * 50.0;
    if (f >= 0.0) return 0.0;
    return 1.0 - std::exp(f / (1.0 * RT));
}

int main()
{
    std::vector<std::string> subs;
    subs.push_back("CH4"); subs.push_back("SO4"); subs.push_back("HS"); subs.push_back("HCO3");
    std::vector<std::string> mics;
    mics.push_back("ANME");

    std::printf("--- the file loads and binds to this run\n");
    write("aom.thm", AOM);
    complab_thermo::Model M;
    std::string err;
    ckTrue("loads", complab_thermo::load(M, "aom.thm", &err));
    if (!err.empty()) std::printf("     %s\n", err.c_str());
    ckTrue("one reaction block", M.specs.size() == 1);
    ckTrue("provenance kept, '#' comment dropped",
           M.provenance.find("He et al. 2021") != std::string::npos);
    ck("RT at 4 C", M.RT(), RGAS * 277.15);
    ckTrue("binds", complab_thermo::bindToRun(M, subs, mics).empty());
    ckTrue("CH4 bound to substrate 0", M.specs[0].subsOfSpecies[0] == 0);
    ckTrue("HCO3 bound to substrate 3", M.specs[0].subsOfSpecies[3] == 3);
    ckTrue("the microbe name resolved to index 0", M.specs[0].microbe == 0);

    /* -------------------------------------------------------------------------------------------
     *  The three regimes the gate exists to tell apart.  Concentrations in mol/L.
     * ----------------------------------------------------------------------------------------- */
    std::printf("--- dG and F_T in the three regimes\n");
    {
        struct Case { const char *what; double c[4]; };
        const Case cases[3] = {
            /* rim of an aggregate: substrate delivered, products swept away */
            { "rim, substrate rich",      { 2.0e-3, 2.0e-2, 1.0e-5, 1.0e-4 } },
            /* interior: the same reaction, its own products accumulated */
            { "interior, products built", { 6.0e-5, 5.0e-4, 6.0e-3, 2.0e-2 } },
            /* dead zone: past the threshold */
            { "dead zone, past threshold",{ 1.0e-6, 1.0e-5, 3.0e-2, 5.0e-2 } }
        };
        for (int i = 0; i < 3; ++i) {
            std::vector<double> c(cases[i].c, cases[i].c + 4);
            double dGa = 0.0, dGb = 0.0;
            const double got  = complab_thermo::factor(M.specs[0], M, c, &dGa);
            const double want = byHand(c[0], c[1], c[2], c[3], &dGb);
            char lbl[128];
            std::sprintf(lbl, "dG   %s", cases[i].what); ck(lbl, dGa, dGb);
            std::sprintf(lbl, "F_T  %s", cases[i].what); ck(lbl, got, want);
            std::printf("       dG = %8.3f kJ/mol   F_T = %.6f\n", dGa, got);
        }
    }

    std::printf("--- the two limits, which are the ones a sign error breaks\n");
    {
        /* Strongly exergonic: F_T must go to 1, not to 0. */
        double c1a[4] = { 1.0, 1.0, 1e-9, 1e-9 };
        std::vector<double> c1(c1a, c1a + 4);
        double dG = 0.0;
        const double F1 = complab_thermo::factor(M.specs[0], M, c1, &dG);
        ckTrue("a strongly exergonic composition gives dG < 0", dG < -50.0);
        ck("and F_T = 1 there", F1, 1.0, 1e-9);

        /* Endergonic: F_T must be exactly zero, not a small positive number. */
        double c2a[4] = { 1e-9, 1e-9, 1.0, 1.0 };
        std::vector<double> c2(c2a, c2a + 4);
        const double F2 = complab_thermo::factor(M.specs[0], M, c2, &dG);
        ckTrue("an endergonic composition gives dG > 0", dG > 0.0);
        ck("and F_T = 0 there", F2, 0.0, 0.0);

        /* Exactly at the threshold, dG = -m dG_ATP = -12.5, F_T must be exactly zero. */
        complab_thermo::Spec s = M.specs[0];
        s.dG0 = -12.5;
        double c3a[4] = { 1.0, 1.0, 1.0, 1.0 };            /* ln Q = 0, so dG = dG0 */
        std::vector<double> c3(c3a, c3a + 4);
        ck("at unit activities dG is exactly dG0", complab_thermo::deltaG(s, M, c3), -12.5);
        ck("F_T is exactly zero at the threshold", complab_thermo::factor(s, M, c3), 0.0, 0.0);
        s.dG0 = -12.5 - 1e-6;
        ckTrue("and strictly positive just past it", complab_thermo::factor(s, M, c3) > 0.0);
    }

    std::printf("--- a zero concentration floors, it does not become an infinity\n");
    {
        double ca[4] = { 0.0, 0.0, 1e-3, 1e-3 };
        std::vector<double> c(ca, ca + 4);
        const double dG = complab_thermo::deltaG(M.specs[0], M, c);
        ckTrue("dG is finite when a reactant is exactly zero", dG == dG && std::fabs(dG) < 1e6);
        /* the floor is cmin, so the answer must equal the one at cmin exactly */
        double cb[4] = { 1e-9, 1e-9, 1e-3, 1e-3 };
        std::vector<double> c2(cb, cb + 4);
        ck("and equals the value at cmin", dG, complab_thermo::deltaG(M.specs[0], M, c2));
        ck("a reactant at zero gives F_T = 0", complab_thermo::factor(M.specs[0], M, c), 0.0, 0.0);
    }

    std::printf("--- chi stretches the transition, it does not move the threshold\n");
    {
        complab_thermo::Spec a = M.specs[0], b = M.specs[0];
        b.chi = 4.0;
        double ca[4] = { 1e-3, 1e-3, 1e-3, 1e-3 };
        std::vector<double> c(ca, ca + 4);
        const double Fa = complab_thermo::factor(a, M, c);
        const double Fb = complab_thermo::factor(b, M, c);
        ckTrue("a larger chi gives a smaller F_T away from the threshold", Fb < Fa);
        /* At the threshold both must be zero whatever chi is. */
        a.dG0 = b.dG0 = -12.5;
        double cu[4] = { 1.0, 1.0, 1.0, 1.0 };
        std::vector<double> c1(cu, cu + 4);
        ck("chi = 1 is zero at the threshold", complab_thermo::factor(a, M, c1), 0.0, 0.0);
        ck("chi = 4 is zero at the same place", complab_thermo::factor(b, M, c1), 0.0, 0.0);
    }

    std::printf("--- Heijnen dissipation and the energy yield\n");
    {
        /* Eq. 3.9 written out for methane: NoC = 1, degree of reduction 8.
         *   200 + 18*(6-1)^1.8 + exp[ ((3.8-8)^2)^0.16 * (3.6+0.4*1) ]              */
        const double a = 200.0 + 18.0 * std::pow(5.0, 1.8);
        const double b = std::exp(std::pow((3.8 - 8.0) * (3.8 - 8.0), 0.16) * (3.6 + 0.4));
        ck("dG_dis for methane", complab_thermo::dissipation(1.0, 8.0), a + b);

        /* Eq. 3.10 at the rim composition. */
        double ca[4] = { 2.0e-3, 2.0e-2, 1.0e-5, 1.0e-4 };
        std::vector<double> c(ca, ca + 4);
        const double dG  = complab_thermo::deltaG(M.specs[0], M, c);
        const double Y   = complab_thermo::yieldOf(M.specs[0], dG);
        const double dis = complab_thermo::dissipation(1.0, 8.0);
        ck("yield at the rim", Y, dG / -(-30.0 + dis));
        ckTrue("the yield is positive and small", Y > 0.0 && Y < 0.2);
        std::printf("       dG_dis = %.2f kJ/C-mol   Y = %.5f C-mol biomass / mol CH4\n", dis, Y);

        /* The cap is a cap. */
        complab_thermo::Spec s = M.specs[0];
        s.yieldMax = 1e-4;
        ck("yieldmax caps the result", complab_thermo::yieldOf(s, dG), 1e-4);

        /* yield off means the caller keeps its own number. */
        s = M.specs[0]; s.yieldOn = false;
        ck("yield off returns zero", complab_thermo::yieldOf(s, dG), 0.0, 0.0);
    }

    std::printf("--- the FBA growth efficiency, Craig Eq. 3.20\n");
    {
        /* 1 % of the methane carbon in biomass: 22.5 C per biomass unit, F_BM 4.4e-4, F_CH4 1.0 */
        ck("gamma = nuC * F_BM / F_CH4",
           complab_thermo::fbaGrowthEfficiency(22.5, 4.4e-4, 1.0), 22.5 * 4.4e-4);
        ck("the sign of the uptake flux does not matter",
           complab_thermo::fbaGrowthEfficiency(22.5, 4.4e-4, -1.0), 22.5 * 4.4e-4);
        ck("a zero substrate flux gives zero, not a division by zero",
           complab_thermo::fbaGrowthEfficiency(22.5, 4.4e-4, 0.0), 0.0, 0.0);
    }

    std::printf("--- the runtime gate, which is what the processors call\n");
    {
        complab_thermo::Runtime &R = complab_thermo::runtime();
        R = complab_thermo::Runtime();
        ck("with no model registered the gate is 1", complab_thermo::gateFor(0, std::vector<double>(4, 1e-3)), 1.0);
        ckTrue("and nothing is counted", R.evaluations == 0);

        R.model = M;
        complab_thermo::registerModel(2);         /* two organisms, only the first has a block */
        ckTrue("the gate is active", complab_thermo::enabled());
        ckTrue("microbe 0 has a spec", R.specOfMicrobe[0] == 0);
        ckTrue("microbe 1 has none",   R.specOfMicrobe[1] == -1);

        double ca[4] = { 2.0e-3, 2.0e-2, 1.0e-5, 1.0e-4 };
        std::vector<double> c(ca, ca + 4);
        const double F0 = complab_thermo::gateFor(0, c);
        ck("microbe 0 is gated", F0, byHand(c[0], c[1], c[2], c[3], 0));
        ck("microbe 1 is not gated at all", complab_thermo::gateFor(1, c), 1.0);
        ckTrue("only the gated organism was counted", R.evaluations == 1);

        double cb[4] = { 1e-9, 1e-9, 1.0, 1.0 };
        std::vector<double> c2(cb, cb + 4);
        complab_thermo::gateFor(0, c2);
        ckTrue("a blocked evaluation is counted as blocked", R.blocked == 1);
        ckTrue("the report mentions the threshold",
               complab_thermo::runtimeReport().find("threshold") != std::string::npos);

        R = complab_thermo::Runtime();
        ckTrue("a run that never evaluated says nothing", complab_thermo::runtimeReport().empty());
    }

    std::printf("--- 'all' applies one block to every organism\n");
    {
        std::string body(AOM);
        body.replace(body.find("microbe   ANME"), 14, "microbe   all ");
        write("all.thm", body);
        complab_thermo::Model A;
        ckTrue("loads", complab_thermo::load(A, "all.thm", &err));
        ckTrue("binds", complab_thermo::bindToRun(A, subs, mics).empty());
        ckTrue("microbe resolves to -1", A.specs[0].microbe == -1);
        complab_thermo::Runtime &R = complab_thermo::runtime();
        R = complab_thermo::Runtime();
        R.model = A;
        complab_thermo::registerModel(3);
        ckTrue("all three organisms get the block",
               R.specOfMicrobe[0] == 0 && R.specOfMicrobe[1] == 0 && R.specOfMicrobe[2] == 0);
        R = complab_thermo::Runtime();
    }

    std::printf("--- units: writing the same file in J must give the same physics\n");
    {
        std::string body(AOM);
        body.replace(body.find("energy_units kJ"), 15, "energy_units J ");
        body.replace(body.find("dG0       -33.12"), 16, "dG0       -33120");
        body.replace(body.find("dGatp     50.0"),   14, "dGatp     50000");
        body.replace(body.find("dGana     -30.0"),  15, "dGana     -30000");
        write("aom_j.thm", body);
        complab_thermo::Model J;
        ckTrue("the joule file loads", complab_thermo::load(J, "aom_j.thm", &err));
        if (!err.empty()) std::printf("     %s\n", err.c_str());
        ckTrue("binds", complab_thermo::bindToRun(J, subs, mics).empty());
        ck("dG0 came back in kJ", J.specs[0].dG0, -33.12, 1e-9);
        double ca[4] = { 2.0e-3, 2.0e-2, 1.0e-5, 1.0e-4 };
        std::vector<double> c(ca, ca + 4);
        ck("and F_T is identical to the kJ file",
           complab_thermo::factor(J.specs[0], J, c), complab_thermo::factor(M.specs[0], M, c));
        ckTrue("the log says which unit was read",
               complab_thermo::describe(J, "aom_j.thm").find("read in J") != std::string::npos);
    }

    std::printf("--- a file that is wrong is refused, not half-loaded\n");
    {
        struct Case { const char *what; std::string body; };
        std::vector<Case> bad;
        std::string b;

        b = AOM; b.erase(b.find("  stoich    CH4 -1  SO4 -1  HS 1  HCO3 1\n"), 41);
        { Case c = { "a reaction with no stoich line", b }; bad.push_back(c); }

        b = AOM; b.replace(b.find("stoich    CH4 -1  SO4 -1  HS 1  HCO3 1"), 38,
                           "stoich    CH4 -1  SO4    HS 1  HCO3 1 ");
        { Case c = { "a stoich line with a missing coefficient", b }; bad.push_back(c); }

        b = AOM; b.replace(b.find("CH4 -1"), 6, "CH4  0");
        { Case c = { "a zero stoichiometric coefficient", b }; bad.push_back(c); }

        b = AOM; b.replace(b.find("chi       1.0"), 13, "chi       0.0");
        { Case c = { "chi of zero, which divides by zero", b }; bad.push_back(c); }

        b = AOM; b.replace(b.find("cmin      1e-9"), 14, "cmin      0e+0");
        { Case c = { "cmin of zero, which lets ln Q reach -inf", b }; bad.push_back(c); }

        b = AOM; b += "  dGloss    -5\n";
        { Case c = { "a negative transfer loss", b }; bad.push_back(c); }

        b = AOM; b += "  wibble    3\n";
        { Case c = { "an unknown keyword", b }; bad.push_back(c); }

        b = AOM; b.replace(b.find("energy_units kJ"), 15, "energy_units eV");
        { Case c = { "an unrecognised energy unit", b }; bad.push_back(c); }

        b = AOM; b.replace(b.find("temperature 277.15"), 18, "temperature -277.1");
        { Case c = { "a negative temperature", b }; bad.push_back(c); }

        b = "dG0 -1\n";
        { Case c = { "a key before any reaction block", b }; bad.push_back(c); }

        b = "provenance nothing here\n";
        { Case c = { "a file with no reaction block at all", b }; bad.push_back(c); }

        b = AOM; b += "reaction AOM\n  stoich CH4 -1 HS 1\n";
        { Case c = { "two reactions with the same name", b }; bad.push_back(c); }

        for (size_t i = 0; i < bad.size(); ++i) {
            write("bad.thm", bad[i].body);
            complab_thermo::Model B;
            std::string e;
            const bool loaded = complab_thermo::load(B, "bad.thm", &e);
            if (loaded) ++fails;
            std::printf("  %-58s %s\n", bad[i].what, loaded ? "** FAIL: accepted **" : "refused ok");
        }

        complab_thermo::Model X;
        ckTrue("a file that is not there is refused",
               !complab_thermo::load(X, "no_such_file.thm", &err));
    }

    std::printf("--- a file that does not match this run is refused at bind time\n");
    {
        complab_thermo::Model B;
        ckTrue("loads", complab_thermo::load(B, "aom.thm", &err));
        std::vector<std::string> wrong;
        wrong.push_back("glucose"); wrong.push_back("o2");
        const std::string m = complab_thermo::bindToRun(B, wrong, mics);
        ckTrue("a species that is not a substrate is refused", !m.empty());
        ckTrue("and the message names the species and lists what the run has",
               m.find("CH4") != std::string::npos && m.find("glucose") != std::string::npos);

        complab_thermo::Model C;
        complab_thermo::load(C, "aom.thm", &err);
        std::vector<std::string> nomic;
        nomic.push_back("Ecoli");
        ckTrue("a microbe name that is not in the run is refused",
               !complab_thermo::bindToRun(C, subs, nomic).empty());

        /* an index is accepted where a name would be */
        complab_thermo::Model D;
        std::string body(AOM);
        body.replace(body.find("microbe   ANME"), 14, "microbe   0   ");
        write("idx.thm", body);
        complab_thermo::load(D, "idx.thm", &err);
        ckTrue("an index binds", complab_thermo::bindToRun(D, subs, mics).empty());
        ckTrue("to the same organism", D.specs[0].microbe == 0);

        /* two blocks on one organism would multiply two gates together */
        complab_thermo::Model E;
        std::string two(AOM);
        two += "\nreaction AOM2\n  microbe ANME\n  dG0 -20\n  stoich CH4 -1 HS 1\n";
        write("two.thm", two);
        ckTrue("loads", complab_thermo::load(E, "two.thm", &err));
        ckTrue("two gates on one organism are refused",
               !complab_thermo::bindToRun(E, subs, mics).empty());
    }

    std::remove("aom.thm");  std::remove("aom_j.thm"); std::remove("all.thm");
    std::remove("bad.thm");  std::remove("idx.thm");   std::remove("two.thm");
    std::printf("\n%s\n", fails ? "SOME CHECKS FAILED" : "all thermodynamics checks passed");
    return fails;
}
