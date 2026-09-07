/* ================================================================================================
 *  test_wiring.cpp  --  the tests that would have caught a rate path that is implemented,
 *                       unit-tested, and never called.
 *
 *  Every other test in this directory constructs a processor directly and checks that it computes
 *  the right thing.  None of them can tell whether anything in the program ever constructs it.
 *  Two rate paths -- symbolic and graphnet -- were implemented, tested and shipped in exactly that
 *  state: correct, and dead.  These are the checks that close that gap.
 *
 *      1. every spelling of every rate path parses to a distinct enum value
 *      2. the enum predicates agree with each other
 *      3. a .sym program binds to substrate and microbe names, and REFUSES a mismatch
 *      4. a .gnn network binds to species names, and REFUSES a mismatch
 *      5. an abiotic file with a growth line is refused, because there is no biomass for it
 *      6. registering a binding makes it visible to the processor's runtime
 *
 *  The seventh check -- that complab.cpp actually constructs each processor -- cannot be done
 *  from C++ without Palabos, so it lives in tests/check_repo.sh as a grep over the source.
 *
 *      g++ -O2 -Wall -Wextra -std=c++11 -I. -I../src -o t test_wiring.cpp && ./t
 * ================================================================================================ */

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "palabos_stub.hh"
namespace plb { PlbOut pcout; }
using namespace plb;

#include "complab3d_metabolic.hh"
#include "complab3d_symbolic.hh"
#include "complab3d_graphnet.hh"

static int fails = 0;

static void ck(const char *what, bool ok)
{
    if (!ok) ++fails;
    std::printf("  %-62s %s\n", what, ok ? "ok" : "** FAIL **");
}
static void ckEq(const char *what, long got, long want)
{
    const bool ok = (got == want);
    if (!ok) ++fails;
    std::printf("  %-62s got %-6ld want %-6ld %s\n", what, got, want, ok ? "ok" : "** FAIL **");
}
static void writeFile(const char *path, const char *body)
{
    std::FILE *f = std::fopen(path, "w");
    std::fputs(body, f);
    std::fclose(f);
}

/* ------------------------------------------------------------------------------------------ */
int main()
{
    std::printf("\ntest_wiring -- the rate paths are reachable, not just implemented\n\n");

    /* ---- 1. every spelling parses ---------------------------------------------------------- */
    std::printf("1. reaction_type spellings\n");
    ckEq("\"kinetics\"",                rxntype::parse("kinetics"),               rxntype::KINETICS);
    ckEq("\"glpk\"",                    rxntype::parse("glpk"),                   rxntype::GLPK);
    ckEq("\"surrogate\"",               rxntype::parse("surrogate"),              rxntype::SURROGATE);
    ckEq("\"cobrapy\"",                 rxntype::parse("cobrapy"),                rxntype::COBRAPY);
    ckEq("\"symbolic\"",                rxntype::parse("symbolic"),               rxntype::SYMBOLIC);
    ckEq("\"sym\"",                     rxntype::parse("sym"),                    rxntype::SYMBOLIC);
    ckEq("\"graphnet\"",                rxntype::parse("graphnet"),               rxntype::GRAPHNET);
    ckEq("\"gnn\"",                     rxntype::parse("gnn"),                    rxntype::GRAPHNET);
    ckEq("\"graph_network\"",           rxntype::parse("graph_network"),          rxntype::GRAPHNET);
    ckEq("\"symbolic_and_kinetics\"",   rxntype::parse("symbolic_and_kinetics"),  rxntype::SYM_KNS);
    ckEq("\"graphnet_and_kinetics\"",   rxntype::parse("graphnet_and_kinetics"),  rxntype::GNN_KNS);
    ckEq("\"SYMBOLIC\" (case)",         rxntype::parse("SYMBOLIC"),               rxntype::SYMBOLIC);
    ckEq("numeric \"8\"",               rxntype::parse("8"),                      rxntype::SYMBOLIC);
    ckEq("numeric \"9\"",               rxntype::parse("9"),                      rxntype::GRAPHNET);
    ckEq("nonsense is rejected",        rxntype::parse("gnnn"),                   -1);

    /* Every enum value must have a name, and no two paths may collide. */
    ck("every enum value names itself",
       rxntype::name(rxntype::SYMBOLIC) == "symbolic" &&
       rxntype::name(rxntype::GRAPHNET) == "graphnet" &&
       rxntype::name(rxntype::SYM_KNS)  == "symbolic_and_kinetics" &&
       rxntype::name(rxntype::GNN_KNS)  == "graphnet_and_kinetics");

    /* ---- 2. the predicates agree ----------------------------------------------------------- */
    std::printf("\n2. the enum predicates\n");
    ck("usesSymbolic(SYMBOLIC) and usesSymbolic(SYM_KNS)",
       rxntype::usesSymbolic(rxntype::SYMBOLIC) && rxntype::usesSymbolic(rxntype::SYM_KNS));
    ck("usesGraphnet(GRAPHNET) and usesGraphnet(GNN_KNS)",
       rxntype::usesGraphnet(rxntype::GRAPHNET) && rxntype::usesGraphnet(rxntype::GNN_KNS));
    ck("the combined forms also use kinetics",
       rxntype::usesKinetics(rxntype::SYM_KNS) && rxntype::usesKinetics(rxntype::GNN_KNS));
    ck("a learned path is not an FBA path",
       !rxntype::usesFBA(rxntype::SYMBOLIC) && !rxntype::usesFBA(rxntype::GRAPHNET));
    ck("a learned path is not the surrogate",
       !rxntype::usesSurrogate(rxntype::SYMBOLIC) && !rxntype::usesSurrogate(rxntype::GRAPHNET));
    ck("usesLearned covers both and nothing else",
       rxntype::usesLearned(rxntype::SYMBOLIC) && rxntype::usesLearned(rxntype::GRAPHNET) &&
       !rxntype::usesLearned(rxntype::SURROGATE) && !rxntype::usesLearned(rxntype::GLPK));
    /* The one that matters most: no path may be silently swallowed by another's predicate,
     * because that is how a dispatch block ends up gated on the wrong count. */
    {
        bool overlap = false;
        const plint all[8] = { rxntype::KINETICS, rxntype::GLPK, rxntype::SURROGATE,
                               rxntype::COBRAPY, rxntype::SYMBOLIC, rxntype::GRAPHNET,
                               rxntype::SYM_KNS, rxntype::GNN_KNS };
        for (int i = 0; i < 8; ++i) {
            int n = 0;
            if (rxntype::usesGlpk     (all[i])) ++n;
            if (rxntype::usesCobrapy  (all[i])) ++n;
            if (rxntype::usesSurrogate(all[i])) ++n;
            if (rxntype::usesSymbolic (all[i])) ++n;
            if (rxntype::usesGraphnet (all[i])) ++n;
            if (n > 1) overlap = true;
        }
        ck("no reaction_type belongs to two rate paths at once", !overlap);
    }

    /* ---- 3. binding a .sym program --------------------------------------------------------- */
    std::printf("\n3. binding a symbolic program\n");
    writeFile("_w.sym",
              "version 1\n"
              "units per_hour\n"
              "vars S P Bug\n"
              "rate growth = 0.5 * S\n"
              "rate S = -2.0 * growth * Bug\n"
              "rate P = 1.0 * growth * Bug\n"
              "range S 1e-6 1e-2\n");

    complab_sym::Program P;
    std::string err;
    ck("the file loads", complab_sym::load(P, "_w.sym", &err));

    std::vector<std::string> subs;
    subs.push_back("S"); subs.push_back("P");

    complab_sym::Binding b;
    std::string berr = complab_sym::bindToSubstrates(P, subs, "Bug", b, false);
    ck("it binds against the right names", berr.empty());
    ckEq("variable S -> substrate 0", b.subsOfVar[0], 0);
    ckEq("variable P -> substrate 1", b.subsOfVar[1], 1);
    ckEq("variable Bug -> biomass (-1)", b.subsOfVar[2], -1);
    ckEq("rate growth -> biomass (-1)", b.subsOfRate[0], -1);
    ckEq("rate S -> substrate 0", b.subsOfRate[1], 0);
    ckEq("rate P -> substrate 1", b.subsOfRate[2], 1);

    /* The mismatch is the whole point of binding by name.  A file written against a different
     * <name_of_substrates> must be refused, not bound to whatever happens to sit at that index. */
    std::vector<std::string> wrong;
    wrong.push_back("acetate"); wrong.push_back("o2");
    complab_sym::Binding bw;
    berr = complab_sym::bindToSubstrates(P, wrong, "Bug", bw, false);
    ck("a substrate-name mismatch is refused", !berr.empty());
    ck("the message names the offending variable",
       berr.find("'S'") != std::string::npos || berr.find("'P'") != std::string::npos);

    complab_sym::Binding bm;
    berr = complab_sym::bindToSubstrates(P, subs, "SomeOtherBug", bm, false);
    ck("a microbe-name mismatch is refused", !berr.empty());

    /* ---- 4. binding a .gnn network --------------------------------------------------------- */
    std::printf("\n4. binding a graph network\n");
    /* Written in the loader's own layout: every matrix on one line, no dimension counts.
     * nS = 2, nR = 1, width = 2, one round, no growth column. */
    writeFile("_w.gnn",
              "version 1\n"
              "units per_hour\n"
              "species A B\n"
              "reactions R1\n"
              "rounds 1\n"
              "width 2\n"
              "growth 0\n"
              "stoich 2 1\n"
              "-1\n"
              "1\n"
              "da 1\n"
              "xoffset 0 0\n"
              "xgain 1 1\n"
              "xymin -1\n"
              "yoffset 0 0\n"
              "ygain 1 1\n"
              "yymin -1\n"
              "trainmin 0 0\n"
              "trainmax 1 1\n"
              "Wenc 1 1\n"
              "benc 0 0\n"
              "layer 0\n"
              "Wsr 1 0 0 1\n"
              "Wda 0 0\n"
              "bR 0 0\n"
              "Wss 1 0 0 1\n"
              "Wrs 1 0 0 1\n"
              "bS 0 0\n"
              "Wout 1 0\n"
              "bout 0\n");

    complab_gnn::Network N;
    std::string gerr;
    const bool gload = complab_gnn::load(N, "_w.gnn", &gerr);
    ck("the network loads", gload);
    if (gload) {
        std::vector<std::string> gsubs;
        gsubs.push_back("A"); gsubs.push_back("B");
        complab_gnn::Binding gb;
        std::string gberr = complab_gnn::bindToSubstrates(N, gsubs, gb, false);
        ck("it binds against the right species names", gberr.empty());
        ckEq("species A -> substrate 0", gb.subsOfSpecies[0], 0);
        ckEq("species B -> substrate 1", gb.subsOfSpecies[1], 1);
        ckEq("no growth column means growthSlot -1", gb.growthSlot, -1);

        std::vector<std::string> gwrong;
        gwrong.push_back("CH4"); gwrong.push_back("SO4");
        complab_gnn::Binding gbw;
        gberr = complab_gnn::bindToSubstrates(N, gwrong, gbw, false);
        ck("a species-name mismatch is refused", !gberr.empty());
        ck("the message names the offending species", gberr.find("'A'") != std::string::npos);

        /* ---- 6. the registry ---------------------------------------------------------------- */
        std::printf("\n5. the runtime registry\n");
        complab_gnn::registerNetwork(1, gb, 3);
        ck("registering for one microbe leaves the others empty",
           complab_gnn::runtime().byMicrobe.size() == 3 &&
           complab_gnn::runtime().byMicrobe[1].net != 0 &&
           complab_gnn::runtime().byMicrobe[0].net == 0);
        ck("no abiotic network until one is registered", !complab_gnn::haveAbiotic());
        complab_gnn::registerAbiotic(gb);
        ck("registering an abiotic network makes haveAbiotic() true", complab_gnn::haveAbiotic());
    }

    /* ---- 5. an abiotic file with a growth line is refused ---------------------------------- */
    std::printf("\n6. abiotic files are held to a different rule\n");
    complab_sym::Binding ba;
    berr = complab_sym::bindToSubstrates(P, subs, std::string(), ba, true);
    ck("an abiotic .sym carrying an organism is refused", !berr.empty());
    /* This particular file names the microbe as a variable, so it is refused one step
     * earlier than the growth line, and the message says which variable. Either refusal is
     * correct; what must not happen is silent acceptance. */
    ck("and the message names what it could not resolve",
       berr.find("Bug") != std::string::npos || berr.find("growth") != std::string::npos);

    /* The growth line on its own, with no biomass variable, is the other refusal. */
    writeFile("_wg.sym",
              "version 1\n"
              "units per_second\n"
              "vars S\n"
              "rate growth = 0.5 * S\n");
    complab_sym::Program PG;
    if (complab_sym::load(PG, "_wg.sym", &err)) {
        complab_sym::Binding bg;
        const std::string gerr2 = complab_sym::bindToSubstrates(PG, subs, std::string(), bg, true);
        ck("a growth line in an abiotic file is refused", !gerr2.empty());
        ck("and that message says growth", gerr2.find("growth") != std::string::npos);
    }
    std::remove("_wg.sym");

    writeFile("_wa.sym",
              "version 1\n"
              "units per_second\n"
              "vars S\n"
              "rate S = -0.01 * S\n"
              "range S 0 1\n");
    complab_sym::Program PA;
    ck("an abiotic file without growth loads", complab_sym::load(PA, "_wa.sym", &err));
    complab_sym::Binding ba2;
    berr = complab_sym::bindToSubstrates(PA, subs, std::string(), ba2, true);
    ck("and binds", berr.empty());
    ck("no abiotic law until one is registered", !complab_sym::haveAbiotic());
    complab_sym::registerAbiotic(ba2);
    ck("registering makes haveAbiotic() true", complab_sym::haveAbiotic());

    std::remove("_w.sym"); std::remove("_w.gnn"); std::remove("_wa.sym");

    std::printf("\n%s\n\n", fails ? "SOME CHECKS FAILED" : "all checks passed");
    return fails ? 1 : 0;
}
