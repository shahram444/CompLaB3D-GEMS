/* ================================================================================================
 *  test_abiotic.cpp  --  the abiotic sweep, run on a stub lattice
 *
 *  The one thing this has to prove is the thing the biotic processors get wrong for abiotic
 *  chemistry: that a reaction happens in a voxel holding NO biomass.  Everything else here is the
 *  same set of checks the biotic tests make, because the two paths share a budget and a deposit
 *  and must not drift apart.
 *
 *      g++ -O2 -Wall -Wextra -std=c++11 -I. -I../src -o t test_abiotic.cpp && ./t
 * ================================================================================================ */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "palabos_stub.hh"
namespace plb { PlbOut pcout; }
using namespace plb;

#include "complab3d_processors_abiotic_learned.hh"

static int fails = 0;
static void ck(const char *what, double got, double want, double tol = 1e-12) {
    const bool ok = std::fabs(got - want) <= tol * (1.0 + std::fabs(want));
    if (!ok) ++fails;
    std::printf("  %-54s got %-15.8g want %-15.8g %s\n", what, got, want, ok ? "ok" : "** FAIL **");
}
static void ckTrue(const char *what, bool ok) {
    if (!ok) ++fails;
    std::printf("  %-54s %s\n", what, ok ? "ok" : "** FAIL **");
}
static void writeFile(const char *path, const char *body) {
    std::FILE *f = std::fopen(path, "w"); std::fputs(body, f); std::fclose(f);
}

/* Sulfide reacting with dissolved iron, which is abiotic and has nothing to do with any organism:
 *      Fe(II) + HS(-)  =  FeS + H(+)
 * Second order, so both concentrations matter and a mis-bound name would change the answer. */
static const char *ABIO =
    "# abiotic iron sulfide formation, second order\n"
    "version 1\n"
    "units per_second\n"
    "vars Fe HS\n"
    "rate Fe = -12.0 * Fe * HS\n"
    "rate HS = -12.0 * Fe * HS\n";

/* substrates ... , dC ... , mask   -- the layout run_abiotic_kinetics uses */
struct Rig {
    static const int subsNum = 2;
    std::vector< BlockLattice3D<double, D3Q7>* > L;
    MetabolicConfig cfg;
    Rig() {
        for (int i = 0; i < 2*subsNum + 1; ++i)
            L.push_back(new BlockLattice3D<double, D3Q7>());
        cfg.fix_concentration.assign(subsNum, false);
        cfg.useTotals = false;
        L[4]->c.rho = 0.0;                      /* fluid */
    }
    ~Rig() { for (size_t i = 0; i < L.size(); ++i) delete L[i]; }
    void set(double Fe, double HS) { L[0]->c.rho = Fe; L[1]->c.rho = HS; }
    double dC(int i) const { return L[2 + i]->deposited(); }

    void runSym(double dt) {
        run_symbolic_abiotic3D<double, D3Q7> p(3, subsNum, dt, &cfg, 1, 2);
        Box3D d; d.x0 = d.x1 = 1; d.y0 = d.y1 = 0; d.z0 = d.z1 = 0;
        p.process(d, L);
    }
    void runGnn(double dt) {
        run_graphnet_abiotic3D<double, D3Q7> p(3, subsNum, dt, &cfg, 1, 2);
        Box3D d; d.x0 = d.x1 = 1; d.y0 = d.y1 = 0; d.z0 = d.z1 = 0;
        p.process(d, L);
    }
};

static void bindSym(complab_sym::Program &P) {
    complab_sym::Runtime &R = complab_sym::runtime();
    R.abiotic = complab_sym::Binding();
    R.abiotic.prog = &P;
    R.abiotic.subsOfVar.assign(2, -1);
    R.abiotic.subsOfVar[0] = 0; R.abiotic.subsOfVar[1] = 1;
    R.abiotic.subsOfRate.assign(2, -1);
    R.abiotic.subsOfRate[0] = 0; R.abiotic.subsOfRate[1] = 1;
    R.evaluations = R.clamped = 0;
}

int main() {
    complab_sym::Program P;
    std::string err;
    writeFile("ab.sym", ABIO);

    std::printf("--- THE POINT OF THIS FILE: no biomass anywhere, and it still reacts\n");
    ckTrue("the abiotic law loads", complab_sym::load(P, "ab.sym", &err));
    if (!err.empty()) std::printf("     %s\n", err.c_str());
    {
        bindSym(P);
        /* dt small enough that the step asks for less than the voxel holds; the case where it
         * asks for more is checked further down, where the clamp is the thing being tested */
        Rig r; r.set(2.0, 3.0); r.runSym(0.01);
        /* there is no biomass lattice on this path at all -- the vector has none */
        ck("Fe reacts at -12*Fe*HS*dt", r.dC(0), -12.0 * 2.0 * 3.0 * 0.01);
        ck("HS reacts with it",          r.dC(1), -12.0 * 2.0 * 3.0 * 0.01);
        ckTrue("and the law was evaluated", complab_sym::runtime().evaluations == 1);
    }

    std::printf("--- the rate really does depend on both, so a mis-bound name would show\n");
    {
        bindSym(P);
        Rig r; r.set(2.0, 6.0); r.runSym(0.01);
        ck("doubling HS doubles the rate", r.dC(0), -12.0 * 2.0 * 6.0 * 0.01);
    }
    {
        bindSym(P);
        Rig r; r.set(2.0, 0.0); r.runSym(0.01);
        ck("no sulfide, no reaction", r.dC(0), 0.0);
    }

    std::printf("--- the budget still cannot take a concentration below zero\n");
    {
        bindSym(P);
        Rig r; r.set(0.5, 3.0); r.runSym(1.0);      /* asks for 18, only 0.5 of Fe is present */
        ck("the draw stops at what the voxel holds", r.dC(0), -0.5);
        ck("and the other one stops at its own",     r.dC(1), -3.0);
    }

    std::printf("--- solids and walls are still skipped\n");
    {
        bindSym(P);
        Rig r; r.set(2.0, 3.0); r.L[4]->c.rho = 1.0; r.runSym(0.01);
        ck("a solid voxel does not react", r.dC(0), 0.0);
    }
    {
        bindSym(P);
        Rig r; r.set(2.0, 3.0); r.L[4]->c.rho = 2.0; r.runSym(0.01);
        ck("a bounce-back voxel does not react", r.dC(0), 0.0);
    }

    std::printf("--- a fixed substrate is held\n");
    {
        bindSym(P);
        Rig r; r.cfg.fix_concentration[0] = true;
        r.set(2.0, 3.0); r.runSym(0.01);
        ck("Fe is held", r.dC(0), 0.0);
        ck("HS is not",  r.dC(1), -0.72);
    }

    std::printf("--- with useTotals the draw is charged once, not twice\n");
    {
        bindSym(P);
        Rig r; r.cfg.useTotals = true;
        r.set(2.0, 3.0); r.runSym(0.01);
        ck("charged once", r.dC(0), -0.72);
    }

    std::printf("--- the D3Q7 deposit conserves what it was handed\n");
    {
        bindSym(P);
        Rig r; r.set(2.0, 3.0); r.runSym(0.01);
        Array<double,7> g; r.L[2]->c.getPopulations(g);
        ck("rest is 1/4",        g[0], -0.72 / 4.0);
        ck("direction 1 is 1/8", g[1], -0.72 / 8.0);
        double s = 0; for (int i = 0; i < 7; ++i) s += g[i];
        ck("the seven sum to the whole", s, -0.72);
    }

    std::printf("--- nothing loaded means nothing happens, rather than a crash\n");
    {
        complab_sym::runtime().abiotic = complab_sym::Binding();
        Rig r; r.set(2.0, 3.0); r.runSym(0.01);
        ck("no law, no increment", r.dC(0), 0.0);
    }

    std::printf("--- the graph network takes the same sweep\n");
    {
        /* identity scaling, all weights zero, so yoffset IS the rate -- same trick as the
         * biotic processor test, and no growth output, which the abiotic path requires */
        const char *NET =
            "version 1\nunits per_second\nspecies Fe HS\nreactions R1\n"
            "rounds 1\nwidth 1\ngrowth 0\nstoich 2 1\n-1\n-1\n"
            "da 1\nxoffset 0 0\nxgain 1 1\nxymin 0\n"
            "yoffset -4 -4\nygain 1 1\nyymin 0\n"
            "trainmin 0 0\ntrainmax 100 100\n"
            "Wenc 0\nbenc 0\nlayer 0\nWsr 0\nWda 0\nbR 0\nWss 0\nWrs 0\nbS 0\n"
            "Wout 0\nbout 0\n";
        writeFile("ab.gnn", NET);
        complab_gnn::Network N;
        ckTrue("the abiotic network loads", complab_gnn::load(N, "ab.gnn", &err));
        if (!err.empty()) std::printf("     %s\n", err.c_str());
        ckTrue("and it has no growth output", !N.hasGrowth);

        complab_gnn::Runtime &R = complab_gnn::runtime();
        R.abiotic = complab_gnn::Binding();
        R.abiotic.net = &N;
        R.abiotic.growthSlot = -1;
        R.abiotic.subsOfSpecies.assign(2, -1);
        R.abiotic.subsOfSpecies[0] = 0; R.abiotic.subsOfSpecies[1] = 1;
        R.evaluations = R.clamped = 0;

        Rig r; r.set(20.0, 30.0); r.runGnn(0.5);
        ck("Fe reacts with no biomass present", r.dC(0), -4.0 * 0.5);
        ck("HS with it",                        r.dC(1), -4.0 * 0.5);
        ckTrue("the network was evaluated", complab_gnn::runtime().evaluations == 1);

        Rig r2; r2.set(1.0, 30.0); r2.runGnn(0.5);
        ck("and its budget clamps too", r2.dC(0), -1.0);

        complab_gnn::runtime().abiotic = complab_gnn::Binding();
        Rig r3; r3.set(20.0, 30.0); r3.runGnn(0.5);
        ck("no network, no increment", r3.dC(0), 0.0);
    }

    std::remove("ab.sym"); std::remove("ab.gnn");
    std::printf("\n%s\n", fails ? "SOME CHECKS FAILED" : "all abiotic checks passed");
    return fails;
}
