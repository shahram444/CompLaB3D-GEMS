/* ================================================================================================
 *  test_graphnet_processor.cpp  --  run_graphnet3D actually RUN, on a stub lattice
 *
 *  test_gnn.cpp checks the network.  This checks the thing around the network: that species names
 *  reach the right substrate, that the shared budget cannot take a concentration below zero, that
 *  a product IS written (the excretion clamp is deliberately absent on this path), that growth
 *  scales with biomass, and that the D3Q7 deposit conserves what it was given.
 *
 *  Those are the five things that would be wrong in a way a simulation would not announce.  A run
 *  with a sign error still produces smooth fields and plausible pictures.
 *
 *      g++ -O2 -Wall -Wextra -std=c++11 -I<stub> -I../src -o t test_graphnet_processor.cpp && ./t
 *
 *  The stub is tests/palabos_stub.hh: just enough of Palabos to instantiate and call the
 *  processor, with cells that hold a value so the arithmetic is real.  It is NOT a Palabos
 *  substitute and is only ever used by this test.
 * ================================================================================================ */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "palabos_stub.hh"
namespace plb { PlbOut pcout; }
using namespace plb;

#include "complab3d_processors_graphnet.hh"

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

/* ------------------------------------------------------------------------------------------------
 *  A network contrived so its answer is known: with all the tanh machinery in place we cannot ask
 *  for a specific rate, so instead we set the OUTPUT scaling to put the answer where we want it.
 *  yoffset and ygain are free parameters of the file, so choosing them fixes the rate exactly.
 *
 *      species  S  P        S is consumed, P is produced
 *      rate S = -R,  rate P = +R,  growth = g
 * ---------------------------------------------------------------------------------------------- */
static std::string tinyNet(double rateS, double rateP, double growth)
{
    /* With Wenc = 0 and every bias 0, every hidden unit is tanh(0) = 0 through both rounds, so the
     * readout is exactly bout[0] = 0 for species and bout[1] = 0 for growth.  The scaled output is
     * therefore 0 for everything, and v = (0 - yymin)/ygain + yoffset = yoffset when yymin = 0.
     * So yoffset IS the rate.  Nothing is hidden and nothing is fitted. */
    char b[2048];
    std::sprintf(b,
        "version 1\n"
        "units per_second\n"
        "species S P\n"
        "reactions R1\n"
        "rounds 1\n"
        "width 1\n"
        "growth 1\n"
        "stoich 2 1\n"
        "-1\n"
        "1\n"
        "da 1\n"
        "xoffset 0 0\n"
        "xgain 1 1\n"
        "xymin 0\n"
        "yoffset %.17g %.17g %.17g\n"
        "ygain 1 1 1\n"
        "yymin 0\n"
        "trainmin 0 0\n"
        "trainmax 1000 1000\n"
        "Wenc 0\n"
        "benc 0\n"
        "layer 0\n"
        "Wsr 0\nWda 0\nbR 0\nWss 0\nWrs 0\nbS 0\n"
        "Wout 0 0\n"
        "bout 0 0\n", rateS, rateP, growth);
    return std::string(b);
}

static void writeFile(const char *path, const std::string &body) {
    std::FILE *f = std::fopen(path, "w"); std::fputs(body.c_str(), f); std::fclose(f);
}

/* ------------------------------------------------------------------------------------------------
 *  One voxel, two substrates, one organism.  The lattice vector layout is the one complab.cpp
 *  builds: substrates, biomass, dC, dB, mask.
 * ---------------------------------------------------------------------------------------------- */
struct Rig {
    static const int subsNum = 2, bioNum = 1;
    std::vector< BlockLattice3D<double, D3Q7>* > L;
    MetabolicConfig cfg;

    Rig() {
        for (int i = 0; i < 2*(subsNum + bioNum) + 1; ++i)
            L.push_back(new BlockLattice3D<double, D3Q7>());
        cfg.fix_concentration.assign(subsNum, false);
        cfg.useTotals = false;
        L[6]->c.rho = 0.0;                     /* the mask: 0 is fluid, solid/bb are 1 and 2 */
    }
    ~Rig() { for (size_t i = 0; i < L.size(); ++i) delete L[i]; }

    void set(double S, double P, double B) { L[0]->c.rho = S; L[1]->c.rho = P; L[2]->c.rho = B; }
    double dC(int i) const { return L[3 + i]->deposited(); }
    double dB()      const { return L[5]->deposited(); }

    void run(double dt) {
        std::vector<plint> gid; gid.push_back(0);
        run_graphnet3D<double, D3Q7> proc(3, subsNum, bioNum, 1, dt, gid, &cfg, 1, 2);
        Box3D dom; dom.x0 = dom.x1 = 1; dom.y0 = dom.y1 = 0; dom.z0 = dom.z1 = 0;
        proc.process(dom, L);
    }
};

static void bind(complab_gnn::Network &N) {
    complab_gnn::Runtime &R = complab_gnn::runtime();
    R.byMicrobe.assign(1, complab_gnn::Binding());
    R.byMicrobe[0].net = &N;
    R.byMicrobe[0].subsOfSpecies.assign(2, -1);
    R.byMicrobe[0].subsOfSpecies[0] = 0;      /* network species S -> substrate 0 */
    R.byMicrobe[0].subsOfSpecies[1] = 1;      /* network species P -> substrate 1 */
    R.byMicrobe[0].growthSlot = 2;
    R.evaluations = R.clamped = 0;
}

int main() {
    complab_gnn::Network N;
    std::string err;

    std::printf("--- the network reaches the right substrate, with the right sign\n");
    writeFile("p.gnn", tinyNet(-2.0, 3.0, 0.5));      /* consume S at 2, make P at 3, grow at 0.5/s */
    ckTrue("network loads", complab_gnn::load(N, "p.gnn", &err));
    if (!err.empty()) std::printf("     %s\n", err.c_str());
    bind(N);
    {
        Rig r; r.set(10.0, 1.0, 4.0); r.run(0.1);
        ck("substrate S falls by rate * dt",              r.dC(0), -0.2);
        ck("substrate P rises by rate * dt",              r.dC(1),  0.3);
        ck("biomass grows by mu * B * dt",                r.dB(),   0.5 * 4.0 * 0.1);
        ckTrue("the evaluation was counted", complab_gnn::runtime().evaluations == 1);
    }

    std::printf("--- production is NOT clamped away (the defect on the surrogate path)\n");
    {
        bind(N);
        Rig r; r.set(10.0, 0.0, 1.0); r.run(1.0);
        ckTrue("a product is written even when the voxel holds none of it", r.dC(1) > 0.0);
        ck("and it is the full amount the network asked for", r.dC(1), 3.0);
    }

    std::printf("--- the budget cannot take a concentration below zero\n");
    {
        bind(N);
        Rig r; r.set(0.5, 0.0, 1.0); r.run(1.0);     /* asks for 2.0 of S, only 0.5 is there */
        ck("the draw is limited to what the voxel holds", r.dC(0), -0.5);
        ck("and the coupled product is unaffected",       r.dC(1),  3.0);
    }
    {
        bind(N);
        Rig r; r.set(0.0, 0.0, 1.0); r.run(1.0);
        ck("an empty voxel gives up nothing", r.dC(0), 0.0);
    }

    std::printf("--- a fixed substrate is held, as <fix_concentration> promises\n");
    {
        bind(N);
        Rig r; r.cfg.fix_concentration[0] = true;
        r.set(10.0, 0.0, 1.0); r.run(1.0);
        ck("the fixed substrate is untouched",   r.dC(0), 0.0);
        ck("the free one still reacts",          r.dC(1), 3.0);
    }


    std::printf("--- with <fba_concentration_basis>total</> the draw is applied ONCE\n");
    {
        /* This is the shape of a bug that was in this processor while it was being written:
         * accumulating into dC and THEN calling drawDown(), which writes dC itself, charges the
         * voxel twice.  With no equilibrium tableau built, drawDown() takes its fallback branch
         * and subtracts the amount once, which is what must come back out. */
        bind(N);
        Rig r; r.cfg.useTotals = true;
        r.set(10.0, 0.0, 1.0); r.run(1.0);
        ck("consumption is charged once, not twice", r.dC(0), -2.0);
        ck("and production is unaffected by the basis", r.dC(1), 3.0);

        bind(N);
        Rig r2; r2.cfg.useTotals = true;
        r2.set(0.5, 0.0, 1.0); r2.run(1.0);
        ck("and it still cannot go below zero", r2.dC(0), -0.5);
    }

    std::printf("--- the D3Q7 deposit conserves what it was handed\n");
    {
        bind(N);
        Rig r; r.set(10.0, 1.0, 4.0); r.run(0.1);
        Array<double,7> g; r.L[3]->c.getPopulations(g);
        ck("rest population is 1/4 of the increment",  g[0], -0.2 / 4.0);
        for (int i = 1; i < 7; ++i)
            ck("each of the six directions is 1/8",    g[i], -0.2 / 8.0);
        double s = 0; for (int i = 0; i < 7; ++i) s += g[i];
        ck("and the seven weights sum to the whole",   s, -0.2);
    }

    std::printf("--- nothing happens where nothing should\n");
    {
        bind(N);
        Rig r; r.set(10.0, 1.0, 0.0); r.run(1.0);     /* no biomass in this voxel */
        ck("no biomass, no reaction (S)", r.dC(0), 0.0);
        ck("no biomass, no reaction (P)", r.dC(1), 0.0);
        ckTrue("and the network was never called", complab_gnn::runtime().evaluations == 0);
    }
    {
        bind(N);
        Rig r; r.set(10.0, 1.0, 4.0); r.L[6]->c.rho = 1.0;   /* solid */
        r.run(1.0);
        ck("a solid voxel does not react", r.dC(0), 0.0);
    }
    {
        bind(N);
        Rig r; r.set(10.0, 1.0, 4.0); r.L[6]->c.rho = 2.0;   /* bounce-back */
        r.run(1.0);
        ck("a bounce-back voxel does not react", r.dC(0), 0.0);
    }

    std::printf("--- a network with no growth output changes chemistry but grows nothing\n");
    {
        std::string body = tinyNet(-2.0, 3.0, 0.0);
        /* strip the growth slot: growth 0, and the output vectors lose their third entry */
        size_t p = body.find("growth 1\n");        body.replace(p, 9, "growth 0\n");
        p = body.find("yoffset");
        size_t e = body.find('\n', p);
        body.replace(p, e - p, "yoffset -2 3");
        p = body.find("ygain");   e = body.find('\n', p); body.replace(p, e - p, "ygain 1 1");
        p = body.find("Wout");    e = body.find('\n', p); body.replace(p, e - p, "Wout 0");
        p = body.find("bout");    e = body.find('\n', p); body.replace(p, e - p, "bout 0");
        writeFile("p2.gnn", body);

        complab_gnn::Network M;
        ckTrue("the growth-free network loads", complab_gnn::load(M, "p2.gnn", &err));
        if (!err.empty()) std::printf("     %s\n", err.c_str());
        complab_gnn::Runtime &R = complab_gnn::runtime();
        R.byMicrobe.assign(1, complab_gnn::Binding());
        R.byMicrobe[0].net = &M;
        R.byMicrobe[0].subsOfSpecies.assign(2, -1);
        R.byMicrobe[0].subsOfSpecies[0] = 0;
        R.byMicrobe[0].subsOfSpecies[1] = 1;
        R.byMicrobe[0].growthSlot = -1;
        R.evaluations = R.clamped = 0;

        Rig r; r.set(10.0, 1.0, 4.0); r.run(1.0);
        ck("the chemistry still moves", r.dC(0), -2.0);
        ck("but the biomass does not",  r.dB(),   0.0);
    }

    std::printf("--- an unbound organism is skipped rather than read out of bounds\n");
    {
        complab_gnn::Runtime &R = complab_gnn::runtime();
        R.byMicrobe.assign(1, complab_gnn::Binding());   /* net == 0 */
        R.evaluations = R.clamped = 0;
        Rig r; r.set(10.0, 1.0, 4.0); r.run(1.0);
        ck("no network, no increment", r.dC(0), 0.0);
    }

    std::remove("p.gnn"); std::remove("p2.gnn");
    std::printf("\n%s\n", fails ? "SOME CHECKS FAILED" : "all processor checks passed");
    return fails;
}
