/* ================================================================================================
 *  test_symbolic_processor.cpp  --  run_symbolic3D actually RUN, on a stub lattice
 *
 *  test_sym.cpp checks the parser and the evaluator.  This checks the thing around them, and it is
 *  the same set of five questions test_graphnet_processor.cpp asks, deliberately: the two
 *  processors share a budget loop and a deposit, so if one of them is ever changed the other has to
 *  be, and the two tests failing together is how that gets noticed.
 *
 *      g++ -O2 -Wall -Wextra -std=c++11 -I. -I../src -o t test_symbolic_processor.cpp && ./t
 * ================================================================================================ */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "palabos_stub.hh"
namespace plb { PlbOut pcout; }
using namespace plb;

#include "complab3d_processors_symbolic.hh"

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

/* Substrate 0 is S, substrate 1 is P, and the organism is called X.  First-order decay of S, a
 * product P at a fixed stoichiometric ratio to it, and growth proportional to S -- so every one of
 * the three outputs depends on something, and a binding that went to the wrong substrate would
 * change the answer rather than leaving it alone. */
static const char *SYM =
    "# a rate law whose answer can be worked out on paper\n"
    "version 1\n"
    "units per_second\n"
    "vars S P X\n"
    "rate growth = 0.05 * S\n"
    "rate S = -0.2 * S\n"
    "rate P = -0.5 * S\n";

struct Rig {
    static const int subsNum = 2, bioNum = 1;
    std::vector< BlockLattice3D<double, D3Q7>* > L;
    MetabolicConfig cfg;

    Rig() {
        for (int i = 0; i < 2*(subsNum + bioNum) + 1; ++i)
            L.push_back(new BlockLattice3D<double, D3Q7>());
        cfg.fix_concentration.assign(subsNum, false);
        cfg.useTotals = false;
        L[6]->c.rho = 0.0;
    }
    ~Rig() { for (size_t i = 0; i < L.size(); ++i) delete L[i]; }

    void set(double S, double P, double B) { L[0]->c.rho = S; L[1]->c.rho = P; L[2]->c.rho = B; }
    double dC(int i) const { return L[3 + i]->deposited(); }
    double dB()      const { return L[5]->deposited(); }

    void run(double dt) {
        std::vector<plint> gid; gid.push_back(0);
        run_symbolic3D<double, D3Q7> proc(3, subsNum, bioNum, 1, dt, gid, &cfg, 1, 2);
        Box3D dom; dom.x0 = dom.x1 = 1; dom.y0 = dom.y1 = 0; dom.z0 = dom.z1 = 0;
        proc.process(dom, L);
    }
};

/* what complab.cpp's binding step would produce for the file above */
static void bind(complab_sym::Program &P) {
    complab_sym::Runtime &R = complab_sym::runtime();
    R.byMicrobe.assign(1, complab_sym::Binding());
    R.byMicrobe[0].prog = &P;
    R.byMicrobe[0].subsOfVar.assign(3, -1);
    R.byMicrobe[0].subsOfVar[0] = 0;      /* S */
    R.byMicrobe[0].subsOfVar[1] = 1;      /* P */
    R.byMicrobe[0].subsOfVar[2] = -1;     /* X is the organism's own biomass */
    R.byMicrobe[0].subsOfRate.assign(3, -1);
    R.byMicrobe[0].subsOfRate[0] = -1;    /* growth */
    R.byMicrobe[0].subsOfRate[1] = 0;     /* S */
    R.byMicrobe[0].subsOfRate[2] = 1;     /* P */
    R.evaluations = R.clamped = 0;
}

int main() {
    complab_sym::Program P;
    std::string err;

    writeFile("s.sym", SYM);
    std::printf("--- each rate reaches the substrate it names\n");
    ckTrue("the file loads", complab_sym::load(P, "s.sym", &err));
    if (!err.empty()) std::printf("     %s\n", err.c_str());
    ckTrue("three rates, three variables", P.rates.size() == 3 && P.vars.size() == 3);
    {
        bind(P);
        Rig r; r.set(10.0, 1.0, 4.0); r.run(0.1);
        ck("S falls at -0.2*S over dt",   r.dC(0), -0.2 * 10.0 * 0.1);
        ck("P falls at -0.5*S over dt",   r.dC(1), -0.5 * 10.0 * 0.1);
        ck("X grows at 0.05*S times B",   r.dB(),   0.05 * 10.0 * 4.0 * 0.1);
        ckTrue("the evaluation was counted", complab_sym::runtime().evaluations == 1);
    }

    std::printf("--- units are applied once, and chained rates stay pure\n");
    {
        std::string body(SYM);
        body.replace(body.find("per_second"), 10, "per_hour");
        writeFile("s_hr.sym", body.c_str());
        complab_sym::Program H;
        ckTrue("the per_hour file loads", complab_sym::load(H, "s_hr.sym", &err));
        bind(H);
        Rig r; r.set(10.0, 1.0, 4.0); r.run(0.1);
        ck("every rate is divided by 3600 exactly once", r.dC(0), -0.2 * 10.0 * 0.1 / 3600.0);
        ck("and the P:S ratio is untouched by the unit", r.dC(1) / r.dC(0), 2.5);
    }

    std::printf("--- production is NOT clamped away\n");
    {
        std::string body =
            "version 1\nunits per_second\nvars S P X\n"
            "rate S = -2\nrate P = 3\nrate growth = 0\n";
        writeFile("s2.sym", body.c_str());
        complab_sym::Program Q;
        ckTrue("loads", complab_sym::load(Q, "s2.sym", &err));
        complab_sym::Runtime &R = complab_sym::runtime();
        R.byMicrobe.assign(1, complab_sym::Binding());
        R.byMicrobe[0].prog = &Q;
        R.byMicrobe[0].subsOfVar.assign(3, -1);
        R.byMicrobe[0].subsOfVar[0] = 0; R.byMicrobe[0].subsOfVar[1] = 1;
        R.byMicrobe[0].subsOfRate.assign(3, -1);
        R.byMicrobe[0].subsOfRate[0] = 0; R.byMicrobe[0].subsOfRate[1] = 1;
        R.evaluations = R.clamped = 0;

        Rig r; r.set(10.0, 0.0, 1.0); r.run(1.0);
        ckTrue("a product is written into a voxel holding none of it", r.dC(1) > 0.0);
        ck("at the full amount asked for", r.dC(1), 3.0);

        Rig r2; r2.set(0.5, 0.0, 1.0); r2.run(1.0);
        ck("but consumption stops at what is there", r2.dC(0), -0.5);
    }

    std::printf("--- a fixed substrate is held\n");
    {
        bind(P);
        Rig r; r.cfg.fix_concentration[0] = true;
        r.set(10.0, 1.0, 4.0); r.run(0.1);
        ck("the fixed substrate is untouched", r.dC(0), 0.0);
        ck("the free one still reacts",        r.dC(1), -0.5);
    }


    std::printf("--- with <fba_concentration_basis>total</> the draw is applied ONCE\n");
    {
        bind(P);
        Rig r; r.cfg.useTotals = true;
        r.set(10.0, 1.0, 4.0); r.run(0.1);
        ck("consumption is charged once, not twice", r.dC(0), -0.2);
        ck("the second consumed substrate too",      r.dC(1), -0.5);
    }

    std::printf("--- the D3Q7 deposit conserves what it was handed\n");
    {
        bind(P);
        Rig r; r.set(10.0, 1.0, 4.0); r.run(0.1);
        Array<double,7> g; r.L[3]->c.getPopulations(g);
        ck("rest is 1/4",   g[0], -0.2 / 4.0);
        ck("direction 1 is 1/8", g[1], -0.2 / 8.0);
        double s = 0; for (int i = 0; i < 7; ++i) s += g[i];
        ck("the seven sum to the whole", s, -0.2);
    }

    std::printf("--- nothing happens where nothing should\n");
    {
        bind(P);
        Rig r; r.set(10.0, 1.0, 0.0); r.run(1.0);
        ck("no biomass, no reaction", r.dC(0), 0.0);
        ckTrue("and the law was never evaluated", complab_sym::runtime().evaluations == 0);
    }
    {
        bind(P);
        Rig r; r.set(10.0, 1.0, 4.0); r.L[6]->c.rho = 1.0; r.run(1.0);
        ck("a solid voxel does not react", r.dC(0), 0.0);
    }
    {
        bind(P);
        Rig r; r.set(10.0, 1.0, 4.0); r.L[6]->c.rho = 2.0; r.run(1.0);
        ck("a bounce-back voxel does not react", r.dC(0), 0.0);
    }
    {
        complab_sym::Runtime &R = complab_sym::runtime();
        R.byMicrobe.assign(1, complab_sym::Binding());   /* prog == 0 */
        R.evaluations = R.clamped = 0;
        Rig r; r.set(10.0, 1.0, 4.0); r.run(1.0);
        ck("an unbound organism is skipped", r.dC(0), 0.0);
    }

    std::printf("--- biomass can never be driven negative\n");
    {
        std::string body = "version 1\nunits per_second\nvars S P X\nrate growth = -1000\n";
        writeFile("s3.sym", body.c_str());
        complab_sym::Program D;
        ckTrue("a strongly negative growth law loads", complab_sym::load(D, "s3.sym", &err));
        complab_sym::Runtime &R = complab_sym::runtime();
        R.byMicrobe.assign(1, complab_sym::Binding());
        R.byMicrobe[0].prog = &D;
        R.byMicrobe[0].subsOfVar.assign(3, -1);
        R.byMicrobe[0].subsOfVar[0] = 0; R.byMicrobe[0].subsOfVar[1] = 1;
        R.byMicrobe[0].subsOfRate.assign(1, -1);
        R.evaluations = R.clamped = 0;
        Rig r; r.set(10.0, 1.0, 2.0); r.run(1.0);
        ck("decay is capped at the biomass present", r.dB(), -2.0);
    }

    std::remove("s.sym"); std::remove("s_hr.sym"); std::remove("s2.sym"); std::remove("s3.sym");
    std::printf("\n%s\n", fails ? "SOME CHECKS FAILED" : "all processor checks passed");
    return fails;
}
