/* ================================================================================================
 *  test_dissolution.cpp  --  the dissolution processor, RUN, with the mineral-decrement fix
 *
 *  The defect this pins: the mineral decrement used to be written into dC at the dissolving
 *  voxel.  That voxel is solid by construction, and the step that applies dC skips solid voxels,
 *  so the decrement was discarded while the products -- written to the pore neighbours -- were
 *  applied.  Mineral never fell; solute appeared from nothing, every step, forever.
 *
 *  So the load-bearing check here is simply: after one step, is there less mineral than before,
 *  and does the amount that left the solid equal the amount that arrived in the pores?
 * ================================================================================================ */
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include "palabos_stub.hh"
namespace plb { PlbOut pcout; }
using namespace plb;

/* the user hook the processor calls; declared before the header that calls it */
void defineDissolutionRate(plint, const std::vector<double>&, double,
                           std::vector<double>&, double&, plint);

#include "dissolutionVOP.hh"

static int fails = 0;
static void ck(const char *what, double got, double want, double tol = 1e-12) {
    const bool ok = std::fabs(got - want) <= tol * (1.0 + std::fabs(want));
    if (!ok) ++fails;
    std::printf("  %-56s got %-14.7g want %-14.7g %s\n", what, got, want, ok ? "ok" : "** FAIL **");
}
static void ckTrue(const char *what, bool ok) {
    if (!ok) ++fails;
    std::printf("  %-56s %s\n", what, ok ? "ok" : "** FAIL **");
}

/* The rate law hook the processor calls.  Kept trivial so the arithmetic is checkable:
 *      FeS(s) + H+  ->  Fe2+ + HS-      rate = k * [H+] * M
 * substrate 0 = Fe, 1 = HS, 2 = H, 3 = FeS(mineral) */
static double g_k = 0.0;
void defineDissolutionRate(plint, const std::vector<double>& C, double mineral,
                           std::vector<double>& subsR, double& mineralR, plint)
{
    const double r = g_k * C[2] * mineral;          // mol/L/s, positive
    mineralR = -r;
    subsR[0] = +r; subsR[1] = +r; subsR[2] = -r; subsR[3] = 0.0;
}

int main() {
    const int S = 4, sM = 3;
    /* layout: [C0..C3][dC0..dC3][mask][phase]  = 2S + 2 */
    std::vector< BlockLattice3D<double,D3Q7>* > L;
    for (int i = 0; i < 2*S + 2; ++i) L.push_back(new BlockLattice3D<double,D3Q7>());

    DissolutionConfig cfg;
    cfg.enabled = true; cfg.surface_only = 1; cfg.reopen_fraction = 0.9;
    SolidPhase ph; ph.id = 1; ph.material_number = 0; ph.substrate = sM;
    ph.full_density = 48.9; ph.initial_fill = 48.9; ph.is_precipitate = true; ph.name = "FeS";
    cfg.phases.push_back(ph);

    /* one solid, phase-carrying voxel holding mineral, with H+ in the pore around it */
    L[2]->c.rho = 2.0e-2;          // H+  (read from the OPEN NEIGHBOURS, see below)
    L[sM]->c.rho = 48.9;           // the mineral inventory at the voxel
    L[2*S]->c.rho = 0.0;           // mask: fluid, so the six faces read as open
    L[2*S+1]->c.rho = 1.0;         // phase id 1
    g_k = 5.0e-2;

    std::vector<plint> pore(1, 2);
    const double dt = 1.2e-2;

    std::printf("--- THE FIX: the mineral inventory actually falls\n");
    {
        surfaceDissolutionKinetics3D<double,D3Q7> p(3, 3, 3, S, dt, 1, 2, pore, &cfg);
        Box3D d; d.x0 = d.x1 = 1; d.y0 = d.y1 = 0; d.z0 = d.z1 = 0;
        const double before = L[sM]->c.computeDensity();
        p.process(d, L);
        const double after = L[sM]->c.computeDensity();

        /* rate = k * H * M = 5e-2 * 2e-2 * 48.9 = 4.89e-2 mol/L/s
         * dM   = -4.89e-2 * 1.2e-2 = -5.868e-4, well under the 50% limiter */
        const double expected = -5.0e-2 * 2.0e-2 * 48.9 * dt;
        ckTrue("the mineral lattice was written at all", after != before);
        ck("and it fell by exactly rate * dt", after - before, expected);
        ck("nothing was left stranded in dC[mineral]", L[sM+S]->deposited(), 0.0);
    }

    std::printf("--- what leaves the solid equals what arrives in the pores\n");
    {
        /* six open faces, so each product is split six ways; the stub gives every neighbour
         * the same cell, so the six deposits land on one lattice and simply add back up */
        const double dM = -5.0e-2 * 2.0e-2 * 48.9 * dt;
        ck("Fe released into the pores  = -dM", L[0+S]->deposited(), -dM);
        ck("HS released into the pores  = -dM", L[1+S]->deposited(), -dM);
        ck("H+ consumed from the pores  = +dM", L[2+S]->deposited(),  dM);
    }

    std::printf("--- the 50%% per-step limiter still binds, and scales products with it\n");
    {
        for (int i = 0; i < 2*S + 2; ++i) { for (int j = 0; j < 7; ++j) L[i]->c.p[j] = 0.0; }
        L[2]->c.rho = 2.0e-2; L[2*S]->c.rho = 0.0; L[2*S+1]->c.rho = 1.0;
        L[sM]->c.rho = 48.9;
        g_k = 5.0e4;                                  /* absurdly fast on purpose */
        surfaceDissolutionKinetics3D<double,D3Q7> p(3, 3, 3, S, dt, 1, 2, pore, &cfg);
        Box3D d; d.x0 = d.x1 = 1; d.y0 = d.y1 = 0; d.z0 = d.z1 = 0;
        const double before = L[sM]->c.computeDensity();
        p.process(d, L);
        const double removed = before - L[sM]->c.computeDensity();
        ck("at most half the inventory goes in one step", removed, 0.5 * 48.9);
        ck("and the products are scaled by the same factor", L[0+S]->deposited(), 0.5 * 48.9);
    }

    std::printf("--- an enclosed voxel does not dissolve when surface_only is on\n");
    {
        for (int i = 0; i < 2*S + 2; ++i) { for (int j = 0; j < 7; ++j) L[i]->c.p[j] = 0.0; }
        L[2]->c.rho = 2.0e-2; L[sM]->c.rho = 48.9; L[2*S+1]->c.rho = 1.0;
        L[2*S]->c.rho = 1.0;                          /* every neighbour reads solid */
        g_k = 5.0e-2;
        surfaceDissolutionKinetics3D<double,D3Q7> p(3, 3, 3, S, dt, 1, 2, pore, &cfg);
        Box3D d; d.x0 = d.x1 = 1; d.y0 = d.y1 = 0; d.z0 = d.z1 = 0;
        const double before = L[sM]->c.computeDensity();
        p.process(d, L);
        ck("mineral untouched", L[sM]->c.computeDensity() - before, 0.0);
    }

    for (size_t i = 0; i < L.size(); ++i) delete L[i];
    std::printf("\n%s\n", fails ? "SOME CHECKS FAILED" : "all dissolution checks passed");
    return fails;
}
