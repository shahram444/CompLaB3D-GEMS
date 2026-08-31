/* ================================================================================================
 *  complab3d_processors_graphnet.hh  --  the per-voxel driver for a bipartite graph network
 *
 *  Part of the CompLaB program.  GNU Affero General Public License v3 or later.
 *  Meile Lab, University of Georgia.  shahram.asgari@uga.edu
 *
 *  ------------------------------------------------------------------------------------------------
 *  This is run_symbolic3D with the evaluator swapped.  The lattice layout, the biomass gathering,
 *  the shared substrate budget and the D3Q7 deposit are identical, and deliberately so: a graph
 *  network organism, a symbolic organism and a surrogate organism sitting in the same voxel must
 *  draw on the same supply under the same rules, or whichever ran second could push a concentration
 *  negative.  If one of these three budget loops is ever changed, all three change.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT COMES BACK FROM THE NETWORK
 *
 *  One rate per species it was trained on, and optionally a growth rate.  That is the difference
 *  from the surrogate, which returns growth alone and leaves the substrate fields to the linear
 *  program.  Here there is no linear program: the network is the whole reaction step.
 *
 *  UNITS.  Network::eval() has already applied the file's units scale, so a substrate rate is
 *  mol/L/s and a growth rate is 1/s -- the defineKinetics.hh convention, which is right because a
 *  .gnn file replaces a hand-written rate law rather than a metabolic model.  There is no
 *  biomass_molar_mass anywhere on this path and no fluxToDeltaC(): a fitted rate law is not posed
 *  per gram dry weight, and putting it through the FBA unit chain would be wrong by a factor M_B.
 *
 *  ------------------------------------------------------------------------------------------------
 *  ONE LINE THIS FILE DELIBERATELY DOES NOT COPY.  runFBA_cobrapy3D and run_surrogate3D both end
 *  their budget loop with
 *
 *      if (draw > T()) draw = T();          // consumption only
 *
 *  which discards net release, so an organism on those paths can consume but cannot excrete.
 *  runFBA_glpk3D has no such line.  A network fitted to a reaction system predicts products as well
 *  as reactants -- that is the point of giving it the stoichiometry -- so the clamp here is
 *  one-sided the other way: draw is limited by what the voxel holds, and production passes through.
 * ================================================================================================ */

#ifndef COMPLAB3D_PROCESSORS_GRAPHNET_HH
#define COMPLAB3D_PROCESSORS_GRAPHNET_HH

#include "complab3d_metabolic.hh"   /* MetabolicConfig: useTotals and the equilibrium tableau */
#include "complab3d_graphnet.hh"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class run_graphnet3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    /* Same argument list as run_surrogate3D and run_symbolic3D, so the call site in complab.cpp is
     * a copy with one word changed.  cfg carries useTotals and the equilibrium tableau, so a
     * network reads the same total concentration an FBA organism would. */
    run_graphnet3D (plint nx_, plint subsNum_, plint bioNum_, plint totalMicrobes_, T dt_,
                    const std::vector<plint> &globalId_,
                    const MetabolicConfig    *cfg_,
                    plint solid_, plint bb_)
        : nx(nx_), subsNum(subsNum_), bioNum(bioNum_), totalMicrobes(totalMicrobes_), dt(dt_),
          globalId(globalId_), cfg(cfg_), solid(solid_), bb(bb_),
          dCloc(subsNum_ + bioNum_),
          dBloc(2*subsNum_ + bioNum_),
          maskLloc(2*(subsNum_ + bioNum_))
    {}

    virtual void process (Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices)
    {
        const T thrd = (T) 1e-12;

        Dot3D absoluteOffset = lattices[0]->getLocation();
        Dot3D maskOffset     = computeRelativeDisplacement(*lattices[0], *lattices[maskLloc]);

        std::vector<Dot3D> off;
        off.reserve(maskLloc);
        for (plint iT = 0; iT < maskLloc; ++iT)
            off.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT]));

        std::vector<T> conc (subsNum, T());
        std::vector<T> avail(subsNum, T());
        std::vector<T> dC   (subsNum, T());
        std::vector<T> req  (subsNum, T());   /* what the networks asked for, before the budget */
        std::vector<T> bmass(bioNum,  T());
        std::vector<T> dB   (bioNum,  T());
        std::vector<plint> bLoc;

        /* scratch for the forward pass, allocated once for the whole sweep rather than per voxel.
         * eval() is called in the innermost loop, so the buffers it works in are hoisted out here;
         * letting it allocate its own would cost more than the arithmetic. */
        std::vector<double> in, out;
        complab_gnn::Network::Scratch scratch;

        complab_gnn::Runtime &RT = complab_gnn::runtime();

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            const plint absX = iX + absoluteOffset.x;
            if (absX <= 0 || absX >= nx-1) continue;          /* the two ghost columns */

            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {

                    const plint mask = util::roundToInt(
                        lattices[maskLloc]->get(iX+maskOffset.x, iY+maskOffset.y, iZ+maskOffset.z).computeDensity());
                    if (mask == solid || mask == bb) continue;

                    /* ---- which organisms are actually here ---------------------------------- */
                    bLoc.clear();
                    for (plint iB = 0; iB < bioNum; ++iB) {
                        const T b = lattices[subsNum+iB]->get(iX+off[subsNum+iB].x,
                                                              iY+off[subsNum+iB].y,
                                                              iZ+off[subsNum+iB].z).computeDensity();
                        bmass[iB] = (b > thrd) ? b : T();
                        dB[iB]    = T();
                        if (bmass[iB] > thrd) bLoc.push_back(iB);
                    }
                    if (bLoc.empty()) continue;

                    /* ---- local chemistry ---------------------------------------------------- */
                    for (plint iS = 0; iS < subsNum; ++iS) {
                        const T c = lattices[iS]->get(iX+off[iS].x, iY+off[iS].y, iZ+off[iS].z).computeDensity();
                        conc[iS]  = (c > thrd) ? c : T();
                        dC[iS]    = T();
                        req[iS]   = T();
                        avail[iS] = T();
                    }
                    for (plint iS = 0; iS < subsNum; ++iS)
                        avail[iS] = cfg->useTotals ? cfg->totals.total(conc, iS) : conc[iS];

                    /* ---- one forward pass per organism present ------------------------------ */
                    for (size_t k = 0; k < bLoc.size(); ++k) {
                        const plint iB = bLoc[k];
                        const plint gM = globalId[(size_t) iB];
                        if (gM < 0 || gM >= (plint) RT.byMicrobe.size()) continue;
                        const complab_gnn::Binding &B = RT.byMicrobe[(size_t) gM];
                        if (B.net == 0 || !B.net->valid()) continue;

                        /* The network's species order is its own; subsOfSpecies maps it onto the
                         * lattice.  A species the file names but the simulation does not carry was
                         * refused at load, so every entry here is a real index. */
                        in.assign((size_t) B.net->nS, 0.0);
                        for (size_t s = 0; s < in.size(); ++s) {
                            const int iS = B.subsOfSpecies[s];
                            if (iS >= 0) in[s] = (double) avail[iS];
                        }

                        B.net->eval(in, out, scratch, true, &RT.clamped);
                        ++RT.evaluations;

                        /* mol/L/s per species, and 1/s for growth, both already unit-scaled */
                        for (size_t s = 0; s < in.size() && s < out.size(); ++s) {
                            const int iS = B.subsOfSpecies[s];
                            if (iS >= 0) req[iS] += (T) out[s] * dt;          /* signed increment */
                        }
                        if (B.growthSlot >= 0 && B.growthSlot < (int) out.size())
                            dB[iB] += (T) out[(size_t) B.growthSlot] * bmass[iB] * dt;
                    }

                    /* ---- one shared substrate budget ---------------------------------------
                     * Structurally identical to run_surrogate3D's and run_symbolic3D's.  dC is
                     * written ONLY here; req holds what the networks asked for.  Accumulating
                     * straight into dC would double-count under useTotals, because drawDown()
                     * writes dC itself. */
                    for (plint iS = 0; iS < subsNum; ++iS) {
                        if (cfg->fix_concentration[iS]) continue;
                        T draw = req[iS];                       /* negative = consumed */
                        if (std::fabs(draw) <= thrd) continue;
                        const T headroom = avail[iS] + (cfg->useTotals ? T() : dC[iS]);
                        if (draw + headroom < T()) draw = -headroom;
                        /* NO "if (draw > T()) draw = T();" here.  See the header comment. */
                        if (draw < T() && cfg->useTotals) cfg->totals.drawDown(conc, iS, -draw, dC);
                        else                              dC[iS] += draw;
                    }

                    /* ---- deposit ------------------------------------------------------------
                     * D3Q7 weights: rest 1/4, six directions 1/8 each, summing to one. */
                    for (plint iS = 0; iS < subsNum; ++iS) {
                        if (std::fabs(dC[iS]) <= thrd) continue;
                        Array<T,7> g;
                        Cell<T,Descriptor> &cell = lattices[iS+dCloc]->get(
                            iX+off[iS+dCloc].x, iY+off[iS+dCloc].y, iZ+off[iS+dCloc].z);
                        cell.getPopulations(g);
                        g[0]+=dC[iS]/4; g[1]+=dC[iS]/8; g[2]+=dC[iS]/8; g[3]+=dC[iS]/8;
                        g[4]+=dC[iS]/8; g[5]+=dC[iS]/8; g[6]+=dC[iS]/8;
                        cell.setPopulations(g);
                    }
                    for (plint iB = 0; iB < bioNum; ++iB) {
                        if (std::fabs(dB[iB]) <= thrd) continue;
                        if (dB[iB] + bmass[iB] < T()) dB[iB] = -bmass[iB];   /* never negative */
                        Array<T,7> g;
                        Cell<T,Descriptor> &cell = lattices[iB+dBloc]->get(
                            iX+off[iB+dBloc].x, iY+off[iB+dBloc].y, iZ+off[iB+dBloc].z);
                        cell.getPopulations(g);
                        g[0]+=dB[iB]/4; g[1]+=dB[iB]/8; g[2]+=dB[iB]/8; g[3]+=dB[iB]/8;
                        g[4]+=dB[iB]/8; g[5]+=dB[iB]/8; g[6]+=dB[iB]/8;
                        cell.setPopulations(g);
                    }
                }
            }
        }
    }

    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulkAndEnvelope; }

    virtual run_graphnet3D<T,Descriptor>* clone() const {
        return new run_graphnet3D<T,Descriptor>(*this);
    }

    void getTypeOfModification (std::vector<modif::ModifT> &modified) const {
        for (size_t i = 0; i < modified.size(); ++i) modified[i] = modif::nothing;
        for (plint iT = dCloc; iT < maskLloc; ++iT) modified[iT] = modif::staticVariables;
    }

private:
    plint nx, subsNum, bioNum, totalMicrobes;
    T dt;
    std::vector<plint> globalId;
    const MetabolicConfig *cfg;
    plint solid, bb;
    plint dCloc, dBloc, maskLloc;
};

}  // namespace plb

#endif  // COMPLAB3D_PROCESSORS_GRAPHNET_HH
