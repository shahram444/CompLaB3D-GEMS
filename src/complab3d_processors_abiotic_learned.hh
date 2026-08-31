/* ================================================================================================
 *  complab3d_processors_abiotic_learned.hh
 *      the abiotic sweep for a .sym rate law and for a .gnn network
 *
 *  Part of the CompLaB program.  GNU Affero General Public License v3 or later.
 *  Meile Lab, University of Georgia.  shahram.asgari@uga.edu
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHY THIS FILE EXISTS AT ALL
 *
 *  run_symbolic3D and run_graphnet3D were built from the FBA processor, and the FBA processor is a
 *  MICROBE processor: the first thing it does in a voxel is gather biomass, and if there is none it
 *  moves on.  That is right for an organism -- a cell that is not there cannot eat -- and wrong for
 *  abiotic chemistry, which happens in open water whether anything is alive nearby or not.
 *
 *  Put an abiotic reaction in a .sym file and run it through the biotic path and it would fire
 *  inside the biofilm and stop dead at its edge.  Nothing would warn you.  You would simply get
 *  reaction where the organisms are and none anywhere else, which looks like a result.
 *
 *  So these two are the same evaluators over a different sweep: every fluid voxel, no biomass
 *  gathering, no growth output, no per-microbe binding.  They are the learned-rate-law counterpart
 *  of run_abiotic_kinetics, and they sit on the same lattice vector it uses:
 *
 *      substrates ... , dC ... , mask            (2 * subsNum + 1 lattices)
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT AN ABIOTIC FILE MAY AND MAY NOT CONTAIN
 *
 *  Every variable must be a substrate, and every rate must name a substrate.  A "growth" rate is
 *  refused at load time rather than ignored: an abiotic file with a growth line means whoever wrote
 *  it believed an organism was involved, and silently dropping that line would hide the mistake.
 *  The same goes for a variable naming a microbe.  Those checks live in complab.cpp, where the
 *  names are known.
 * ================================================================================================ */

#ifndef COMPLAB3D_PROCESSORS_ABIOTIC_LEARNED_HH
#define COMPLAB3D_PROCESSORS_ABIOTIC_LEARNED_HH

#include "complab3d_metabolic.hh"   /* MetabolicConfig: useTotals and the equilibrium tableau */
#include "complab3d_symbolic.hh"
#include "complab3d_graphnet.hh"

namespace plb {

/* ================================================================================================
 *  A .sym rate law, applied everywhere
 * ================================================================================================ */
template<typename T, template<typename U> class Descriptor>
class run_symbolic_abiotic3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    run_symbolic_abiotic3D (plint nx_, plint subsNum_, T dt_,
                            const MetabolicConfig *cfg_, plint solid_, plint bb_)
        : nx(nx_), subsNum(subsNum_), dt(dt_), cfg(cfg_), solid(solid_), bb(bb_),
          dCloc(subsNum_), maskLloc(2*subsNum_)
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

        std::vector<T> conc (subsNum, T()), avail(subsNum, T());
        std::vector<T> dC   (subsNum, T()), req  (subsNum, T());
        std::vector<double> vars, out;

        complab_sym::Runtime &RT = complab_sym::runtime();
        const complab_sym::Binding &B = RT.abiotic;
        if (B.prog == 0 || !B.prog->valid()) return;

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            const plint absX = iX + absoluteOffset.x;
            if (absX <= 0 || absX >= nx-1) continue;          /* the two ghost columns */

            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {

                    const plint mask = util::roundToInt(
                        lattices[maskLloc]->get(iX+maskOffset.x, iY+maskOffset.y, iZ+maskOffset.z).computeDensity());
                    if (mask == solid || mask == bb) continue;

                    /* NO biomass test here.  That single missing line is the entire difference
                     * between this file and complab3d_processors_symbolic.hh. */

                    for (plint iS = 0; iS < subsNum; ++iS) {
                        const T c = lattices[iS]->get(iX+off[iS].x, iY+off[iS].y, iZ+off[iS].z).computeDensity();
                        conc[iS]  = (c > thrd) ? c : T();
                        dC[iS]    = T();
                        req[iS]   = T();
                        avail[iS] = T();
                    }
                    for (plint iS = 0; iS < subsNum; ++iS)
                        avail[iS] = cfg->useTotals ? cfg->totals.total(conc, iS) : conc[iS];

                    vars.assign(B.prog->vars.size(), 0.0);
                    for (size_t v = 0; v < vars.size(); ++v) {
                        const int iS = B.subsOfVar[v];
                        vars[v] = (iS >= 0) ? (double) avail[iS] : 0.0;   /* no biomass on this path */
                    }

                    complab_sym::evaluate(*B.prog, vars, out, true, &RT.clamped);
                    ++RT.evaluations;

                    for (size_t r = 0; r < out.size() && r < B.subsOfRate.size(); ++r) {
                        const int iS = B.subsOfRate[r];
                        if (iS >= 0) req[iS] += (T) out[r] * dt;
                    }

                    /* the same budget the biotic paths use, so a substrate shared between an
                     * abiotic law and an organism obeys one rule, not two */
                    for (plint iS = 0; iS < subsNum; ++iS) {
                        if (cfg->fix_concentration[iS]) continue;
                        T draw = req[iS];
                        if (std::fabs(draw) <= thrd) continue;
                        const T headroom = avail[iS] + (cfg->useTotals ? T() : dC[iS]);
                        if (draw + headroom < T()) draw = -headroom;
                        if (draw < T() && cfg->useTotals) cfg->totals.drawDown(conc, iS, -draw, dC);
                        else                              dC[iS] += draw;
                    }

                    /* D3Q7 weights: rest 1/4, six directions 1/8 each */
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
                }
            }
        }
    }

    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulkAndEnvelope; }
    virtual run_symbolic_abiotic3D<T,Descriptor>* clone() const {
        return new run_symbolic_abiotic3D<T,Descriptor>(*this);
    }
    void getTypeOfModification (std::vector<modif::ModifT> &modified) const {
        for (size_t i = 0; i < modified.size(); ++i) modified[i] = modif::nothing;
        for (plint iT = dCloc; iT < maskLloc; ++iT) modified[iT] = modif::staticVariables;
    }

private:
    plint nx, subsNum;
    T dt;
    const MetabolicConfig *cfg;
    plint solid, bb;
    plint dCloc, maskLloc;
};


/* ================================================================================================
 *  A .gnn network, applied everywhere
 * ================================================================================================ */
template<typename T, template<typename U> class Descriptor>
class run_graphnet_abiotic3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    run_graphnet_abiotic3D (plint nx_, plint subsNum_, T dt_,
                            const MetabolicConfig *cfg_, plint solid_, plint bb_)
        : nx(nx_), subsNum(subsNum_), dt(dt_), cfg(cfg_), solid(solid_), bb(bb_),
          dCloc(subsNum_), maskLloc(2*subsNum_)
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

        std::vector<T> conc (subsNum, T()), avail(subsNum, T());
        std::vector<T> dC   (subsNum, T()), req  (subsNum, T());
        std::vector<double> in, out;
        complab_gnn::Network::Scratch scratch;

        complab_gnn::Runtime &RT = complab_gnn::runtime();
        const complab_gnn::Binding &B = RT.abiotic;
        if (B.net == 0 || !B.net->valid()) return;

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            const plint absX = iX + absoluteOffset.x;
            if (absX <= 0 || absX >= nx-1) continue;

            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {

                    const plint mask = util::roundToInt(
                        lattices[maskLloc]->get(iX+maskOffset.x, iY+maskOffset.y, iZ+maskOffset.z).computeDensity());
                    if (mask == solid || mask == bb) continue;

                    for (plint iS = 0; iS < subsNum; ++iS) {
                        const T c = lattices[iS]->get(iX+off[iS].x, iY+off[iS].y, iZ+off[iS].z).computeDensity();
                        conc[iS]  = (c > thrd) ? c : T();
                        dC[iS]    = T();
                        req[iS]   = T();
                        avail[iS] = T();
                    }
                    for (plint iS = 0; iS < subsNum; ++iS)
                        avail[iS] = cfg->useTotals ? cfg->totals.total(conc, iS) : conc[iS];

                    in.assign((size_t) B.net->nS, 0.0);
                    for (size_t s = 0; s < in.size(); ++s) {
                        const int iS = B.subsOfSpecies[s];
                        if (iS >= 0) in[s] = (double) avail[iS];
                    }

                    B.net->eval(in, out, scratch, true, &RT.clamped);
                    ++RT.evaluations;

                    /* Only the species rates are read.  A growth slot, if the file has one, is
                     * ignored here because there is nothing on this path for it to grow -- and
                     * complab.cpp refuses such a file at load time anyway. */
                    for (size_t s = 0; s < in.size() && s < out.size(); ++s) {
                        const int iS = B.subsOfSpecies[s];
                        if (iS >= 0) req[iS] += (T) out[s] * dt;
                    }

                    for (plint iS = 0; iS < subsNum; ++iS) {
                        if (cfg->fix_concentration[iS]) continue;
                        T draw = req[iS];
                        if (std::fabs(draw) <= thrd) continue;
                        const T headroom = avail[iS] + (cfg->useTotals ? T() : dC[iS]);
                        if (draw + headroom < T()) draw = -headroom;
                        if (draw < T() && cfg->useTotals) cfg->totals.drawDown(conc, iS, -draw, dC);
                        else                              dC[iS] += draw;
                    }

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
                }
            }
        }
    }

    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulkAndEnvelope; }
    virtual run_graphnet_abiotic3D<T,Descriptor>* clone() const {
        return new run_graphnet_abiotic3D<T,Descriptor>(*this);
    }
    void getTypeOfModification (std::vector<modif::ModifT> &modified) const {
        for (size_t i = 0; i < modified.size(); ++i) modified[i] = modif::nothing;
        for (plint iT = dCloc; iT < maskLloc; ++iT) modified[iT] = modif::staticVariables;
    }

private:
    plint nx, subsNum;
    T dt;
    const MetabolicConfig *cfg;
    plint solid, bb;
    plint dCloc, maskLloc;
};

}  // namespace plb

#endif  // COMPLAB3D_PROCESSORS_ABIOTIC_LEARNED_HH
