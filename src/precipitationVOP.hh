/* ============================================================================
 * precipitationVOP.hh  —  Single-mineral VOP (Volume-Of-Pixel) precipitation
 *                         with pore-clogging flow feedback.
 * ----------------------------------------------------------------------------
 * Drop-in header for the CompLaB3D (Palabos) solver. Adds Kang-style pore-scale
 * mineral precipitation on top of the existing immobile-solid machinery.
 *
 * Method (Yoon, Kang & Valocchi 2015, RiMG 80; Kang et al. 2003/2006):
 *   1. SURFACE reaction: precipitation A -> P is heterogeneous — it only occurs
 *      at fluid–solid INTERFACE voxels (a pore voxel touching a grain or an
 *      already-precipitated voxel).  => precipitate nucleates on the walls and
 *      grows inward, instead of filling the pore homogeneously.
 *   2. VOP node conversion: each interface voxel accumulates solid P (the
 *      existing immobile-substrate density).  When P >= max_precipRho the voxel
 *      is "full" and its mask flips pore -> solid.  The new solid voxel then
 *      counts as a wall for its neighbours, so the interface (and growth) steps
 *      one voxel inward — the moving mineral front.
 *   3. Flow feedback: a converted (now-solid) voxel is handed to the EXISTING
 *      updateNsLatticesDynamics3D, which swaps its flow dynamics to BounceBack
 *      (impermeable) or reduced-omega IncBGK (partially permeable), then the NS
 *      lattice re-converges and the flow reroutes.  Clogging is detected by the
 *      existing outlet-velocity / percolation check.
 *
 * This header only provides the two NEW processors (surface reaction + node
 * conversion).  The flow-feedback step reuses updateNsLatticesDynamics3D that is
 * already in complab3d_processors_part2.hh.  See the companion integration notes
 * for where to call these in the soiltestcomplab.cpp main loop.
 *
 * Lattice-layout conventions match run_abiotic_kinetics exactly:
 *   substrate reaction lattices : [C0..C(n-1), dC0..dC(n-1), mask]  (2n+1 total)
 *   node-conversion lattices     : [P_solid_substrate, mask]        (2 total)
 * Mask encoding matches updateLocalMaskNtotalLattices3D:
 *   density(mask) stored as g[0]=(m-1)/4, g[1..6]=(m-1)/8  (=> sum = m)
 * ============================================================================ */
#ifndef PRECIPITATION_VOP_HH
#define PRECIPITATION_VOP_HH

#include <vector>
#include <cmath>

// The abiotic rate law A -> P lives in the test's defineAbioticKinetics.hh,
// exactly as for run_abiotic_kinetics:  void defineAbioticRxnKinetics(conc, subs_rate, mask);

/* ----------------------------------------------------------------------------
 * (1) surfaceAbioticKinetics3D
 *     Same dC plumbing as run_abiotic_kinetics, but the reaction is applied
 *     ONLY at fluid–solid interface voxels when surfaceOnly != 0.
 *     Interface = a non-solid/non-bb voxel with >=1 face-neighbour that is
 *     `solid` or `bb` (grain wall OR already-precipitated voxel, since a
 *     converted precipitate carries the `solid` mask code).
 *     Set surfaceOnly = 0 to recover the original bulk (homogeneous) behaviour.
 * -------------------------------------------------------------------------- */
template<typename T, template<typename U> class Descriptor>
class surfaceAbioticKinetics3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    surfaceAbioticKinetics3D(plint nx_, plint ny_, plint nz_, plint subsNum_, T dt_,
                             plint solid_, plint bb_, int surfaceOnly_)
        : nx(nx_), ny(ny_), nz(nz_), subsNum(subsNum_), dt(dt_),
          solid(solid_), bb(bb_), surfaceOnly(surfaceOnly_),
          dCloc(subsNum_), maskLloc(2*subsNum_)
    {}

    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        Dot3D absoluteOffset = lattices[0]->getLocation();
        Dot3D maskOffset = computeRelativeDisplacement(*lattices[0], *lattices[maskLloc]);

        // 6 face-neighbour displacements (D3Q7 style)
        const plint dloc[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};

        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;
            if (absX > 0 && absX < nx-1) {                    // interior in flow direction (as in run_abiotic_kinetics)
                for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                    plint absY = iY + absoluteOffset.y;
                    for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                        plint absZ = iZ + absoluteOffset.z;
                        plint iXm=iX+maskOffset.x, iYm=iY+maskOffset.y, iZm=iZ+maskOffset.z;
                        plint mask = util::roundToInt(lattices[maskLloc]->get(iXm,iYm,iZm).computeDensity());

                        if (mask != solid && mask != bb) {
                            // ---- interface gate ----
                            bool isInterface = (surfaceOnly == 0);   // surfaceOnly==0 => react everywhere
                            if (surfaceOnly != 0) {
                                for (int d=0; d<6 && !isInterface; ++d) {
                                    plint ax=absX+dloc[d][0], ay=absY+dloc[d][1], az=absZ+dloc[d][2];
                                    if (ax<0||ax>nx-1||ay<0||ay>ny-1||az<0||az>nz-1) continue;
                                    plint nb = util::roundToInt(lattices[maskLloc]->get(
                                        iXm+dloc[d][0], iYm+dloc[d][1], iZm+dloc[d][2]).computeDensity());
                                    if (nb==solid || nb==bb) isInterface = true;
                                }
                            }
                            if (!isInterface) continue;

                            // ---- reaction (identical dC handling to run_abiotic_kinetics) ----
                            std::vector<Dot3D> vec_offset;
                            for (plint iT=0; iT<maskLloc; ++iT)
                                vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT]));

                            std::vector<T> conc;
                            for (plint iS=0; iS<subsNum; ++iS) {
                                T c0 = lattices[iS]->get(iX+vec_offset[iS].x, iY+vec_offset[iS].y, iZ+vec_offset[iS].z).computeDensity();
                                if (c0 < thrd) c0 = 0;
                                conc.push_back(c0);
                            }

                            std::vector<T> subs_rate(subsNum, 0.0);
                            defineAbioticRxnKinetics(conc, subs_rate, mask);

                            for (plint iS=0; iS<subsNum; ++iS) {
                                T dC = subs_rate[iS] * dt;
                                if (dC > thrd || dC < -thrd) {
                                    plint iXd=iX+vec_offset[iS+dCloc].x, iYd=iY+vec_offset[iS+dCloc].y, iZd=iZ+vec_offset[iS+dCloc].z;
                                    Array<T,7> g;
                                    lattices[iS+dCloc]->get(iXd,iYd,iZd).getPopulations(g);
                                    g[0]+=(T)(dC)/4; g[1]+=(T)(dC)/8; g[2]+=(T)(dC)/8;
                                    g[3]+=(T)(dC)/8; g[4]+=(T)(dC)/8; g[5]+=(T)(dC)/8; g[6]+=(T)(dC)/8;
                                    lattices[iS+dCloc]->get(iXd,iYd,iZd).setPopulations(g);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // MUST be bulk (not bulkAndEnvelope): this processor reads face-neighbour mask
    // cells, so it may only run on bulk cells — a neighbour of an envelope cell is
    // outside the allocated block and dereferencing it segfaults.
    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }
    virtual surfaceAbioticKinetics3D<T,Descriptor>* clone() const {
        return new surfaceAbioticKinetics3D<T,Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint i=0; i<(plint)modified.size(); ++i) modified[i] = modif::nothing;
        for (plint i=dCloc; i<2*subsNum; ++i) modified[i] = modif::staticVariables;   // dC lattices written
    }
private:
    plint nx, ny, nz, subsNum;
    T dt;
    plint solid, bb;
    int surfaceOnly;
    plint dCloc, maskLloc;
};

/* ----------------------------------------------------------------------------
 * (2) precipNodeConversion3D
 *     VOP solid-node update.  When the accumulated immobile solid P at a pore
 *     voxel reaches max_precipRho, flip that voxel's mask pore -> precipCode.
 *     Pass precipCode = the solver's `solid` (no_dynamics) mask code so that:
 *       - the reaction/skip logic (mask==solid) stops reacting there,
 *       - it counts as a wall for the interface neighbour check (front advances),
 *       - updateNsLatticesDynamics3D turns its flow node into BounceBack / low-ω.
 *     Lattices: [ P_solid_substrate_lattice , maskLattice ].
 *     The main loop detects the change via its existing old/new total-mask sum.
 * -------------------------------------------------------------------------- */
template<typename T, template<typename U> class Descriptor>
class precipNodeConversion3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    /* precipPhaseId: when dissolution is enabled, the converted voxel is also
     * stamped with this phase id on lattices[2] (the phase lattice), so
     * dissolutionVOP.hh knows the voxel is precipitate rather than original rock
     * and may reopen it later.  Pass 0 and omit lattices[2] to keep the original
     * two-lattice behaviour. */
    precipNodeConversion3D(plint nx_, plint ny_, plint nz_, T max_precipRho_,
                           plint solid_, plint bb_, plint precipCode_,
                           plint precipPhaseId_ = 0)
        : nx(nx_), ny(ny_), nz(nz_), max_precipRho(max_precipRho_),
          solid(solid_), bb(bb_), precipCode(precipCode_), precipPhaseId(precipPhaseId_)
    {}

    // lattices[0] = immobile solid (P) substrate lattice ; lattices[1] = mask lattice
    //             ; lattices[2] = phase lattice (optional, only when precipPhaseId > 0)
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        Dot3D maskOffset = computeRelativeDisplacement(*lattices[0], *lattices[1]);
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint iXm=iX+maskOffset.x, iYm=iY+maskOffset.y, iZm=iZ+maskOffset.z;
                    plint mask = util::roundToInt(lattices[1]->get(iXm,iYm,iZm).computeDensity());
                    if (mask != solid && mask != bb) {                 // only pore/biofilm voxels
                        T P = lattices[0]->get(iX,iY,iZ).computeDensity();
                        if (P >= max_precipRho) {                       // voxel is "full" -> convert to solid
                            Array<T,7> g;
                            g[0]=(T)(precipCode-1)/4;
                            g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T)(precipCode-1)/8;
                            lattices[1]->get(iXm,iYm,iZm).setPopulations(g);

                            // Stamp the phase so dissolution can find it later.
                            if (precipPhaseId > 0 && lattices.size() > 2) {
                                Dot3D pOff = computeRelativeDisplacement(*lattices[0], *lattices[2]);
                                Array<T,7> h;
                                h[0]=(T)(precipPhaseId-1)/4;
                                h[1]=h[2]=h[3]=h[4]=h[5]=h[6]=(T)(precipPhaseId-1)/8;
                                lattices[2]->get(iX+pOff.x,iY+pOff.y,iZ+pOff.z).setPopulations(h);
                            }
                        }
                    }
                }
            }
        }
    }

    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }
    virtual precipNodeConversion3D<T,Descriptor>* clone() const {
        return new precipNodeConversion3D<T,Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;          // P substrate read-only
        modified[1] = modif::staticVariables;  // mask written
        if (modified.size() > 2) modified[2] = modif::staticVariables;  // phase written
    }
private:
    plint nx, ny, nz;
    T max_precipRho;
    plint solid, bb, precipCode, precipPhaseId;
};

#endif // PRECIPITATION_VOP_HH
