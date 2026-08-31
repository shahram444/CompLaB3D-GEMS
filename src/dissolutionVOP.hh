/* This file is part of the CompLB3D library.
 *
 * CompLB3D is developed since 2022 by the University of Georgia (United States)
 * and Chungnam National University (South Korea).
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more
 * details.
 * ============================================================================
 *
 *  dissolutionVOP.hh
 *  --------------------------------------------------------------------------
 *  MINERAL DISSOLUTION AND PORE RE-OPENING
 *
 *  The reverse of precipitationVOP.hh.  Where that one lets a growing mineral
 *  fill a pore voxel and seal it, this one lets a mineral be eaten away until
 *  its voxel opens back up to flow.
 *
 *  It is optional and independent: <precipitation> and <dissolution> each have
 *  their own <enabled>, so you can run
 *      precipitation only   -- exactly the behaviour that shipped before
 *      dissolution only     -- start from a mineral-filled rock and leach it
 *      both                 -- a full precipitate / dissolve cycle
 *      neither              -- the solver is untouched
 *
 *  --------------------------------------------------------------------------
 *  THE IDEA: A "SOLID PHASE"
 *
 *  There are two different kinds of solid voxel in a CompLB3D domain, and until
 *  now the code could not tell them apart:
 *
 *    1. Original grains, read from geometry.dat.  The rock you started with.
 *    2. Voxels that BECAME solid during the run because a mineral precipitated
 *       and filled them.
 *
 *  That distinction matters the moment dissolution exists.  If the rule were
 *  simply "a solid voxel whose mineral has run out turns back into pore", then
 *  an inert quartz grain -- which contains no mineral and never did -- is
 *  already below any threshold you pick, and the entire grain pack would
 *  dissolve on the first step.
 *
 *  So dissolution works on declared SOLID PHASES.  A phase is four numbers:
 *
 *      material_number   which mask code this phase occupies
 *      substrate         which <immobile> substrate holds its inventory, mol/L
 *      full_density      the inventory of a voxel that is completely full,
 *                        mol/L; derive it from the mineral's molar volume,
 *                        1000 / (molar volume in cm3/mol)
 *      initial_fill      how full a voxel of this phase starts, mol/L
 *
 *  Precipitated mineral and original rock then become the SAME mechanism,
 *  differing only in how the inventory is seeded:
 *
 *      precipitated FeS   material 7, substrate 9,  full 48.9, initial 0
 *                         (starts empty, grows by the reaction you wrote)
 *      calcite grain      material 0, substrate 13, full 27.1, initial 27.1
 *                         (starts full, only shrinks)
 *      inert quartz       not declared as a phase at all
 *
 *  Anything not declared can never dissolve.  Quartz is safe by construction
 *  rather than by a special case, which is the property worth having.
 *
 *  --------------------------------------------------------------------------
 *  HOW A VOXEL IS TRACKED
 *
 *  A new lattice, `phaseLattice`, sits alongside the existing maskLattice,
 *  ageLattice and distLattice -- the same pattern the code already uses for
 *  carrying an integer field on a D3Q7 lattice.  Its value is the phase id
 *  occupying that voxel, or 0 for none.
 *
 *      seeded at start-up      from geometry, for the grain phases
 *      stamped                 by precipNodeConversion3D when a voxel seals
 *      cleared                 by dissolNodeConversion3D when a voxel reopens
 *
 *  The MASK value of a sealed voxel stays `solid`, exactly as before.  That is
 *  deliberate: every existing processor tests `mask == solid` to decide what to
 *  skip, and introducing a new mask code would have meant auditing all of them.
 *  The phase lattice is purely additive, and invisible to code that does not
 *  ask for it.
 *
 *  --------------------------------------------------------------------------
 *  WHERE THE CHEMISTRY LIVES
 *
 *  Not here.  Exactly as precipitation takes its rate from whatever you wrote
 *  in defineAbioticKinetics.hh, dissolution calls one hook in the same file:
 *
 *      void defineDissolutionRate(plint phaseId, const std::vector<double>& C,
 *                                 double mineral, std::vector<double>& subsR,
 *                                 double& mineralR, plb::plint mask);
 *
 *  You are handed the water in contact with the mineral and how much mineral is
 *  left; you return how fast it dissolves and what that releases.  This file
 *  only handles the geometry: which voxels are wetted, where the products go,
 *  and when a voxel reopens.
 *
 *  --------------------------------------------------------------------------
 *  WHERE THE PRODUCTS GO
 *
 *  A dissolving voxel is solid, so its own solute lattices are not transported
 *  -- anything released there would sit in a dead cell.  The products are
 *  therefore deposited into the OPEN face neighbours, split evenly among them.
 *  A voxel with no open neighbour is fully enclosed, cannot be reached by water,
 *  and does not dissolve at all, which is the physically right answer.
 *
 *  This is the one piece with no precipitation counterpart, and it is where
 *  mass conservation lives:  mineral removed  ==  product released, summed over
 *  the receiving neighbours.  The driver in the test notes verifies exactly that
 *  over a seal / dissolve / reopen cycle.
 *
 *  --------------------------------------------------------------------------
 *  HYSTERESIS -- why there are two thresholds and not one
 *
 *  A voxel seals at   inventory >= full_density
 *  and reopens at     inventory <  reopen_fraction * full_density
 *
 *  with reopen_fraction < 1 (0.9 by default).  With a single threshold a voxel
 *  sitting right at it would flip every update interval, and because any mask
 *  change triggers a full Navier-Stokes re-solve, that is both physically
 *  meaningless chatter and ruinously expensive.  The gap between the two
 *  thresholds prevents it.
 * ============================================================================
 */

#ifndef DISSOLUTION_VOP_HH
#define DISSOLUTION_VOP_HH

#include <vector>
#include <string>
#include <cmath>


/* ============================================================================
 *  SolidPhase -- one declared dissolvable solid
 * ============================================================================
 */
struct SolidPhase {
    plint       id;              // 1-based; 0 is reserved for "no phase"
    plint       material_number; // which mask code this phase occupies
    plint       substrate;       // which <immobile> substrate holds its inventory
    T           full_density;    // mol/L in a completely full voxel
    T           initial_fill;    // mol/L a voxel of this phase starts with
    bool        is_precipitate;  // true for the phase precipitation converts into
    std::string name;

    SolidPhase()
      : id(0), material_number(-1), substrate(-1),
        full_density((T) 0), initial_fill((T) 0), is_precipitate(false) {}
};


/* ============================================================================
 *  DissolutionConfig
 * ============================================================================
 */
struct DissolutionConfig {
    bool  enabled;
    T     reopen_fraction;    // reopen below this fraction of full_density
    int   surface_only;       // 1 = only wetted voxels react (the physical case)
    plint update_interval;    // steps between geometry re-checks
    std::vector<SolidPhase> phases;

    DissolutionConfig()
      : enabled(false), reopen_fraction((T) 0.9), surface_only(1), update_interval(200) {}

    /* Which declared phase, if any, occupies a voxel of this material number? */
    plint phaseOfMaterial (plint material) const {
        for (size_t k = 0; k < phases.size(); ++k)
            if (phases[k].material_number == material) return phases[k].id;
        return 0;
    }
    const SolidPhase* byId (plint pid) const {
        for (size_t k = 0; k < phases.size(); ++k)
            if (phases[k].id == pid) return &phases[k];
        return 0;
    }
};


/* ============================================================================
 *  initPhaseLattice3D
 *
 *  Seeds the phase lattice from the geometry at start-up: every voxel whose
 *  material number matches a declared phase gets that phase's id, everything
 *  else gets 0.
 *
 *  Lattices/fields: lattice = phaseLattice (written), field = geometry (read).
 * ============================================================================
 */
template<typename T1, template<typename U> class Descriptor, typename T2>
class initPhaseLattice3D : public BoxProcessingFunctional3D_LS<T1,Descriptor,T2>
{
public:
    initPhaseLattice3D (const std::vector<plint> &materials_, const std::vector<plint> &ids_)
        : materials(materials_), ids(ids_) {}

    virtual void process (Box3D domain, BlockLattice3D<T1,Descriptor> &lattice, ScalarField3D<T2> &field) {
        Dot3D offset = computeRelativeDisplacement(lattice, field);
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    const plint mat = (plint) field.get(iX+offset.x, iY+offset.y, iZ+offset.z);
                    T1 pid = (T1) 0;
                    for (size_t k = 0; k < materials.size(); ++k)
                        if (materials[k] == mat) { pid = (T1) ids[k]; break; }
                    Array<T1,7> g;
                    g[0]=(pid-(T1)1)/4;
                    g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(pid-(T1)1)/8;
                    lattice.get(iX,iY,iZ).setPopulations(g);
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulkAndEnvelope; }
    virtual initPhaseLattice3D<T1,Descriptor,T2>* clone() const {
        return new initPhaseLattice3D<T1,Descriptor,T2>(*this);
    }
    void getTypeOfModification (std::vector<modif::ModifT> &modified) const {
        modified[0] = modif::staticVariables;   // phase lattice written
        modified[1] = modif::nothing;           // geometry read only
    }
private:
    std::vector<plint> materials, ids;
};


/* ============================================================================
 *  surfaceDissolutionKinetics3D
 *
 *  Dissolves mineral at the WETTED surface of every declared phase, and puts
 *  the products into the open face neighbours.
 *
 *  For each voxel that carries a phase id:
 *      1. find its open (pore or biofilm) face neighbours -- that is the water
 *         actually touching the mineral.  No open neighbour => enclosed => no
 *         reaction.
 *      2. average the solute concentrations over those neighbours: this is the
 *         water the mineral sees.
 *      3. call defineDissolutionRate() with that water, the phase id and the
 *         remaining mineral.
 *      4. reduce the mineral inventory, and split the released products evenly
 *         among the same open neighbours.
 *
 *  Limiter: a step may not remove more than the mineral present, and may not
 *  remove more than MAX_DISSOL_FRACTION of it at once.  Everything released is
 *  scaled by the same factor as the mineral removed, so the element balance is
 *  exact whatever the limiter does.
 *
 *  Lattice layout   (length 2*subsNum + 2):
 *      [0            .. subsNum-1]   substrate concentrations   C
 *      [subsNum      .. 2*subsNum-1] substrate increments       dC
 *      [2*subsNum]                   mask lattice
 *      [2*subsNum+1]                 phase lattice
 *
 *  appliesTo() is bulk, NOT bulkAndEnvelope: this processor reads and WRITES
 *  face neighbours, and doing that from an envelope cell would run off the end
 *  of the block.  Same reasoning as surfaceAbioticKinetics3D.
 * ============================================================================
 */
template<typename T, template<typename U> class Descriptor>
class surfaceDissolutionKinetics3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    surfaceDissolutionKinetics3D (plint nx_, plint ny_, plint nz_, plint subsNum_, T dt_,
                                  plint solid_, plint bb_,
                                  const std::vector<plint> &pore_,
                                  const DissolutionConfig *cfg_)
        : nx(nx_), ny(ny_), nz(nz_), subsNum(subsNum_), dt(dt_),
          solid(solid_), bb(bb_), pore(pore_), cfg(cfg_),
          dCloc(subsNum_), maskLloc(2*subsNum_), phaseLloc(2*subsNum_ + 1)
    {}

    virtual void process (Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices)
    {
        Dot3D absoluteOffset = lattices[0]->getLocation();

        std::vector<Dot3D> off;
        off.reserve(phaseLloc + 1);
        for (plint iT = 0; iT <= phaseLloc; ++iT)
            off.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT]));

        /* the six face directions */
        static const plint dloc[6][3] = { {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1} };

        std::vector<T> water(subsNum, T());     // solutes averaged over the wetted neighbours
        std::vector<T> subsR(subsNum, T());     // what the user's rate law releases
        std::vector<plint> openDir;             // which of the six faces are open
        openDir.reserve(6);

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            const plint absX = iX + absoluteOffset.x;
            if (absX <= 0 || absX >= nx-1) continue;       // inlet / outlet planes

            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                const plint absY = iY + absoluteOffset.y;
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    const plint absZ = iZ + absoluteOffset.z;

                    /* ---- does this voxel carry a dissolvable phase? -------- */
                    const plint pid = util::roundToInt(
                        lattices[phaseLloc]->get(iX+off[phaseLloc].x, iY+off[phaseLloc].y,
                                                 iZ+off[phaseLloc].z).computeDensity());
                    if (pid <= 0) continue;
                    const SolidPhase *ph = cfg->byId(pid);
                    if (ph == 0 || ph->substrate < 0 || ph->substrate >= subsNum) continue;

                    /* ---- how much mineral is left? ------------------------- */
                    const plint sM = ph->substrate;
                    const T mineral = lattices[sM]->get(iX+off[sM].x, iY+off[sM].y,
                                                        iZ+off[sM].z).computeDensity();
                    if (!(mineral > thrd)) continue;

                    /* ---- which faces are open to water? -------------------- */
                    openDir.clear();
                    for (int d = 0; d < 6; ++d) {
                        const plint aX = absX + dloc[d][0];
                        const plint aY = absY + dloc[d][1];
                        const plint aZ = absZ + dloc[d][2];
                        if (aX < 1 || aX > nx-2 || aY < 0 || aY > ny-1 || aZ < 0 || aZ > nz-1) continue;

                        const plint jX = iX + dloc[d][0], jY = iY + dloc[d][1], jZ = iZ + dloc[d][2];
                        const plint nbrMask = util::roundToInt(
                            lattices[maskLloc]->get(jX+off[maskLloc].x, jY+off[maskLloc].y,
                                                    jZ+off[maskLloc].z).computeDensity());
                        if (nbrMask == solid || nbrMask == bb) continue;   // still rock
                        openDir.push_back((plint) d);
                    }
                    /* Fully enclosed: no water can reach it, so it does not
                     * dissolve.  With surface_only turned off we let it react
                     * anyway, which is unphysical but useful for testing. */
                    if (openDir.empty() && cfg->surface_only != 0) continue;

                    /* ---- the water in contact with the mineral ------------- */
                    const plint nOpen = (plint) openDir.size();
                    for (plint iS = 0; iS < subsNum; ++iS) water[iS] = T();
                    if (nOpen > 0) {
                        for (plint k = 0; k < nOpen; ++k) {
                            const int d = (int) openDir[k];
                            const plint jX = iX + dloc[d][0], jY = iY + dloc[d][1], jZ = iZ + dloc[d][2];
                            for (plint iS = 0; iS < subsNum; ++iS) {
                                water[iS] += lattices[iS]->get(jX+off[iS].x, jY+off[iS].y,
                                                               jZ+off[iS].z).computeDensity();
                            }
                        }
                        for (plint iS = 0; iS < subsNum; ++iS) water[iS] /= (T) nOpen;
                    }
                    for (plint iS = 0; iS < subsNum; ++iS) if (water[iS] < T()) water[iS] = T();

                    /* ---- ask the user's rate law -------------------------- */
                    const plint mask = util::roundToInt(
                        lattices[maskLloc]->get(iX+off[maskLloc].x, iY+off[maskLloc].y,
                                                iZ+off[maskLloc].z).computeDensity());
                    for (plint iS = 0; iS < subsNum; ++iS) subsR[iS] = T();
                    T mineralR = T();
                    defineDissolutionRate(pid, water, (double) mineral, subsR, mineralR, mask);

                    /* mineralR is a RATE (mol/L/s) and must be negative or zero:
                     * dissolution removes mineral.  A positive value would be
                     * precipitation, which is the other file's job. */
                    if (!(mineralR < T())) continue;

                    T dMineral = mineralR * dt;                       // negative
                    /* Never remove more than is there, and never more than a
                     * fixed fraction in one step.  Anything that scales the
                     * mineral scales the products by the SAME factor, so the
                     * element balance survives the limiter exactly. */
                    const T MAX_DISSOL_FRACTION = (T) 0.5;
                    T scale = (T) 1;
                    const T maxRemove = MAX_DISSOL_FRACTION * mineral;
                    if (-dMineral > maxRemove) scale = maxRemove / (-dMineral);
                    dMineral *= scale;

                    if (!(-dMineral > thrd)) continue;

                    /* ---- take the mineral out -----------------------------
                     * [FIX] This used to accumulate into dC[sM] and leave the apply to
                     * update_abiotic_rxnLattices.  That step skips any voxel whose mask is solid
                     * or bounce-back -- and a dissolving mineral voxel is solid BY CONSTRUCTION,
                     * which is the whole reason it carries a phase id.  So the decrement went
                     * into dC, was never applied, and was wiped when dC was reset at the top of
                     * the next step.  The mineral inventory never fell.
                     *
                     * The products meanwhile go into the PORE neighbours, which the apply step
                     * does visit, so they WERE applied.  Net effect: solute created from nothing
                     * every step, forever, no mineral consumed, and nothing ever reaching the
                     * reopen threshold.
                     *
                     * The mineral lattice is written DIRECTLY here instead.  Safe, and only for
                     * this field: a phase substrate must be declared <immobile>, so it is never
                     * collided and never streamed -- storage, not transport. */
                    {
                        Array<T,7> g;
                        Cell<T,Descriptor> &cell = lattices[sM]->get(
                            iX+off[sM].x, iY+off[sM].y, iZ+off[sM].z);
                        cell.getPopulations(g);
                        g[0]+=dMineral/4; g[1]+=dMineral/8; g[2]+=dMineral/8; g[3]+=dMineral/8;
                        g[4]+=dMineral/8; g[5]+=dMineral/8; g[6]+=dMineral/8;
                        cell.setPopulations(g);
                    }

                    /* ---- put the products into the open neighbours --------- */
                    if (nOpen > 0) {
                        for (plint iS = 0; iS < subsNum; ++iS) {
                            if (iS == sM) continue;                  // the mineral itself, already handled
                            const T released = subsR[iS] * dt * scale;
                            if (std::fabs(released) <= thrd) continue;
                            const T perFace = released / (T) nOpen;

                            for (plint k = 0; k < nOpen; ++k) {
                                const int d = (int) openDir[k];
                                const plint jX = iX + dloc[d][0], jY = iY + dloc[d][1], jZ = iZ + dloc[d][2];
                                Array<T,7> g;
                                Cell<T,Descriptor> &cell = lattices[iS+dCloc]->get(
                                    jX+off[iS+dCloc].x, jY+off[iS+dCloc].y, jZ+off[iS+dCloc].z);
                                cell.getPopulations(g);
                                g[0]+=perFace/4; g[1]+=perFace/8; g[2]+=perFace/8; g[3]+=perFace/8;
                                g[4]+=perFace/8; g[5]+=perFace/8; g[6]+=perFace/8;
                                cell.setPopulations(g);
                            }
                        }
                    }
                }
            }
        }
    }

    /* bulk only: this processor reads AND writes face neighbours */
    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }

    virtual surfaceDissolutionKinetics3D<T,Descriptor>* clone() const {
        return new surfaceDissolutionKinetics3D<T,Descriptor>(*this);
    }
    void getTypeOfModification (std::vector<modif::ModifT> &modified) const {
        for (size_t i = 0; i < modified.size(); ++i) modified[i] = modif::nothing;
        for (plint iT = dCloc; iT < maskLloc; ++iT) modified[iT] = modif::staticVariables;
        /* [FIX] The mineral inventory is now written straight into its own substrate lattice
         * rather than into dC, so those lattices are modified too.  Only the phase substrates
         * are ever touched, so only those are marked. */
        if (cfg != 0) {
            for (size_t k = 0; k < cfg->phases.size(); ++k) {
                const plint sM = cfg->phases[k].substrate;
                if (sM >= 0 && sM < subsNum) modified[(size_t) sM] = modif::staticVariables;
            }
        }
    }

private:
    plint nx, ny, nz, subsNum;
    T dt;
    plint solid, bb;
    std::vector<plint> pore;
    const DissolutionConfig *cfg;
    plint dCloc, maskLloc, phaseLloc;
};


/* ============================================================================
 *  dissolNodeConversion3D
 *
 *  The mirror of precipNodeConversion3D.  When a phase-carrying voxel's mineral
 *  inventory falls below  reopen_fraction * full_density, the voxel becomes
 *  pore again and its phase id is cleared.
 *
 *  The reopen threshold is deliberately LOWER than the seal threshold; see the
 *  note on hysteresis at the top of this file.
 *
 *  Lattices: [ mineral substrate 0..subsNum-1 , mask , phase ]  -- the same
 *  leading block as the dissolution processor so the caller can reuse the
 *  vector, minus the dC block.
 *
 *  Layout (length subsNum + 2):
 *      [0 .. subsNum-1]  substrate concentrations (only the phase substrates
 *                        are read)
 *      [subsNum]         mask lattice   (written)
 *      [subsNum+1]       phase lattice  (written)
 * ============================================================================
 */
template<typename T, template<typename U> class Descriptor>
class dissolNodeConversion3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    dissolNodeConversion3D (plint subsNum_, plint poreCode_, plint solid_, plint bb_,
                            const DissolutionConfig *cfg_)
        : subsNum(subsNum_), poreCode(poreCode_), solid(solid_), bb(bb_), cfg(cfg_),
          maskLloc(subsNum_), phaseLloc(subsNum_ + 1)
    {}

    virtual void process (Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices)
    {
        std::vector<Dot3D> off;
        off.reserve(phaseLloc + 1);
        for (plint iT = 0; iT <= phaseLloc; ++iT)
            off.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT]));

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {

                    const plint pid = util::roundToInt(
                        lattices[phaseLloc]->get(iX+off[phaseLloc].x, iY+off[phaseLloc].y,
                                                 iZ+off[phaseLloc].z).computeDensity());
                    if (pid <= 0) continue;
                    const SolidPhase *ph = cfg->byId(pid);
                    if (ph == 0 || ph->substrate < 0 || ph->substrate >= subsNum) continue;
                    if (!(ph->full_density > T())) continue;

                    const plint sM = ph->substrate;
                    const T mineral = lattices[sM]->get(iX+off[sM].x, iY+off[sM].y,
                                                        iZ+off[sM].z).computeDensity();

                    if (mineral < cfg->reopen_fraction * ph->full_density) {
                        /* open the voxel back up ... */
                        Array<T,7> g;
                        g[0]=(T)(poreCode-1)/4;
                        g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T)(poreCode-1)/8;
                        lattices[maskLloc]->get(iX+off[maskLloc].x, iY+off[maskLloc].y,
                                                iZ+off[maskLloc].z).setPopulations(g);
                        /* ... and clear its phase, so it is an ordinary pore voxel
                         * again.  If mineral later re-precipitates there,
                         * precipNodeConversion3D stamps the phase back on. */
                        Array<T,7> h;
                        h[0]=(T)(0-1)/4;
                        h[1]=h[2]=h[3]=h[4]=h[5]=h[6]=(T)(0-1)/8;
                        lattices[phaseLloc]->get(iX+off[phaseLloc].x, iY+off[phaseLloc].y,
                                                 iZ+off[phaseLloc].z).setPopulations(h);
                    }
                }
            }
        }
    }

    virtual BlockDomain::DomainT appliesTo() const { return BlockDomain::bulk; }
    virtual dissolNodeConversion3D<T,Descriptor>* clone() const {
        return new dissolNodeConversion3D<T,Descriptor>(*this);
    }
    void getTypeOfModification (std::vector<modif::ModifT> &modified) const {
        for (size_t i = 0; i < modified.size(); ++i) modified[i] = modif::nothing;
        modified[maskLloc]  = modif::staticVariables;
        modified[phaseLloc] = modif::staticVariables;
    }

private:
    plint subsNum, poreCode, solid, bb;
    const DissolutionConfig *cfg;
    plint maskLloc, phaseLloc;
};

#endif // DISSOLUTION_VOP_HH
