/* This file is a part of the CompLaB program.
 *
 * The CompLaB3D software (3-D pore-scale extension) is developed since 2024
 * by the University of Georgia (United States, Meile Lab, Department of
 * Marine Sciences). The original 2-D CompLaB v1.0 was a collaboration of
 * the University of Georgia and Chungnam National University (South Korea).
 * 
 * The CompLaB3D extension (3-D pore-scale model with biotic/abiotic
 * kinetics, Anderson-Accelerated Newton-Raphson PCF equilibrium chemistry,
 * and the CompLaB Studio GUI) is developed by Shahram Asgari and
 * Christof Meile (Meile Lab) at the University of Georgia. The 2-D
 * predecessor (CompLaB v1.0) was developed by Heewon Jung et al.
 *
 * Contact:
 * Shahram Asgari
 * Department of Marine Sciences (Meile Lab)
 * University of Georgia
 * Athens, GA 30602, USA
 * shahram.asgari@uga.edu
 *
 * Christof Meile
 * Department of Marine Sciences
 * University of Georgia
 * Athens, GA 30602, USA
 *
 * The most recent release of CompLaB can be downloaded at 
 * https://bitbucket.org/MeileLab/complab/downloads/
 *
 * CompLaB is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef COMPLAB3D_PROCESSORS_PART2_HH
#define COMPLAB3D_PROCESSORS_PART2_HH

#include <random>
using namespace plb;
typedef double T;

/* ===============================================================================================================
   ======================================= MASK AND DYNAMICS UPDATE PROCESSORS ===================================
   =============================================================================================================== */

/* -----------------------------------------------------------------------------------------------------------------
 * updateLocalMaskNtotalLattices3D — Voxel Reclassification (Appendix B, Algorithm B.5)
 *
 * After the redistribution CA (Algorithms B.2/B.3) has finished and every voxel is within
 * capacity, this processor updates two things at every non-solid, non-boundary voxel:
 *   1. The total biomass density lattice (sum of all species at the voxel).
 *   2. The mask (voxel identity) lattice, switching voxels between pore and biofilm
 *      based on the reclassification threshold:
 *          thrdbMassRho = thrd_bFilmFrac × max_bMassRho        (Eq. in Appendix B.5)
 *
 *      - Biofilm → pore:  if total biomass < thrdbMassRho      (biofilm has eroded away)
 *      - Pore → biofilm:  if total biomass ≥ thrdbMassRho      (new biofilm has grown)
 *
 * When a biofilm voxel is reclassified to pore, its new pore-type label is copied from the
 * nearest face-connected pore neighbor so that pore labels stay spatially consistent.
 * When a pore voxel becomes biofilm, the new mask is the sum of species identifiers present.
 * ----------------------------------------------------------------------------------------------------------------- */
template<typename T, template<typename U> class Descriptor>
class updateLocalMaskNtotalLattices3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    /* Constructor parameters:
     *   nx_, ny_, nz_    — global domain dimensions (used for boundary checks)
     *   length_          — total number of lattices passed to process()
     *   bb_              — mask value for domain boundary (bounce-back) voxels
     *   solid_           — mask value for solid-grain voxels
     *   bio_             — vector of [speciesID, ...] for each microbial species
     *   pore_            — vector of all valid pore-type mask values
     *   bMassFrac_       — thrd_bFilmFrac (fraction, 0–1), from XML config
     *   maxbMassRho_     — max_bMassRho, from XML config
     *   thrdbMassRho is precomputed as bMassFrac_ × maxbMassRho_ (reclassification threshold) */
    updateLocalMaskNtotalLattices3D(plint nx_, plint ny_, plint nz_, plint length_, plint bb_, plint solid_, std::vector<std::vector<plint>> bio_, std::vector<plint> pore_, T bMassFrac_, T maxbMassRho_)
    : nx(nx_), ny(ny_), nz(nz_), length(length_), bb(bb_), solid(solid_), bio(bio_), pore(pore_), thrdbMassRho(bMassFrac_*maxbMassRho_)
    {}
    /* Lattice layout (D3Q7 populations encode a scalar density):
     *   lattices[0 .. #ofbM-1]     = original biomass lattices (one per microbial species)
     *   lattices[#ofbM .. len-4]   = copy biomass lattices (read-only snapshots, unused here)
     *   lattices[len-3]            = total biomass lattice (sum of all species)
     *   lattices[len-2]            = mask lattice (voxel identity: pore / biofilm / solid / bb)
     *   lattices[len-1]            = age lattice (distance-to-pore, not modified here)         */
    // lattices[0~(#ofbM-1)] = original biomass lattices
    // lattices[#ofbM~(len-3)] = copy biomass lattices
    // lattices[len-3] = total biomass lattice
    // lattices[len-2] = mask lattice
    // lattices[len-1] = age lattice
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        plint maskLloc = length-2, bMtLloc = length-3;  // shorthand indices for mask and total-biomass lattices
        std::vector<Dot3D> vec_offset;
        Dot3D absoluteOffset = lattices[0]->getLocation();  // offset from local to global coordinates
        // Precompute relative displacements between lattice[0] and every other lattice
        // (needed because Palabos may store lattices in different memory blocks with offsets)
        for (plint iL=0; iL<length; ++iL) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0],*lattices[iL])); }
        // ---- Triple loop over every voxel in the local subdomain ----
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint iXm = iX+vec_offset[maskLloc].x;  // mask-lattice x-index for this voxel
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint iYm = iY+vec_offset[maskLloc].y;
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint iZm = iZ+vec_offset[maskLloc].z;
                    // Read the current mask (voxel identity) — stored as D3Q7 density
                    plint mask = util::roundToInt( lattices[maskLloc]->get(iXm,iYm,iZm).computeDensity() );
                    // Skip domain-boundary and solid-grain voxels (they never hold biomass)
                    if (mask != bb && mask != solid) {
                        T bmass = 0.; plint newmask = 0;  // bmass = running total; newmask = species-based label
                        // Read the previously stored total biomass at this voxel
                        T bMt = lattices[bMtLloc]->get(iX+vec_offset[bMtLloc].x,iY+vec_offset[bMtLloc].y,iZ+vec_offset[bMtLloc].z).computeDensity();
                        // Sum biomass across all microbial species at this voxel
                        for (size_t iM=0; iM<bio.size(); ++iM) {
                            T bm = lattices[iM]->get(iX+vec_offset[iM].x,iY+vec_offset[iM].y,iZ+vec_offset[iM].z).computeDensity();
                            // Only count species with non-negligible biomass (above numerical threshold)
                            if ( bm > COMPLAB_THRD) { bmass += bm; newmask += bio[iM][0];}
                        }
                        // --- Step 1: Update the total-biomass lattice if the sum has changed ---
                        // (After redistribution, individual species densities may have shifted)
                        if ( std::abs(bmass-bMt) > COMPLAB_THRD) {
                            // Encode the scalar 'bmass' into D3Q7 populations:
                            // density = g[0] + g[1..6] = (bmass-1)/4 + 6×(bmass-1)/8 = bmass
                            Array<T,7> g;
                            g[0]=(T) (bmass-1)/4; g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T) (bmass-1)/8;
                            plint iXt = iX+vec_offset[bMtLloc].x, iYt = iY+vec_offset[bMtLloc].y, iZt = iZ+vec_offset[bMtLloc].z;
                            lattices[bMtLloc]->get(iXt,iYt,iZt).setPopulations(g);
                        }
                        // --- Step 2: Check if this voxel needs to switch identity (Algorithm B.5) ---
                        // Determine whether the current mask is a pore type
                        bool poreflag = 0;
                        for (size_t iP=0; iP<pore.size(); ++iP) {
                            if (mask == pore[iP]) { poreflag = 1; break; }
                        }
                        // Case A: Biofilm → Pore  (total biomass dropped below threshold)
                        if (poreflag == 0 && bmass < thrdbMassRho) {
                            // Default to the first pore type; will be overwritten if a neighboring pore exists
                            newmask = pore[0];
                            // If multiple pore types exist, copy the label from the nearest pore neighbor
                            // to keep pore labels spatially consistent (Algorithm B.5)
                            if (pore.size()>1) {
                                plint absX = iX + absoluteOffset.x, absY = iY + absoluteOffset.y, absZ = iZ + absoluteOffset.z;
                                plint nbrs = 0;
                                std::vector<plint> delXYZ;

                                // Determine neighbors based on boundary position (3D version)
                                // Build list of up to 6 face-connected neighbors (D3Q7 lattice: ±x, ±y, ±z)
                                // excluding directions that would go outside the domain
                                bool atXmin = (absX == 0);
                                bool atXmax = (absX == nx-1);
                                bool atYmin = (absY == 0);
                                bool atYmax = (absY == ny-1);
                                bool atZmin = (absZ == 0);
                                bool atZmax = (absZ == nz-1);

                                // Add valid neighbors (each neighbor stored as 3 consecutive entries: dx, dy, dz)
                                if (!atXmax) { delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); ++nbrs; }
                                if (!atXmin) { delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); ++nbrs; }
                                if (!atYmax) { delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); ++nbrs; }
                                if (!atYmin) { delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0); ++nbrs; }
                                if (!atZmax) { delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); ++nbrs; }
                                if (!atZmin) { delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1); ++nbrs; }

                                // Search face neighbors for the first pore-type voxel; adopt its mask label
                                for (plint iT=0; iT<nbrs; ++iT) {
                                    plint delx = delXYZ[iT*3], dely = delXYZ[iT*3+1], delz = delXYZ[iT*3+2];
                                    plint nbrmask = util::roundToInt( lattices[maskLloc]->get(iXm+delx,iYm+dely,iZm+delz).computeDensity() );
                                    if (nbrmask != bb && nbrmask != solid) {
                                        bool nbrporeflag = 0;
                                        for (size_t iP=0; iP<pore.size(); ++iP) {
                                            if (nbrmask == pore[iP]) { nbrporeflag = 1; break; }
                                        }
                                        // Found a neighboring pore voxel — copy its pore-type label
                                        if (nbrporeflag == 1) { newmask = nbrmask; break; }
                                    }
                                }
                            }
                            // Write the new pore mask into the mask lattice (encoded as D3Q7 density)
                            Array<T,7> g;
                            g[0]=(T) (newmask-1)/4; g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T) (newmask-1)/8;
                            lattices[maskLloc]->get(iXm,iYm,iZm).setPopulations(g);
                        }
                        // Case B: Pore → Biofilm  (total biomass reached or exceeded threshold)
                        else if (poreflag == 1 && bmass >= thrdbMassRho) {
                            // newmask was accumulated from species IDs above — identifies the biofilm type
                            if (newmask > 0) {
                                // Write the new biofilm mask into the mask lattice
                                Array<T,7> g;
                                g[0]=(T) (newmask-1)/4; g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T) (newmask-1)/8;
                                lattices[maskLloc]->get(iXm,iYm,iZm).setPopulations(g);
                            }
                            else {
                                // Safety check: should never happen (biomass ≥ threshold but no species found)
                                std::cout << "Error: Updating mask failed.\n";
                                /* [ROBUST BUG5] do not abort the MPI job here; leave the voxel identity unchanged and continue */ ;
                            }
                        }
                        // If neither case applies, the voxel keeps its current identity (no action needed)
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulk;  // operates on interior voxels only; envelope handled by MPI exchange
    }
    virtual updateLocalMaskNtotalLattices3D<T,Descriptor>* clone() const {
        return new updateLocalMaskNtotalLattices3D<T,Descriptor>(*this);
    }
    /* Declare which lattices are modified:
     *   biomass lattices        — not modified (read only)
     *   total biomass (len-3)   — modified (staticVariables: populations updated)
     *   mask (len-2)            — modified (staticVariables: populations updated)
     *   age (len-1)             — not modified                                    */
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iL=0; iL<length-3; ++iL) {
            modified[iL] = modif::nothing;
        }
        modified[length-3] = modif::staticVariables;
        modified[length-2] = modif::staticVariables;
        modified[length-1] = modif::nothing;
    }
private:
    plint nx, ny, nz, length, bb, solid;
    std::vector<std::vector<plint>> bio;   // bio[species][0] = species mask ID
    std::vector<plint> pore;               // list of valid pore-type mask values
    T thrdbMassRho;                        // reclassification threshold = thrd_bFilmFrac × max_bMassRho
};

/* -----------------------------------------------------------------------------------------------------------------
 * fdDiffusion3D — Finite-Difference Biomass Diffusion (3D Laplacian)
 *
 * Diffuses planktonic (suspended) biomass through the pore space using an explicit
 * finite-difference approximation of the 3D Laplacian:
 *
 *     B_new = B + ν × [ (B_x+ − 2B + B_x−) + (B_y+ − 2B + B_y−) + (B_z+ − 2B + B_z−) ]
 *
 * where ν is the diffusion coefficient (in lattice units) and B_x+, B_x−, etc. are the
 * biomass densities at the six face-connected neighbors.
 *
 * Boundary conditions:
 *   - Domain edges and solid/boundary interfaces use zero-flux (mirror/Neumann) BC:
 *     the missing neighbor value is replaced by the center value (b0), so the second
 *     derivative in that direction becomes zero and no flux crosses the boundary.
 *   - The bdryGap parameter excludes a strip of voxels at the x-boundaries (inlet/outlet).
 *
 * The processor reads from COPY lattices (snapshot from before diffusion) and writes
 * to ORIGINAL lattices, avoiding read-after-write conflicts within a single pass.
 * ----------------------------------------------------------------------------------------------------------------- */
template<typename T, template<typename U> class Descriptor>
class fdDiffusion3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    /* Constructor parameters:
     *   nx_, ny_, nz_  — global domain dimensions
     *   length_        — total number of lattices: 2×numbM originals/copies + 1 mask
     *   bdryGap_       — number of voxels to skip at x-min and x-max (inlet/outlet buffer)
     *   nu_            — diffusion coefficient in lattice units (Δt/Δx²-scaled)            */
    fdDiffusion3D(plint nx_, plint ny_, plint nz_, plint length_, plint bdryGap_, T nu_)
    : nx(nx_), ny(ny_), nz(nz_), length(length_), bdryGap(bdryGap_), nu(nu_)
    {}
    /* Lattice layout:
     *   lattices[0 .. numbM-1]         = ORIGINAL biomass lattices (write targets)
     *   lattices[numbM .. 2×numbM-1]   = COPY biomass lattices (read-only snapshots)
     *   lattices[length-1]             = mask lattice (read-only, for boundary detection) */
    // lattices[0~(#ofbM-1)] = original biomass lattices
    // lattices[#ofbM~(len-2)] = copy biomass lattices
    // lattices[len-1] = mask lattice
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        std::vector<Dot3D> vec_offset;
        plint maskLloc = length-1, numbM = (length-1)/2;  // numbM = number of biomass species
        for (plint iL=0; iL<length; ++iL) {
            vec_offset.push_back(computeRelativeDisplacement(*lattices[0],*lattices[iL]));
        }
        Dot3D absoluteOffset = lattices[0]->getLocation();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;
            // Skip voxels within the inlet/outlet boundary gap (no diffusion there)
            if ( absX >= bdryGap && absX <= nx-1-bdryGap ) {
                for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                    plint absY = iY + absoluteOffset.y;
                    for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                        plint absZ = iZ + absoluteOffset.z;
                        // Read the mask to determine if this voxel is pore or biofilm (mask > 1)
                        // mask ≤ 1 means solid or boundary — skip (no biomass lives there)
                        plint iXm = iX + vec_offset[maskLloc].x, iYm = iY + vec_offset[maskLloc].y, iZm = iZ + vec_offset[maskLloc].z;
                        plint mask = util::roundToInt( lattices[maskLloc]->get(iXm,iYm,iZm).computeDensity() );
                        if (mask > 1) {
                            // --- Collect neighbor mask values for boundary detection ---
                            // mxp/mxn = mask of +x/−x neighbor, etc. (0 if at domain edge)
                            // A neighbor mask < 2 means solid or boundary → apply zero-flux BC
                            plint mxp=0,mxn=0,myp=0,myn=0,mzp=0,mzn=0;
                            // x-direction: handle domain edges and interior separately
                            if ( absX == bdryGap ) { mxp=lattices[maskLloc]->get(iXm+1,iYm,iZm).computeDensity(); }
                            else if ( absX == nx-1-bdryGap ) { mxn=lattices[maskLloc]->get(iXm-1,iYm,iZm).computeDensity(); }
                            else { mxp=lattices[maskLloc]->get(iXm+1,iYm,iZm).computeDensity(); mxn=lattices[maskLloc]->get(iXm-1,iYm,iZm).computeDensity(); }
                            // y-direction
                            if ( absY == 0 ) { myp=lattices[maskLloc]->get(iXm,iYm+1,iZm).computeDensity(); }
                            else if ( absY == ny-1 ) { myn=lattices[maskLloc]->get(iXm,iYm-1,iZm).computeDensity(); }
                            else { myp=lattices[maskLloc]->get(iXm,iYm+1,iZm).computeDensity(); myn=lattices[maskLloc]->get(iXm,iYm-1,iZm).computeDensity(); }
                            // z-direction
                            if ( absZ == 0 ) { mzp=lattices[maskLloc]->get(iXm,iYm,iZm+1).computeDensity(); }
                            else if ( absZ == nz-1 ) { mzn=lattices[maskLloc]->get(iXm,iYm,iZm-1).computeDensity(); }
                            else { mzp=lattices[maskLloc]->get(iXm,iYm,iZm+1).computeDensity(); mzn=lattices[maskLloc]->get(iXm,iYm,iZm-1).computeDensity(); }

                            // --- Collect biomass densities from COPY lattices ---
                            // bxp/bxn = biomass at +x/−x neighbor; b0 = biomass at center
                            // All reads come from copy lattices [numbM .. 2*numbM-1] to avoid
                            // read-after-write conflicts (writes go to originals [0 .. numbM-1])
                            std::vector<T> bxp(numbM,0), bxn(numbM,0), byp(numbM,0), byn(numbM,0), bzp(numbM,0), bzn(numbM,0), b0(numbM,0);
                            // Read center value for each species from the copy lattice
                            for (plint iM=0; iM<numbM; ++iM) {
                                plint iXb = iX + vec_offset[iM+numbM].x, iYb = iY + vec_offset[iM+numbM].y, iZb = iZ + vec_offset[iM+numbM].z;
                                b0[iM] = lattices[iM+numbM]->get(iXb,iYb,iZb).computeDensity();
                            }
                            // x-direction neighbor biomass with zero-flux (mirror) BC:
                            // At domain edge or if neighbor is solid (mask < 2), set missing
                            // neighbor = center value (b0), making d²f/dx² = 0 in that direction
                            if ( absX == bdryGap || mxn < 2) {
                                for (plint iM=0; iM<numbM; ++iM) {
                                    plint iXb = iX + vec_offset[iM+numbM].x, iYb = iY + vec_offset[iM+numbM].y, iZb = iZ + vec_offset[iM+numbM].z;
                                    bxp[iM] = lattices[iM+numbM]->get(iXb+1,iYb,iZb).computeDensity();
                                    bxn[iM] = b0[iM];  // mirror BC: −x neighbor = center
                                }
                            }
                            else if ( absX == nx-1-bdryGap || mxp < 2 ) {
                                for (plint iM=0; iM<numbM; ++iM) {
                                    plint iXb = iX + vec_offset[iM+numbM].x, iYb = iY + vec_offset[iM+numbM].y, iZb = iZ + vec_offset[iM+numbM].z;
                                    bxp[iM] = b0[iM];  // mirror BC: +x neighbor = center
                                    bxn[iM] = lattices[iM+numbM]->get(iXb-1,iYb,iZb).computeDensity();
                                }
                            }
                            else {
                                // Interior: read both +x and −x neighbors normally
                                for (plint iM=0; iM<numbM; ++iM) {
                                    plint iXb = iX + vec_offset[iM+numbM].x, iYb = iY + vec_offset[iM+numbM].y, iZb = iZ + vec_offset[iM+numbM].z;
                                    bxp[iM] = lattices[iM+numbM]->get(iXb+1,iYb,iZb).computeDensity();
                                    bxn[iM] = lattices[iM+numbM]->get(iXb-1,iYb,iZb).computeDensity();
                                }
                            }
                            // y-direction neighbor biomass (same zero-flux BC logic)
                            if ( absY == 0 || myn < 2 ) {
                                for (plint iM=0; iM<numbM; ++iM) {
                                    plint iXb = iX + vec_offset[iM+numbM].x, iYb = iY + vec_offset[iM+numbM].y, iZb = iZ + vec_offset[iM+numbM].z;
                                    byp[iM] = lattices[iM+numbM]->get(iXb,iYb+1,iZb).computeDensity();
                                    byn[iM] = b0[iM];  // mirror BC
                                }
                            }
                            else if ( absY == ny-1 || myp < 2 ) {
                                for (plint iM=0; iM<numbM; ++iM) {
                                    plint iXb = iX + vec_offset[iM+numbM].x, iYb = iY + vec_offset[iM+numbM].y, iZb = iZ + vec_offset[iM+numbM].z;
                                    byp[iM] = b0[iM];  // mirror BC
                                    byn[iM] = lattices[iM+numbM]->get(iXb,iYb-1,iZb).computeDensity();
                                }
                            }
                            else {
                                for (plint iM=0; iM<numbM; ++iM) {
                                    plint iXb = iX + vec_offset[iM+numbM].x, iYb = iY + vec_offset[iM+numbM].y, iZb = iZ + vec_offset[iM+numbM].z;
                                    byp[iM] = lattices[iM+numbM]->get(iXb,iYb+1,iZb).computeDensity();
                                    byn[iM] = lattices[iM+numbM]->get(iXb,iYb-1,iZb).computeDensity();
                                }
                            }
                            // z-direction neighbor biomass (same zero-flux BC logic)
                            if ( absZ == 0 || mzn < 2 ) {
                                for (plint iM=0; iM<numbM; ++iM) {
                                    plint iXb = iX + vec_offset[iM+numbM].x, iYb = iY + vec_offset[iM+numbM].y, iZb = iZ + vec_offset[iM+numbM].z;
                                    bzp[iM] = lattices[iM+numbM]->get(iXb,iYb,iZb+1).computeDensity();
                                    bzn[iM] = b0[iM];  // mirror BC
                                }
                            }
                            else if ( absZ == nz-1 || mzp < 2 ) {
                                for (plint iM=0; iM<numbM; ++iM) {
                                    plint iXb = iX + vec_offset[iM+numbM].x, iYb = iY + vec_offset[iM+numbM].y, iZb = iZ + vec_offset[iM+numbM].z;
                                    bzp[iM] = b0[iM];  // mirror BC
                                    bzn[iM] = lattices[iM+numbM]->get(iXb,iYb,iZb-1).computeDensity();
                                }
                            }
                            else {
                                for (plint iM=0; iM<numbM; ++iM) {
                                    plint iXb = iX + vec_offset[iM+numbM].x, iYb = iY + vec_offset[iM+numbM].y, iZb = iZ + vec_offset[iM+numbM].z;
                                    bzp[iM] = lattices[iM+numbM]->get(iXb,iYb,iZb+1).computeDensity();
                                    bzn[iM] = lattices[iM+numbM]->get(iXb,iYb,iZb-1).computeDensity();
                                }
                            }

                            // --- Apply the 3D Laplacian for each biomass species ---
                            for (plint iM=0; iM<numbM; ++iM) {
                                // Clamp negligible values to zero to prevent numerical noise propagation
                                if (b0[iM] < COMPLAB_THRD ) b0[iM]=0;
                                if (bxp[iM] < COMPLAB_THRD ) bxp[iM]=0;
                                if (bxn[iM] < COMPLAB_THRD ) bxn[iM]=0;
                                if (byp[iM] < COMPLAB_THRD ) byp[iM]=0;
                                if (byn[iM] < COMPLAB_THRD ) byn[iM]=0;
                                if (bzp[iM] < COMPLAB_THRD ) bzp[iM]=0;
                                if (bzn[iM] < COMPLAB_THRD ) bzn[iM]=0;
                                // Explicit forward-Euler step: B_new = B + ν × ∇²B
                                // where ∇²B = (B_x+ − 2B + B_x−) + (B_y+ − 2B + B_y−) + (B_z+ − 2B + B_z−)
                                T new_bM = b0[iM] + nu*( (bxp[iM]-2*b0[iM]+bxn[iM]) + (byp[iM]-2*b0[iM]+byn[iM]) + (bzp[iM]-2*b0[iM]+bzn[iM]) );
                                // Only write back if result is above threshold (avoid storing numerical noise)
                                if (new_bM > COMPLAB_THRD ) {
                                    Array<T,7> g;
                                    // Write to the ORIGINAL lattice (index iM), not the copy
                                    plint iXb = iX + vec_offset[iM].x, iYb = iY + vec_offset[iM].y, iZb = iZ + vec_offset[iM].z;
                                    // Encode scalar new_bM into D3Q7 populations
                                    g[0]=(T) (new_bM-1)/4; g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T) (new_bM-1)/8;
                                    lattices[iM]->get(iXb,iYb,iZb).setPopulations(g);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulk;
    }
    virtual fdDiffusion3D<T,Descriptor>* clone() const {
        return new fdDiffusion3D<T,Descriptor>(*this);
    }
    /* Modification flags:
     *   original biomass lattices [0..numbM-1]  — modified (new diffused values written)
     *   copy biomass lattices [numbM..2*numbM-1] — not modified (read-only snapshot)
     *   mask lattice [length-1]                  — not modified (read-only)             */
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        plint numbM = (length-1)/2;
        for (plint iB=0; iB<numbM; ++iB) {
            modified[iB] = modif::staticVariables;
            modified[iB+numbM] = modif::nothing;
        }
        modified[length-1] = modif::nothing;
    }
private:
    plint nx, ny, nz, length, bdryGap;  // bdryGap = inlet/outlet buffer width in voxels
    T nu;                                // diffusion coefficient in lattice units
};

/* -----------------------------------------------------------------------------------------------------------------
 * updateSoluteDynamics3D — Solute LBM Relaxation Rate Update
 *
 * When voxel reclassification (Algorithm B.5) changes a voxel from pore to biofilm or
 * vice versa, the effective diffusivity of dissolved substrates (nutrients, electron
 * acceptors, etc.) changes. This processor synchronizes the LBM relaxation rate (omega)
 * on each solute lattice with the current voxel identity:
 *
 *   - Pore → Biofilm:  omega switches from substrOMEGAinPore to substrOMEGAinbMass
 *                       (solute diffuses more slowly through the biofilm matrix)
 *   - Biofilm → Pore:  omega switches from substrOMEGAinbMass to substrOMEGAinPore
 *                       (solute diffuses freely in open pore space)
 *
 * The omega values are precomputed from the user-specified diffusivities via the
 * standard LBM relation: D = (1/omega − 0.5) × cs² × Δt.
 * ----------------------------------------------------------------------------------------------------------------- */
template<typename T, template<typename U> class Descriptor>
class updateSoluteDynamics3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    /* Constructor parameters:
     *   subsNum_           — number of dissolved substrate species
     *   bb_, solid_        — mask values for domain boundary and solid grains
     *   pore_              — list of valid pore-type mask values
     *   substrOMEGAinbMass_ — omega for each substrate inside biofilm (lower diffusivity)
     *   substrOMEGAinPore_  — omega for each substrate in open pore (higher diffusivity) */
    updateSoluteDynamics3D(plint subsNum_, plint bb_, plint solid_, std::vector<plint> pore_, std::vector<T> substrOMEGAinbMass_, std::vector<T> substrOMEGAinPore_)
    : subsNum(subsNum_), bb(bb_), solid(solid_), pore(pore_), substrOMEGAinbMass(substrOMEGAinbMass_), substrOMEGAinPore(substrOMEGAinPore_)
    {}
    /* Lattice layout:
     *   lattices[0 .. subsNum-1]  = solute (substrate) LBM lattices
     *   lattices[subsNum]         = mask lattice (read-only)          */
    // lattices[0~(#ofSubs-1)] = substrate lattices
    // lattices[#ofSubs] = mask lattice
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        std::vector<Dot3D> vec_offset;
        for (plint iT=0; iT<subsNum+1; ++iT) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0],*lattices[iT])); }
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint iXm = iX+vec_offset[subsNum].x;
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint iYm = iY+vec_offset[subsNum].y;
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint iZm = iZ+vec_offset[subsNum].z;
                    plint mask = util::roundToInt( lattices[subsNum]->get(iXm,iYm,iZm).computeDensity() );
                    // Only process pore and biofilm voxels (skip boundaries and solids)
                    if (mask != bb && mask != solid) {
                        // Determine if this voxel is currently classified as pore
                        bool poreflag = 0;
                        for (size_t iP=0; iP<pore.size(); ++iP) { if (mask == pore[iP]) {poreflag = 1; break;} }
                        // Check and update omega for each substrate species
                        for (plint iS=0; iS<subsNum; ++iS) {
                            plint iXs = iX + vec_offset[iS].x, iYs = iY + vec_offset[iS].y, iZs = iZ + vec_offset[iS].z;
                            T omega = lattices[iS]->get(iXs,iYs,iZs).getDynamics().getOmega();
                            T bMassOmega = substrOMEGAinbMass[iS];
                            // Voxel is biofilm but omega still set to pore value → switch to biofilm omega
                            if (poreflag == 0 && std::abs(omega-bMassOmega) > COMPLAB_THRD ) {
                                lattices[iS]->get(iXs,iYs,iZs).getDynamics().setOmega(bMassOmega); // pore to bfilm
                            }
                            // Voxel is pore but omega still set to biofilm value → switch to pore omega
                            else if (poreflag == 1 && std::abs(omega-bMassOmega) < COMPLAB_THRD ) {
                                lattices[iS]->get(iXs,iYs,iZs).getDynamics().setOmega(substrOMEGAinPore[iS]); // bfilm to pore
                            }
                        }
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;  // must include envelope for correct MPI communication
    }
    virtual updateSoluteDynamics3D<T,Descriptor>* clone() const {
        return new updateSoluteDynamics3D<T,Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iS=0; iS<subsNum; ++iS) {
            modified[iS] = modif::dynamicVariables;  // omega (dynamics parameter) may change
        }
        modified[subsNum] = modif::nothing;  // mask lattice is read-only
    }
private:
    plint subsNum, bb, solid;
    std::vector<plint> pore;
    std::vector<T> substrOMEGAinbMass, substrOMEGAinPore;  // omega values for biofilm vs. pore
};

/* -----------------------------------------------------------------------------------------------------------------
 * updateBiomassDynamics3D — Planktonic Biomass LBM Relaxation Rate Update
 *
 * Analogous to updateSoluteDynamics3D above, but for planktonic (suspended) biomass
 * lattices rather than dissolved substrates. When a voxel switches between pore and
 * biofilm, the effective diffusivity of planktonic cells changes:
 *
 *   - Pore → Biofilm:  omega switches to bmassOMEGAinbMass (slower diffusion in biofilm)
 *   - Biofilm → Pore:  omega switches to bmassOMEGAinPore (free diffusion in pore)
 *
 * This ensures planktonic cells move at the correct rate through biofilm vs. open pore.
 * ----------------------------------------------------------------------------------------------------------------- */
template<typename T, template<typename U> class Descriptor>
class updateBiomassDynamics3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    /* Constructor parameters:
     *   bioNum_            — number of planktonic biomass species
     *   bb_, solid_        — mask values for boundary and solid
     *   pore_              — list of pore-type mask values
     *   bmassOMEGAinbMass_ — omega for each biomass species inside biofilm
     *   bmassOMEGAinPore_  — omega for each biomass species in open pore   */
    updateBiomassDynamics3D(plint bioNum_, plint bb_, plint solid_, std::vector<plint> pore_, std::vector<T> bmassOMEGAinbMass_, std::vector<T> bmassOMEGAinPore_)
    : bioNum(bioNum_), bb(bb_), solid(solid_), pore(pore_), bmassOMEGAinbMass(bmassOMEGAinbMass_), bmassOMEGAinPore(bmassOMEGAinPore_)
    {}
    /* Lattice layout:
     *   lattices[0 .. bioNum-1]  = planktonic biomass LBM lattices
     *   lattices[bioNum]         = mask lattice (read-only)          */
    // lattices[0~(#ofbMs-1)] = planktonic biomass lattices
    // lattices[#ofbMs] = mask lattice
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        std::vector<Dot3D> vec_offset;
        for (plint iT=0; iT<bioNum+1; ++iT) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0],*lattices[iT])); }
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint iXm = iX+vec_offset[bioNum].x;
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint iYm = iY+vec_offset[bioNum].y;
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint iZm = iZ+vec_offset[bioNum].z;
                    plint mask = util::roundToInt( lattices[bioNum]->get(iXm,iYm,iZm).computeDensity() );
                    if (mask != bb && mask != solid) {
                        bool poreflag = 0;
                        for (size_t iP=0; iP<pore.size(); ++iP) { if (mask == pore[iP]) {poreflag = 1; break;} }
                        // Check and update omega for each planktonic biomass species
                        for (plint iB=0; iB<bioNum; ++iB) {
                            plint iXb = iX + vec_offset[iB].x, iYb = iY + vec_offset[iB].y, iZb = iZ + vec_offset[iB].z;
                            T omega = lattices[iB]->get(iXb,iYb,iZb).getDynamics().getOmega();
                            T bMassOmega = bmassOMEGAinbMass[iB];
                            // Voxel is biofilm but omega still at pore value → switch to biofilm omega
                            if (poreflag == 0 && std::abs(omega-bMassOmega) > COMPLAB_THRD ) {
                                lattices[iB]->get(iXb,iYb,iZb).getDynamics().setOmega(bMassOmega); // pore to bfilm
                            }
                            // Voxel is pore but omega still at biofilm value → switch to pore omega
                            else if (poreflag == 1 && std::abs(omega-bMassOmega) < COMPLAB_THRD ) {
                                lattices[iB]->get(iXb,iYb,iZb).getDynamics().setOmega(bmassOMEGAinPore[iB]); // bfilm to pore
                            }
                        }
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    virtual updateBiomassDynamics3D<T,Descriptor>* clone() const {
        return new updateBiomassDynamics3D<T,Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iB=0; iB<bioNum; ++iB) {
            modified[iB] = modif::dynamicVariables;  // omega may change
        }
        modified[bioNum] = modif::nothing;  // mask lattice is read-only
    }
private:
    plint bioNum, bb, solid;
    std::vector<plint> pore;
    std::vector<T> bmassOMEGAinbMass, bmassOMEGAinPore;  // omega values for biofilm vs. pore
};

/* -----------------------------------------------------------------------------------------------------------------
 * updateNsLatticesDynamics3D — Navier-Stokes Flow Field Dynamics Update
 *
 * When voxel reclassification (Algorithm B.5) changes the pore/biofilm geometry, the
 * flow field must be updated accordingly. This processor switches the LBM dynamics
 * at each voxel in the Navier-Stokes lattice:
 *
 *   - Pore → Biofilm:
 *       • If bioX ≤ 0 (impermeable biofilm): assign BounceBack dynamics (no-flow wall)
 *       • If bioX > 0 (permeable biofilm):   assign IncBGKdynamics with bioOmega
 *         (reduced permeability via Brinkman-penalized viscosity)
 *   - Biofilm → Pore:
 *       • Restore IncBGKdynamics with the original pore omega (nsOmega) for free flow
 *
 * bioOmega is derived from nsOmega and the permeability ratio bioX using:
 *     bioOmega = 1 / (bioX × (1/nsOmega − 0.5) + 0.5)
 *
 * After this processor runs, the steady-state flow field is recomputed (Section 3.1).
 * ----------------------------------------------------------------------------------------------------------------- */
template<typename T1, template<typename U> class Descriptor1, typename T2, template<typename U> class Descriptor2>
class updateNsLatticesDynamics3D : public BoxProcessingFunctional3D_LL<T1,Descriptor1,T2,Descriptor2>
{
public:
    /* Constructor parameters:
     *   nsOmega_  — LBM relaxation rate for free pore flow (from fluid viscosity)
     *   bioX_     — biofilm permeability ratio (0 = impermeable → BounceBack; >0 = permeable)
     *   pore_     — list of pore-type mask values
     *   solid_    — mask value for solid grains
     *   bb_       — mask value for domain boundaries                                       */
    updateNsLatticesDynamics3D(T nsOmega_, T bioX_, std::vector<plint> pore_, plint solid_, plint bb_)
    : nsOmega(nsOmega_), bioX(bioX_), pore(pore_), solid(solid_), bb(bb_)
    {}
    /* Two-lattice processor (Navier-Stokes + mask):
     *   lattice0 = Navier-Stokes flow field lattice (dynamics will be changed)
     *   lattice1 = mask field lattice (read-only, provides voxel identity)      */
    // lattice0 = flow field lattice
    // lattice1 = mask field lattice
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor1>& lattice0, BlockLattice3D<T2,Descriptor2>& lattice1) {
        // Precompute the biofilm omega from nsOmega and permeability ratio bioX
        // Formula: bioOmega = 1 / (bioX × (1/nsOmega − 0.5) + 0.5)
        T bioOmega = 1/(bioX*(1/nsOmega-.5)+.5);
        Dot3D offset_12 = computeRelativeDisplacement(lattice0,lattice1);
        for (plint iX0=domain.x0; iX0<=domain.x1; ++iX0) {
            plint iX1 = iX0 + offset_12.x;
            for (plint iY0=domain.y0; iY0<=domain.y1; ++iY0) {
                plint iY1 = iY0 + offset_12.y;
                for (plint iZ0=domain.z0; iZ0<=domain.z1; ++iZ0) {
                    plint iZ1 = iZ0 + offset_12.z;
                    plint mask = util::roundToInt( lattice1.get(iX1,iY1,iZ1).computeDensity() );
                    T currentOmega = lattice0.get(iX0,iY0,iZ0).getDynamics().getOmega();
                    if (mask != bb || mask != solid) {
                        bool poreflag = 0;
                        for (size_t iP=0; iP<pore.size(); ++iP) { if (mask == pore[iP]) {poreflag = 1; break;} }
                        // Pore → Biofilm: voxel is now biofilm but still has pore-flow dynamics
                        if ( poreflag == 0 && std::abs(currentOmega-nsOmega) < COMPLAB_THRD) {
                            // bioX ≤ 0: impermeable biofilm → assign BounceBack (acts as solid wall)
                            if (bioX <= COMPLAB_THRD ) { lattice0.attributeDynamics(iX0, iY0, iZ0, new BounceBack<T1,Descriptor1>() ); }
                            // bioX > 0: permeable biofilm → assign IncBGKdynamics with reduced permeability
                            else { lattice0.attributeDynamics(iX0, iY0, iZ0, new IncBGKdynamics<T1,Descriptor1>(bioOmega) ); }
                        }
                        // Biofilm → Pore: voxel is now pore but still has biofilm dynamics
                        else if ( poreflag == 1 && std::abs(currentOmega-nsOmega) > COMPLAB_THRD) {
                            // Restore free-flow dynamics with original pore omega
                            lattice0.attributeDynamics(iX0, iY0, iZ0, new IncBGKdynamics<T1,Descriptor1>(nsOmega) );
                        }
                    }
                }
            }
        }
    }
    virtual updateNsLatticesDynamics3D<T1,Descriptor1,T2,Descriptor2>* clone() const {
        return new updateNsLatticesDynamics3D<T1,Descriptor1,T2,Descriptor2>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::dataStructure;  // dynamics object itself may be replaced (BGK ↔ BounceBack)
        modified[1] = modif::nothing;        // mask lattice is read-only
    }
private :
    T nsOmega, bioX;           // nsOmega = pore relaxation rate; bioX = biofilm permeability ratio
    std::vector<plint> pore;   // list of pore-type mask values
    plint solid, bb;           // mask codes for solid grains and domain boundaries
};

/* -----------------------------------------------------------------------------------------------------------------
 * updateAgeDistance3D — Age/Distance Lattice Update for Redistribution Stage 2
 *
 * The "age" lattice stores a wavefront counter that propagates outward from the
 * pore-biofilm boundary into the biofilm interior. It is used by Stage 2 of the
 * redistribution algorithm (Algorithms B.2/B.3) to determine the distance-to-nearest-pore
 * for each biofilm voxel. Voxels with a LOWER age value are closer to open pore space,
 * so Stage 2 pushes excess biomass toward neighbors with smaller age (i.e., toward the
 * nearest pore boundary).
 *
 * The age propagation works as follows:
 *   - age = 0: voxel has biomass but has not yet been initialized → set to 1 (pore boundary)
 *   - age = 1 and voxel is at or above Bmax: if ALL valid neighbors already have age ≥ 1,
 *     advance to age = 2 (one layer deeper into biofilm)
 *   - age ≥ 2: advance to age+1 only if ALL valid neighbors have age ≥ current age AND
 *     biomass ≥ Bmax (they are all fully interior biofilm)
 *
 * The "dist" lattice is a companion that marks which voxels participate in the wavefront
 * (dist > 0 means the voxel is part of the biofilm region where age propagation applies).
 * ----------------------------------------------------------------------------------------------------------------- */
template<typename T, template<typename U> class Descriptor>
class updateAgeDistance3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    /* Constructor parameters:
     *   Bmax_          — max_bMassRho (maximum biomass density, redistribution trigger)
     *   nx_, ny_, nz_  — global domain dimensions (for boundary checks)              */
    updateAgeDistance3D(T Bmax_, plint nx_, plint ny_, plint nz_)
    : Bmax(Bmax_), nx(nx_), ny(ny_), nz(nz_)
    {}
    /* Lattice layout:
     *   lattices[0] = age lattice    (wavefront counter, modified)
     *   lattices[1] = dist lattice   (distance field, read-only — marks biofilm region)
     *   lattices[2] = total biomass lattice (read-only)                                */
    // lattices[0] = age 0 lattices
    // lattices[1] = dist lattices
    // lattices[2] = total biomass lattice
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        std::vector<Dot3D> vec_offset;
        Dot3D absoluteOffset = lattices[0]->getLocation();
        plint ageLloc = 0, distLloc = 1, bMtLloc = 2;  // lattice index aliases
        for (plint iL=0; iL<3; ++iL) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0],*lattices[iL])); }
        for (plint iX0=domain.x0; iX0<=domain.x1; ++iX0) {
            plint iXb = iX0+vec_offset[bMtLloc].x;
            for (plint iY0=domain.y0; iY0<=domain.y1; ++iY0) {
                plint iYb = iY0+vec_offset[bMtLloc].y;
                for (plint iZ0=domain.z0; iZ0<=domain.z1; ++iZ0) {
                    plint iZb = iZ0+vec_offset[bMtLloc].z;
                    // Read total biomass at this voxel
                    T B = lattices[bMtLloc]->get(iXb,iYb,iZb).computeDensity();
                    // Only process voxels that contain biomass (skip empty pore/solid)
                    if (B > COMPLAB_THRD ) {
                        plint absX = iX0 + absoluteOffset.x, absY = iY0 + absoluteOffset.y, absZ = iZ0 + absoluteOffset.z;
                        plint iXa = iX0 + vec_offset[ageLloc].x, iYa = iY0 + vec_offset[ageLloc].y, iZa = iZ0 + vec_offset[ageLloc].z;
                        plint iXd = iX0 + vec_offset[distLloc].x, iYd = iY0 + vec_offset[distLloc].y, iZd = iZ0 + vec_offset[distLloc].z;
                        std::vector<plint> delXYZ; plint nbrs=0;

                        // Build the list of face-connected neighbors (D3Q7: ±x, ±y, ±z)
                        // excluding neighbors that fall outside the domain
                        bool atXmin = (absX == 0);
                        bool atXmax = (absX == nx-1);
                        bool atYmin = (absY == 0);
                        bool atYmax = (absY == ny-1);
                        bool atZmin = (absZ == 0);
                        bool atZmax = (absZ == nz-1);

                        // Add valid neighbors (stored as consecutive dx,dy,dz triplets)
                        if (!atXmax) { delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); ++nbrs; }
                        if (!atXmin) { delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); ++nbrs; }
                        if (!atYmax) { delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); ++nbrs; }
                        if (!atYmin) { delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0); ++nbrs; }
                        if (!atZmax) { delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); ++nbrs; }
                        if (!atZmin) { delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1); ++nbrs; }

                        // Filter to only keep neighbors that are in the biofilm region (dist > 0)
                        std::vector<plint> del_iX, del_iY, del_iZ;
                        for (plint iT=0; iT<nbrs; ++iT) {
                            plint delx = delXYZ[iT*3], dely = delXYZ[iT*3+1], delz = delXYZ[iT*3+2];
                            plint dist = util::roundToInt( lattices[distLloc]->get(iXd+delx,iYd+dely,iZd+delz).computeDensity() );
                            if (dist > 0) { del_iX.push_back(delx); del_iY.push_back(dely); del_iZ.push_back(delz); }
                        }
                        T newAge = 0; bool updateflag = 1;
                        // Read the current age (wavefront counter) at this voxel
                        plint age = util::roundToInt( lattices[ageLloc]->get(iXa,iYa,iZa).computeDensity() );

                        // --- Age propagation rules ---
                        // Age 0 → 1: initial seeding — this voxel just entered the wavefront
                        if (age==0 ) { newAge = 1; }
                        else if ( age==1 && (B-Bmax)>-COMPLAB_THRD) {
                            // Age 1 and biomass ≈ Bmax: ready to advance to age 2
                            // But only if NO valid neighbor still has age 0 (all neighbors initialized)
                            for (size_t iT=0; iT<del_iX.size(); ++iT) {
                                plint delx = del_iX[iT], dely = del_iY[iT], delz = del_iZ[iT];
                                plint nbrAge = util::roundToInt( lattices[ageLloc]->get(iXa+delx,iYa+dely,iZa+delz).computeDensity() );
                                plint nbrDist = util::roundToInt( lattices[distLloc]->get(iXd+delx,iYd+dely,iZd+delz).computeDensity() );
                                // If any biofilm neighbor (dist>0) still has age 0, we can't advance yet
                                if ( nbrDist > 0 && nbrAge == 0 ) { updateflag = 0; break; }
                            }
                            if (updateflag == 1) { newAge = 2; }
                        }
                        else {
                            // Age ≥ 2: advance to age+1 only if ALL valid neighbors have
                            // age ≥ current age AND biomass ≥ Bmax (fully interior biofilm)
                            // This ensures the wavefront propagates layer by layer into the interior
                            plint count = 0;
                            for (size_t iT=0; iT<del_iX.size(); ++iT) {
                                plint delx = del_iX[iT], dely = del_iY[iT], delz = del_iZ[iT];
                                plint nbrAge = util::roundToInt( lattices[ageLloc]->get(iXa+delx,iYa+dely,iZa+delz).computeDensity() );
                                T nbrB = lattices[bMtLloc]->get(iXb+delx,iYb+dely,iZb+delz).computeDensity();
                                // Count neighbors that are at least as old and at capacity
                                if ( nbrAge >= age && (nbrB-Bmax)>-COMPLAB_THRD) { ++count; }
                            }
                            // Only advance if ALL valid neighbors qualify
                            if (count == del_iX.size()) { newAge = (T) (age+1); }
                            else { updateflag = 0; }
                        }
                        // Write the updated age value into the age lattice (D3Q7 density encoding)
                        if (updateflag == 1) {
                            Array<T,7> g;
                            g[0]=(T) (newAge-1)/4; g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T) (newAge-1)/8;
                            lattices[ageLloc]->get(iXa,iYa,iZa).setPopulations(g);
                        }
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulk;
    }
    virtual updateAgeDistance3D<T,Descriptor>* clone() const {
        return new updateAgeDistance3D<T,Descriptor>(*this);
    }
    /* Modification flags:
     *   age lattice [0]            — modified (wavefront counter updated)
     *   dist lattice [1]           — not modified (read-only reference)
     *   total biomass lattice [2]  — not modified (read-only)            */
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
        modified[2] = modif::nothing;
    }
private:
    T Bmax;                // max_bMassRho — maximum biomass density (redistribution trigger)
    plint nx, ny, nz;      // global domain dimensions
};

#endif // COMPLAB3D_PROCESSORS_PART2_HH
