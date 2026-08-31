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

#ifndef COMPLAB3D_PROCESSORS_PART1_HH
#define COMPLAB3D_PROCESSORS_PART1_HH

#include "../defineKinetics.hh"
#include "../defineAbioticKinetics.hh"  // For abiotic (substrate-only) reactions
#include <random>
using namespace plb;
typedef double T;

// [UPD_RXN] debug-print toggle. Set from CompLaB.xml <parameters><IO><debug_updRxn> at startup.
//   0 = OFF (clean output, default), 1 = ON. Global so the reaction data processor can read it.
bool g_updRxn_debug = false;

#define NSDES descriptors::D3Q19Descriptor // Cs2 = 1/3
#define RXNDES descriptors::AdvectionDiffusionD3Q7Descriptor // Cs2 = 1/3

/* ===============================================================================================================
   =============================================== DATA PROCESSORS ===============================================
   =============================================================================================================== */

template<typename T, template<typename U> class Descriptor>
class run_kinetics : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    run_kinetics( plint nx_, plint subsNum_, plint bioNum_, T dt_ , std::vector<std::vector<T>> vec2_Kc_kns_, std::vector<T> vec1_mu_kns_, plint solid_, plint bb_)
    : nx(nx_), subsNum(subsNum_), bioNum(bioNum_), dt(dt_), vec2_Kc_kns(vec2_Kc_kns_), vec1_mu_kns(vec1_mu_kns_), solid(solid_), bb(bb_), dCloc(subsNum_+bioNum_), dBloc(subsNum_*2+bioNum_), maskLloc(2*(subsNum_+bioNum_))
    {}
    // substrate lattices and then bio-lattices. the mask number lattice is always the last.
    // dt in seconds, dx in meters
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        Dot3D absoluteOffset = lattices[0]->getLocation();
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint absX = iX+absoluteOffset.x;
            if (absX > 0 && absX < nx-1) {
                for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                    for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                        Dot3D maskOffset = computeRelativeDisplacement(*lattices[0],*lattices[maskLloc]);
                        plint mask = util::roundToInt( lattices[maskLloc]->get(iX+maskOffset.x,iY+maskOffset.y,iZ+maskOffset.z).computeDensity() );
                        if ( mask != solid && mask != bb) {
                            std::vector<Dot3D> vec_offset;
                            for (plint iT=0; iT<maskLloc; ++iT) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0],*lattices[iT])); }
                            std::vector<T> bmass, conc, subs_rate(subsNum,0), bio_rate(bioNum,0);
                            // construct the concentration vector
                            for (plint iS=0; iS<subsNum; ++iS) {
                                plint iXs = iX+vec_offset[iS].x, iYs = iY+vec_offset[iS].y, iZs = iZ+vec_offset[iS].z;
                                T c0 = lattices[iS]->get(iXs,iYs,iZs).computeDensity(); // mmol/L
                                if (c0 < thrd) { c0 = 0; }
                                conc.push_back(c0); // mmol/L
                            }
                            // construct the biomass vector
                            for (plint iM=0; iM<bioNum; ++iM) {
                                plint iXb = iX+vec_offset[subsNum+iM].x, iYb = iY+vec_offset[subsNum+iM].y, iZb = iZ+vec_offset[subsNum+iM].z;
                                T b0 = lattices[subsNum+iM]->get(iXb,iYb,iZb).computeDensity(); // kgDW/m3 = gDW/L
                                if (b0 < thrd) { b0 = 0; }
                                bmass.push_back(b0); // vector storing biomass (unit: kgDW/m3)
                            }

                            defineRxnKinetics( bmass, conc, subs_rate, bio_rate, mask );

                            // update dC
                            for (plint iS=0; iS<subsNum; ++iS) {
                                T dC = subs_rate[iS]*dt;
                                if (dC > thrd || dC < -thrd) {
                                    Array<T,7> g;
                                    lattices[iS+dCloc]->get(iX+vec_offset[iS+dCloc].x,iY+vec_offset[iS+dCloc].y,iZ+vec_offset[iS+dCloc].z).getPopulations(g);
                                    g[0]+=(T) (dC)/4; g[1]+=(T) (dC)/8; g[2]+=(T) (dC)/8; g[3]+=(T) (dC)/8; g[4]+=(T) (dC)/8; g[5]+=(T) (dC)/8; g[6]+=(T) (dC)/8;
                                    lattices[iS+dCloc]->get(iX+vec_offset[iS+dCloc].x,iY+vec_offset[iS+dCloc].y,iZ+vec_offset[iS+dCloc].z).setPopulations(g);
                                }
                            }
                            // update dB
                            for (plint iB=0; iB<bioNum; ++iB) {
                                T dB = bio_rate[iB]*dt; // unit: kgDW/m3 (=gDW/L)
                                if (dB > thrd || dB < -thrd) {
                                    Array<T,7> g;
                                    lattices[dBloc+iB]->get(iX+vec_offset[dBloc+iB].x,iY+vec_offset[dBloc+iB].y,iZ+vec_offset[dBloc+iB].z).getPopulations(g);
                                    g[0]+=(T) (dB)/4; g[1]+=(T) (dB)/8; g[2]+=(T) (dB)/8; g[3]+=(T) (dB)/8; g[4]+=(T) (dB)/8; g[5]+=(T) (dB)/8; g[6]+=(T) (dB)/8;
                                    lattices[dBloc+iB]->get(iX+vec_offset[dBloc+iB].x,iY+vec_offset[dBloc+iB].y,iZ+vec_offset[dBloc+iB].z).setPopulations(g);
                                }
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
    virtual run_kinetics<T,Descriptor>* clone() const {
        return new run_kinetics<T,Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iT = dCloc; iT < maskLloc; ++iT) {
            modified[iT] = modif::staticVariables;
        }
    }
private:
    plint nx, subsNum, bioNum;
    T dt;
    std::vector<std::vector<T>> vec2_Kc_kns;
    std::vector<T> vec1_mu_kns;
    plint solid, bb;
    plint dCloc, dBloc, maskLloc;
};

/* ============================================================================
 * ABIOTIC KINETICS PROCESSOR
 * ============================================================================
 * Processes chemical reactions between substrates WITHOUT biomass.
 * Used when biotic_mode=false but enable_abiotic_kinetics=true.
 *
 * Lattice order: [C0, C1, ..., Cn, dC0, dC1, ..., dCn, mask]
 * Total lattices: 2*subsNum + 1
 * ============================================================================
 */
template<typename T, template<typename U> class Descriptor>
class run_abiotic_kinetics : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    run_abiotic_kinetics(plint nx_, plint subsNum_, T dt_, plint solid_, plint bb_)
    : nx(nx_), subsNum(subsNum_), dt(dt_), solid(solid_), bb(bb_),
      dCloc(subsNum_), maskLloc(2*subsNum_)
    {}

    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        Dot3D absoluteOffset = lattices[0]->getLocation();

        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;
            if (absX > 0 && absX < nx-1) {
                for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                    for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                        // Get mask value
                        Dot3D maskOffset = computeRelativeDisplacement(*lattices[0], *lattices[maskLloc]);
                        plint mask = util::roundToInt(lattices[maskLloc]->get(
                            iX+maskOffset.x, iY+maskOffset.y, iZ+maskOffset.z).computeDensity());

                        // Skip solid and bounce-back cells
                        if (mask != solid && mask != bb) {
                            // Get offsets for all lattices
                            std::vector<Dot3D> vec_offset;
                            for (plint iT=0; iT<maskLloc; ++iT) {
                                vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT]));
                            }

                            // Construct concentration vector
                            std::vector<T> conc;
                            for (plint iS=0; iS<subsNum; ++iS) {
                                plint iXs = iX + vec_offset[iS].x;
                                plint iYs = iY + vec_offset[iS].y;
                                plint iZs = iZ + vec_offset[iS].z;
                                T c0 = lattices[iS]->get(iXs, iYs, iZs).computeDensity();
                                if (c0 < thrd) c0 = 0;
                                conc.push_back(c0);
                            }

                            // Calculate abiotic reaction rates
                            std::vector<T> subs_rate(subsNum, 0.0);
                            defineAbioticRxnKinetics(conc, subs_rate, mask);

                            // Update dC lattices
                            for (plint iS=0; iS<subsNum; ++iS) {
                                T dC = subs_rate[iS] * dt;
                                if (dC > thrd || dC < -thrd) {
                                    plint iXd = iX + vec_offset[iS+dCloc].x;
                                    plint iYd = iY + vec_offset[iS+dCloc].y;
                                    plint iZd = iZ + vec_offset[iS+dCloc].z;
                                    Array<T,7> g;
                                    lattices[iS+dCloc]->get(iXd, iYd, iZd).getPopulations(g);
                                    g[0] += (T)(dC)/4;
                                    g[1] += (T)(dC)/8;
                                    g[2] += (T)(dC)/8;
                                    g[3] += (T)(dC)/8;
                                    g[4] += (T)(dC)/8;
                                    g[5] += (T)(dC)/8;
                                    g[6] += (T)(dC)/8;
                                    lattices[iS+dCloc]->get(iXd, iYd, iZd).setPopulations(g);
                                }
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

    virtual run_abiotic_kinetics<T,Descriptor>* clone() const {
        return new run_abiotic_kinetics<T,Descriptor>(*this);
    }

    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        // dC lattices are modified
        for (plint iT = dCloc; iT < maskLloc; ++iT) {
            modified[iT] = modif::staticVariables;
        }
    }

private:
    plint nx, subsNum;
    T dt;
    plint solid, bb;
    plint dCloc, maskLloc;
};

/* ============================================================================
 * ABIOTIC KINETICS UPDATE PROCESSOR
 * ============================================================================
 * Applies the accumulated dC changes to substrate lattices.
 * Called after run_abiotic_kinetics to update concentrations.
 * ============================================================================
 */
template<typename T, template<typename U> class Descriptor>
class update_abiotic_rxnLattices : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    update_abiotic_rxnLattices(plint nx_, plint subsNum_, plint solid_, plint bb_)
    : nx(nx_), subsNum(subsNum_), solid(solid_), bb(bb_),
      dCloc(subsNum_), maskLloc(2*subsNum_)
    {}

    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        Dot3D absoluteOffset = lattices[0]->getLocation();

        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;
            if (absX > 0 && absX < nx-1) {
                for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                    for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                        // Get mask value
                        Dot3D maskOffset = computeRelativeDisplacement(*lattices[0], *lattices[maskLloc]);
                        plint mask = util::roundToInt(lattices[maskLloc]->get(
                            iX+maskOffset.x, iY+maskOffset.y, iZ+maskOffset.z).computeDensity());

                        // Skip solid and bounce-back cells
                        if (mask != solid && mask != bb) {
                            std::vector<Dot3D> vec_offset;
                            for (plint iT=0; iT<maskLloc; ++iT) {
                                vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT]));
                            }

                            // Apply dC to substrate concentrations
                            for (plint iS=0; iS<subsNum; ++iS) {
                                plint iXd = iX + vec_offset[iS+dCloc].x;
                                plint iYd = iY + vec_offset[iS+dCloc].y;
                                plint iZd = iZ + vec_offset[iS+dCloc].z;
                                T dC = lattices[iS+dCloc]->get(iXd, iYd, iZd).computeDensity();

                                if (dC > thrd || dC < -thrd) {
                                    plint iXs = iX + vec_offset[iS].x;
                                    plint iYs = iY + vec_offset[iS].y;
                                    plint iZs = iZ + vec_offset[iS].z;
                                    Array<T,7> g;
                                    lattices[iS]->get(iXs, iYs, iZs).getPopulations(g);
                                    g[0] += (T)(dC)/4;
                                    g[1] += (T)(dC)/8;
                                    g[2] += (T)(dC)/8;
                                    g[3] += (T)(dC)/8;
                                    g[4] += (T)(dC)/8;
                                    g[5] += (T)(dC)/8;
                                    g[6] += (T)(dC)/8;
                                    lattices[iS]->get(iXs, iYs, iZs).setPopulations(g);
                                }
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

    virtual update_abiotic_rxnLattices<T,Descriptor>* clone() const {
        return new update_abiotic_rxnLattices<T,Descriptor>(*this);
    }

    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        // Substrate lattices are modified
        for (plint iT = 0; iT < subsNum; ++iT) {
            modified[iT] = modif::staticVariables;
        }
    }

private:
    plint nx, subsNum;
    plint solid, bb;
    plint dCloc, maskLloc;
};

template<typename T, template<typename U> class Descriptor>
class update_rxnLattices : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    update_rxnLattices( plint nx_, plint subsNum_, plint bioNum_, plint solid_, plint bb_ )
    : nx(nx_), subsNum(subsNum_), bioNum(bioNum_), solid(solid_), bb(bb_), dCloc(subsNum_+bioNum_), dBloc(subsNum_*2+bioNum_), maskLloc(2*(subsNum_+bioNum_))
    {}
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        Dot3D absoluteOffset = lattices[0]->getLocation();
        
        // DEBUG counters
        static long long total_calls = 0;
        static long long dB_applied = 0;
        static double max_dB_seen = 0;
        static double sum_dB_applied = 0;
        static bool header_printed = false;
        
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint absX = iX+absoluteOffset.x;
            if (absX > 0 && absX < nx-1) {
                for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                    for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                        Dot3D maskOffset = computeRelativeDisplacement(*lattices[0],*lattices[maskLloc]);
                        plint mask = util::roundToInt( lattices[maskLloc]->get(iX+maskOffset.x,iY+maskOffset.y,iZ+maskOffset.z).computeDensity() );
                        if ( mask != solid && mask != bb) {
                            std::vector<Dot3D> vec_offset;
                            for (plint iT=0; iT<maskLloc; ++iT) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0],*lattices[iT])); }

                            // update concentration
                            for (plint iS=0; iS<subsNum; ++iS) {
                                plint iXs = iX+vec_offset[iS+dCloc].x, iYs = iY+vec_offset[iS+dCloc].y, iZs = iZ+vec_offset[iS+dCloc].z;
                                T dC = lattices[iS+dCloc]->get(iXs,iYs,iZs).computeDensity();
                                if (dC > thrd || dC < -thrd) {
                                    Array<T,7> g;
                                    lattices[iS]->get(iX+vec_offset[iS].x,iY+vec_offset[iS].y,iZ+vec_offset[iS].z).getPopulations(g);
                                    g[0]+=(T) (dC)/4; g[1]+=(T) (dC)/8; g[2]+=(T) (dC)/8; g[3]+=(T) (dC)/8; g[4]+=(T) (dC)/8; g[5]+=(T) (dC)/8; g[6]+=(T) (dC)/8;
                                    lattices[iS]->get(iX+vec_offset[iS].x,iY+vec_offset[iS].y,iZ+vec_offset[iS].z).setPopulations(g);
                                }
                            }

                            // update biomass
                            for (plint iB=0; iB<bioNum; ++iB) {
                                plint iXb = iX+vec_offset[iB+dBloc].x, iYb = iY+vec_offset[iB+dBloc].y, iZb = iZ+vec_offset[iB+dBloc].z;
                                T dB = lattices[iB+dBloc]->get(iXb,iYb,iZb).computeDensity();
                                
                                total_calls++;
                                if (std::abs(dB) > max_dB_seen) max_dB_seen = std::abs(dB);
                                
                                if (dB > thrd || dB < -thrd) {
                                    dB_applied++;
                                    sum_dB_applied += dB;
                                    
                                    Array<T,7> g;
                                    lattices[subsNum+iB]->get(iX+vec_offset[subsNum+iB].x,iY+vec_offset[subsNum+iB].y,iZ+vec_offset[subsNum+iB].z).getPopulations(g);
                                    g[0]+=(T) (dB)/4; g[1]+=(T) (dB)/8; g[2]+=(T) (dB)/8; g[3]+=(T) (dB)/8; g[4]+=(T) (dB)/8; g[5]+=(T) (dB)/8; g[6]+=(T) (dB)/8;
                                    lattices[subsNum+iB]->get(iX+vec_offset[subsNum+iB].x,iY+vec_offset[subsNum+iB].y,iZ+vec_offset[subsNum+iB].z).setPopulations(g);
                                }
                            }

                        }
                    }
                }
            }
        }
        
        // DEBUG output every 10M calls -- OFF by default.
        // Toggle from CompLaB.xml:  <IO><debug_updRxn> 1 </debug_updRxn>   (0 = off)
        if (g_updRxn_debug && (total_calls % 10000000 < 50000)) {
            if (!header_printed) {
                std::cout << "\n[UPD_RXN DEBUG] ============================================\n";
                header_printed = true;
            }
            std::cout << "[UPD_RXN] total=" << total_calls 
                      << " dB_applied=" << dB_applied 
                      << " max_dB=" << max_dB_seen 
                      << " sum_dB=" << sum_dB_applied << "\n";
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    virtual update_rxnLattices<T,Descriptor>* clone() const {
        return new update_rxnLattices<T,Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iT = 0; iT < subsNum+bioNum; ++iT) {
            modified[iT] = modif::staticVariables;
        }
    }
private:
    plint nx, subsNum, bioNum, solid, bb, dCloc, dBloc, maskLloc;
};

// redistribute the excessive biomass (push)
template<typename T, template<typename U> class Descriptor>
class pushExcessBiomass3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    pushExcessBiomass3D(T Bmax_, plint nx_, plint ny_, plint nz_, plint bdryGap_, plint length_, plint solid_, plint bb_, std::vector<plint> pore_)
    : Bmax(Bmax_), nx(nx_), ny(ny_), nz(nz_), bdryGap(bdryGap_), length(length_), solid(solid_), bb(bb_), pore(pore_)
    {}
    // lattices[0~(#ofbM-1)] = original biomass lattices
    // lattices[#ofbM~(len-3)] = copy biomass lattices
    // lattices[len-3] = total biomass lattice
    // lattices[len-2] = mask lattice
    // lattices[len-1] = age lattice
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        std::vector<Dot3D> vec_offset;
        plint distLloc = length-1, maskLloc = length-2, bMtLloc = length-3, numbM = (length-3)/2;
        for (plint iL=0; iL<length; ++iL) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0],*lattices[iL])); }
        Dot3D absoluteOffset = lattices[0]->getLocation();
        for (plint iX0=domain.x0; iX0<=domain.x1; ++iX0) {
            plint iXm = iX0+vec_offset[maskLloc].x;
            for (plint iY0=domain.y0; iY0<=domain.y1; ++iY0) {
                plint iYm = iY0+vec_offset[maskLloc].y;
                for (plint iZ0=domain.z0; iZ0<=domain.z1; ++iZ0) {
                    plint iZm = iZ0+vec_offset[maskLloc].z;
                    plint mask = util::roundToInt( lattices[maskLloc]->get(iXm,iYm,iZm).computeDensity() );
                    if (mask != bb && mask != solid) {
                        plint iXt = iX0 + vec_offset[bMtLloc].x, iYt = iY0 + vec_offset[bMtLloc].y, iZt = iZ0 + vec_offset[bMtLloc].z;
                        T bMt = lattices[bMtLloc]->get(iXt,iYt,iZt).computeDensity();
                        if (bMt > Bmax) {
                            T bMd = bMt-Bmax;
                            if ( bMd > thrd ) {
                                plint absX = iX0 + absoluteOffset.x, absY = iY0 + absoluteOffset.y, absZ = iZ0 + absoluteOffset.z;
                                std::vector<plint> delXYZ; plint nbrs = 0;
                                
                                // Determine neighbors based on boundary position (3D version)
                                // Interior point: 6 neighbors
                                // Face: 5 neighbors, Edge: 4 neighbors, Corner: 3 neighbors
                                bool atXmin = (absX == bdryGap);
                                bool atXmax = (absX == nx-1-bdryGap);
                                bool atYmin = (absY == 0);
                                bool atYmax = (absY == ny-1);
                                bool atZmin = (absZ == 0);
                                bool atZmax = (absZ == nz-1);
                                
                                // Add valid neighbors
                                if (!atXmax) { delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); ++nbrs; }
                                if (!atXmin) { delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); ++nbrs; }
                                if (!atYmax) { delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); ++nbrs; }
                                if (!atYmin) { delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0); ++nbrs; }
                                if (!atZmax) { delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); ++nbrs; }
                                if (!atZmin) { delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1); ++nbrs; }

                                // check if iX, iY, or iZ is located on the subdomain boundaries (MPI)
                                plint bdry_dir = 0;
                                if ( absX >= bdryGap && absX <= (nx-1-bdryGap) ) {
                                    if (iX0 == domain.x1) { bdry_dir = 1; }      // +x boundary
                                    else if (iX0 == domain.x0) { bdry_dir = 2; } // -x boundary
                                    else if (iY0 == domain.y1) { bdry_dir = 3; } // +y boundary
                                    else if (iY0 == domain.y0) { bdry_dir = 4; } // -y boundary
                                    else if (iZ0 == domain.z1) { bdry_dir = 5; } // +z boundary
                                    else if (iZ0 == domain.z0) { bdry_dir = 6; } // -z boundary
                                }

                                // search neighbors for redistribution
                                std::vector<plint> nbrsLocMask; plint nbrTlen = 0;
                                for (plint iT=0; iT<nbrs; ++iT) {
                                    plint delx = delXYZ[iT*3], dely = delXYZ[iT*3+1], delz = delXYZ[iT*3+2];
                                    plint nbrmask = util::roundToInt( lattices[maskLloc]->get(iXm+delx,iYm+dely,iZm+delz).computeDensity() );
                                    // exclude wall boundaries
                                    if (nbrmask != bb && nbrmask != solid) {
                                        ++nbrTlen;
                                        nbrsLocMask.push_back(delx); nbrsLocMask.push_back(dely); nbrsLocMask.push_back(delz); nbrsLocMask.push_back(nbrmask);
                                    }
                                }
                                // shuffle a number vector to randomly select a neighbor
                                std::vector<plint> randArray;
                                if (nbrTlen > 1) {
                                    for (plint iR=0; iR<nbrTlen; ++iR) { randArray.push_back(iR); }
                                    std::random_shuffle(randArray.begin(), randArray.end());
                                }
                                // first trial to redistribute excess biomass
                                bool chk = 0;
                                if (nbrTlen > 0) {
                                    for (plint iT=0; iT<nbrTlen; ++iT) {
                                        // select a neighboring grid cell
                                        plint randLoc = 0;
                                        if (nbrTlen > 1) { randLoc = randArray[iT]; }
                                        plint delx = nbrsLocMask[4*randLoc], dely = nbrsLocMask[4*randLoc+1], delz = nbrsLocMask[4*randLoc+2], nbrmask = nbrsLocMask[4*randLoc+3];
                                        T nbrbMt = 0;
                                        bool nbrmaskflag = 0;
                                        if (nbrmask == bb || nbrmask == solid) { nbrmaskflag = 1; }
                                        else { for (size_t iP=0; iP<pore.size(); ++iP) { if ( nbrmask == pore[iP] ) { nbrmaskflag = 1; break; } } }
                                        if (nbrmaskflag == 0) { nbrbMt = lattices[bMtLloc]->get(iXt+delx,iYt+dely,iZt+delz).computeDensity(); }
                                        if (nbrbMt < Bmax) { // else, select a new neighbor
                                            // redefine the push direction indicator based on the direction of the chosen neighbor
                                            plint push_dir = 0;
                                            if (bdry_dir > 0) {
                                                if ( delx == 1 && bdry_dir == 1 ) { push_dir = 1; }
                                                else if ( delx == -1 && bdry_dir == 2 ) { push_dir = 2; }
                                                else if ( dely == 1 && bdry_dir == 3 ) { push_dir = 3; }
                                                else if ( dely == -1 && bdry_dir == 4 ) { push_dir = 4; }
                                                else if ( delz == 1 && bdry_dir == 5 ) { push_dir = 5; }
                                                else if ( delz == -1 && bdry_dir == 6 ) { push_dir = 6; }
                                                else { push_dir = 0; }
                                            }
                                            T hold_capacity = Bmax - nbrbMt;
                                            T partial_bMassT = bMd;
                                            if (bMd > hold_capacity) { partial_bMassT = hold_capacity; bMd -= hold_capacity; }
                                            else { bMd = 0; chk = 1; }
                                            for (plint iM=0; iM<numbM; ++iM) {
                                                Array<T,7> g;
                                                plint iXb1 = iX0+vec_offset[iM].x, iYb1 = iY0+vec_offset[iM].y, iZb1 = iZ0+vec_offset[iM].z; // original biomass lattice
                                                T shove_bmass = ( lattices[iM]->get(iXb1,iYb1,iZb1).computeDensity() ) / bMt * partial_bMassT;
                                                if (shove_bmass > thrd) {
                                                    lattices[iM]->get(iXb1,iYb1,iZb1).getPopulations(g);
                                                    g[0]-=(T) (shove_bmass)/4; g[1]-=(T) (shove_bmass)/8; g[2]-=(T) (shove_bmass)/8; g[3]-=(T) (shove_bmass)/8; g[4]-=(T) (shove_bmass)/8; g[5]-=(T) (shove_bmass)/8; g[6]-=(T) (shove_bmass)/8;
                                                    lattices[iM]->get(iXb1,iYb1,iZb1).setPopulations(g);
                                                    // Update storage lattices if the selected neighbor location is outside of subdomain boundaries
                                                    if ( push_dir > 0 ) {
                                                        plint iXb2 = iX0+vec_offset[iM+numbM].x, iYb2 = iY0+vec_offset[iM+numbM].y, iZb2 = iZ0+vec_offset[iM+numbM].z;
                                                        g[0]=(T) (shove_bmass-1)/4; g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T) (shove_bmass-1)/8;
                                                        lattices[iM+numbM]->get(iXb2,iYb2,iZb2).setPopulations(g);
                                                        lattices[iM+numbM]->get(iXb2,iYb2,iZb2).getDynamics().setOmega(push_dir);
                                                    }
                                                    // Directly update biomass lattices if the selected neighbor location is inside of subdomain boundaries
                                                    else {
                                                        lattices[iM]->get(iXb1+delx,iYb1+dely,iZb1+delz).getPopulations(g);
                                                        g[0]+=(T) (shove_bmass)/4; g[1]+=(T) (shove_bmass)/8; g[2]+=(T) (shove_bmass)/8; g[3]+=(T) (shove_bmass)/8; g[4]+=(T) (shove_bmass)/8; g[5]+=(T) (shove_bmass)/8; g[6]+=(T) (shove_bmass)/8;
                                                        lattices[iM]->get(iXb1+delx,iYb1+dely,iZb1+delz).setPopulations(g);
                                                    }
                                                }
                                            }
                                            bMt -= partial_bMassT;
                                        }
                                        // break out of the for loop if the excess biomass has been redistributed
                                        if (chk == 1) { break; }
                                    }
                                }
                                else {
                                    // [CA-ROBUST] no non-wall neighbour (cell fully encased by wall/solid, e.g. sealed by
                                    // precipitated FeS): the excess cannot be redistributed. Leave this cell's biomass in
                                    // place at the cap and move on to the next voxel instead of aborting the whole run.
                                }

                                // redistribute the remaining biomass (this is the most time consuming part)
                                if (chk == 0) {
                                    /* use the age lattice */
                                    std::vector<T> tmp1Vector;
                                    plint tmp1Len = 0;
                                    plint iXd = iX0+vec_offset[distLloc].x, iYd = iY0+vec_offset[distLloc].y, iZd = iZ0+vec_offset[distLloc].z;
                                    plint id0 = util::roundToInt( lattices[distLloc]->get(iXd,iYd,iZd).computeDensity() );
                                    for (plint iT=0; iT<nbrTlen; ++iT) {
                                        plint delx = nbrsLocMask[iT*4], dely = nbrsLocMask[iT*4+1], delz = nbrsLocMask[iT*4+2];
                                        plint id1 = util::roundToInt( lattices[distLloc]->get(iXd+delx,iYd+dely,iZd+delz).computeDensity() );
                                        if ( id0 > id1 ) { tmp1Vector.push_back(delx); tmp1Vector.push_back(dely); tmp1Vector.push_back(delz); ++tmp1Len; }
                                    }
                                    plint delx, dely, delz, push_dir=0;
                                    if ( tmp1Len > 1) {
                                        plint randLoc = std::rand() % tmp1Len;
                                        delx = tmp1Vector[randLoc*3]; dely = tmp1Vector[randLoc*3+1]; delz = tmp1Vector[randLoc*3+2];
                                    }
                                    else if (tmp1Len == 1) { delx = tmp1Vector[0]; dely = tmp1Vector[1]; delz = tmp1Vector[2]; }
                                    else { // tmp1Len == 0, purely random redistribution
                                        plint randLoc=0;
                                        if (nbrTlen > 1) { randLoc = std::rand() % nbrTlen; }
                                        delx = nbrsLocMask[randLoc*4]; dely = nbrsLocMask[randLoc*4+1]; delz = nbrsLocMask[randLoc*4+2];
                                    }
                                    if (bdry_dir > 0) {
                                        if ( delx == 1 && bdry_dir == 1 ) { push_dir = 1; }
                                        else if ( delx == -1 && bdry_dir == 2 ) { push_dir = 2; }
                                        else if ( dely == 1 && bdry_dir == 3 ) { push_dir = 3; }
                                        else if ( dely == -1 && bdry_dir == 4 ) { push_dir = 4; }
                                        else if ( delz == 1 && bdry_dir == 5 ) { push_dir = 5; }
                                        else if ( delz == -1 && bdry_dir == 6 ) { push_dir = 6; }
                                        else { push_dir = 0; }
                                    }

                                    for (plint iM=0; iM<numbM; ++iM) {
                                        Array<T,7> g;
                                        plint iXb1 = iX0+vec_offset[iM].x, iYb1 = iY0+vec_offset[iM].y, iZb1 = iZ0+vec_offset[iM].z; // original biomass lattice
                                        T shove_bmass = ( lattices[iM]->get(iXb1,iYb1,iZb1).computeDensity() ) / bMt * bMd;
                                        lattices[iM]->get(iXb1,iYb1,iZb1).getPopulations(g);
                                        g[0]-=(T) (shove_bmass)/4; g[1]-=(T) (shove_bmass)/8; g[2]-=(T) (shove_bmass)/8; g[3]-=(T) (shove_bmass)/8; g[4]-=(T) (shove_bmass)/8; g[5]-=(T) (shove_bmass)/8; g[6]-=(T) (shove_bmass)/8;
                                        lattices[iM]->get(iXb1,iYb1,iZb1).setPopulations(g);
                                        // Update storage lattices if the selected neighbor location is outside of subdomain boundaries
                                        if ( push_dir > 0 ) {
                                            plint iXb2 = iX0+vec_offset[iM+numbM].x, iYb2 = iY0+vec_offset[iM+numbM].y, iZb2 = iZ0+vec_offset[iM+numbM].z;
                                            g[0]=(T) (shove_bmass-1)/4; g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T) (shove_bmass-1)/8;
                                            lattices[iM+numbM]->get(iXb2,iYb2,iZb2).setPopulations(g);
                                            lattices[iM+numbM]->get(iXb2,iYb2,iZb2).getDynamics().setOmega(push_dir);
                                        }
                                        // Directly update biomass lattices if the selected neighbor location is inside of subdomain boundaries
                                        else {
                                            lattices[iM]->get(iXb1+delx,iYb1+dely,iZb1+delz).getPopulations(g);
                                            g[0]+=(T) (shove_bmass)/4; g[1]+=(T) (shove_bmass)/8; g[2]+=(T) (shove_bmass)/8; g[3]+=(T) (shove_bmass)/8; g[4]+=(T) (shove_bmass)/8; g[5]+=(T) (shove_bmass)/8; g[6]+=(T) (shove_bmass)/8;
                                            lattices[iM]->get(iXb1+delx,iYb1+dely,iZb1+delz).setPopulations(g);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Don't apply to envelope, because nearest neighbors need to be accessed.
        return BlockDomain::bulk;
    }
    virtual pushExcessBiomass3D<T,Descriptor>* clone() const {
        return new pushExcessBiomass3D<T,Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        plint numbM = (length-3)/2;
        for (plint iB=0; iB<numbM; ++iB) {
            modified[iB] = modif::staticVariables;
            modified[iB+numbM] = modif::allVariables;
        }
        modified[length-3] = modif::nothing;
        modified[length-2] = modif::nothing;
        modified[length-1] = modif::nothing;
    }
private:
    T Bmax;
    plint nx, ny, nz, bdryGap, length, solid, bb;
    std::vector<plint> pore;
};

// redistribute the excessive biomass (push half)
template<typename T, template<typename U> class Descriptor>
class halfPushExcessBiomass3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    halfPushExcessBiomass3D(T Bmax_, plint nx_, plint ny_, plint nz_, plint bdryGap_, plint length_, plint solid_, plint bb_, std::vector<plint> pore_)
    : Bmax(Bmax_), nx(nx_), ny(ny_), nz(nz_), bdryGap(bdryGap_), length(length_), solid(solid_), bb(bb_), pore(pore_)
    {}
    // lattices[0~(#ofbM-1)] = original biomass lattices
    // lattices[#ofbM~(len-3)] = copy biomass lattices
    // lattices[len-3] = total biomass lattice
    // lattices[len-2] = mask lattice
    // lattices[len-1] = dist lattice
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        std::vector<Dot3D> vec_offset;
        plint distLloc = length-1, maskLloc = length-2, bMtLloc = length-3, numbM = (length-3)/2;
        for (plint iL=0; iL<length; ++iL) {
            vec_offset.push_back(computeRelativeDisplacement(*lattices[0],*lattices[iL]));
        }
        Dot3D absoluteOffset = lattices[0]->getLocation();
        for (plint iX0=domain.x0; iX0<=domain.x1; ++iX0) {
            plint iXm = iX0+vec_offset[maskLloc].x;
            for (plint iY0=domain.y0; iY0<=domain.y1; ++iY0) {
                plint iYm = iY0+vec_offset[maskLloc].y;
                for (plint iZ0=domain.z0; iZ0<=domain.z1; ++iZ0) {
                    plint iZm = iZ0+vec_offset[maskLloc].z;
                    plint mask = util::roundToInt( lattices[maskLloc]->get(iXm,iYm,iZm).computeDensity() );
                    bool maskflag = 0;
                    if (mask == bb || mask == solid) { maskflag = 1; }
                    else {
                        for (size_t iP=0; iP<pore.size(); ++iP) {
                            if ( mask == pore[iP] ) {
                                maskflag = 1;
                                break;
                            }
                        }
                    }
                    if ( maskflag == 0 ) {
                        plint iXt = iX0 + vec_offset[bMtLloc].x, iYt = iY0 + vec_offset[bMtLloc].y, iZt = iZ0 + vec_offset[bMtLloc].z;
                        T bMt = lattices[bMtLloc]->get(iXt,iYt,iZt).computeDensity();
                        if (bMt > Bmax) {
                            T bMd = bMt*0.5;
                            plint absX = iX0 + absoluteOffset.x, absY = iY0 + absoluteOffset.y, absZ = iZ0 + absoluteOffset.z;
                            std::vector<plint> delXYZ; plint nbrs = 0;
                            
                            // Determine neighbors based on boundary position (3D version)
                            bool atXmin = (absX == bdryGap);
                            bool atXmax = (absX == nx-1-bdryGap);
                            bool atYmin = (absY == 0);
                            bool atYmax = (absY == ny-1);
                            bool atZmin = (absZ == 0);
                            bool atZmax = (absZ == nz-1);
                            
                            // Add valid neighbors
                            if (!atXmax) { delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); ++nbrs; }
                            if (!atXmin) { delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); ++nbrs; }
                            if (!atYmax) { delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); ++nbrs; }
                            if (!atYmin) { delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0); ++nbrs; }
                            if (!atZmax) { delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); ++nbrs; }
                            if (!atZmin) { delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1); ++nbrs; }

                            // check if iX, iY, or iZ is located on the subdomain boundaries (MPI)
                            plint bdry_dir = 0;
                            if ( absX >= bdryGap && absX <= (nx-1-bdryGap) ) {
                                if (iX0 == domain.x1) { bdry_dir = 1; }      // +x boundary
                                else if (iX0 == domain.x0) { bdry_dir = 2; } // -x boundary
                                else if (iY0 == domain.y1) { bdry_dir = 3; } // +y boundary
                                else if (iY0 == domain.y0) { bdry_dir = 4; } // -y boundary
                                else if (iZ0 == domain.z1) { bdry_dir = 5; } // +z boundary
                                else if (iZ0 == domain.z0) { bdry_dir = 6; } // -z boundary
                            }

                            // search neighbors for redistribution
                            std::vector<plint> nbrsLocMask; plint nbrTlen = 0;
                            for (plint iT=0; iT<nbrs; ++iT) {
                                plint delx = delXYZ[iT*3], dely = delXYZ[iT*3+1], delz = delXYZ[iT*3+2];
                                plint nbrmask = util::roundToInt( lattices[maskLloc]->get(iXm+delx,iYm+dely,iZm+delz).computeDensity() );
                                // exclude wall boundaries
                                if (nbrmask != bb && nbrmask != solid) {
                                    ++nbrTlen;
                                    nbrsLocMask.push_back(delx); nbrsLocMask.push_back(dely); nbrsLocMask.push_back(delz); nbrsLocMask.push_back(nbrmask);
                                }
                            }
                            // shuffle a number vector to randomly select a neighbor
                            std::vector<plint> randArray;
                            if (nbrTlen > 1) {
                                for (plint iR=0; iR<nbrTlen; ++iR) { randArray.push_back(iR); }
                                std::random_shuffle(randArray.begin(), randArray.end());
                            }
                            // first trial to redistribute excess biomass
                            bool chk = 0;
                            if (nbrTlen > 0) {
                                for (plint iT=0; iT<nbrTlen; ++iT) {
                                    // select a neighboring grid cell
                                    plint randLoc = 0;
                                    if (nbrTlen > 1) { randLoc = randArray[iT]; }
                                    plint delx = nbrsLocMask[4*randLoc], dely = nbrsLocMask[4*randLoc+1], delz = nbrsLocMask[4*randLoc+2], nbrmask = nbrsLocMask[4*randLoc+3];
                                    T nbrbMt = 0;
                                    bool nbrmaskflag = 0;
                                    if (nbrmask == bb || nbrmask == solid) { nbrmaskflag = 1; }
                                    else { for (size_t iP=0; iP<pore.size(); ++iP) { if ( nbrmask == pore[iP] ) { nbrmaskflag = 1; break; } } }
                                    if (nbrmaskflag == 0) { nbrbMt = lattices[bMtLloc]->get(iXt+delx,iYt+dely,iZt+delz).computeDensity(); }
                                    if (nbrbMt < Bmax) { // else, select a new neighbor
                                        // redefine the push direction indicator based on the direction of the chosen neighbor
                                        plint push_dir = 0;
                                        if (bdry_dir > 0) {
                                            if ( delx == 1 && bdry_dir == 1 ) { push_dir = 1; }
                                            else if ( delx == -1 && bdry_dir == 2 ) { push_dir = 2; }
                                            else if ( dely == 1 && bdry_dir == 3 ) { push_dir = 3; }
                                            else if ( dely == -1 && bdry_dir == 4 ) { push_dir = 4; }
                                            else if ( delz == 1 && bdry_dir == 5 ) { push_dir = 5; }
                                            else if ( delz == -1 && bdry_dir == 6 ) { push_dir = 6; }
                                            else { push_dir = 0; }
                                        }
                                        T hold_capacity = Bmax - nbrbMt;
                                        T partial_bMassT = bMd;
                                        if (bMd > hold_capacity) { partial_bMassT = hold_capacity; bMd -= hold_capacity; }
                                        else { bMd = 0; chk = 1; }
                                        for (plint iM=0; iM<numbM; ++iM) {
                                            Array<T,7> g;
                                            plint iXb1 = iX0+vec_offset[iM].x, iYb1 = iY0+vec_offset[iM].y, iZb1 = iZ0+vec_offset[iM].z;
                                            T shove_bmass = ( lattices[iM]->get(iXb1,iYb1,iZb1).computeDensity() ) / bMt * partial_bMassT;
                                            if (shove_bmass > thrd) {
                                                lattices[iM]->get(iXb1,iYb1,iZb1).getPopulations(g);
                                                g[0]-=(T) (shove_bmass)/4; g[1]-=(T) (shove_bmass)/8; g[2]-=(T) (shove_bmass)/8; g[3]-=(T) (shove_bmass)/8; g[4]-=(T) (shove_bmass)/8; g[5]-=(T) (shove_bmass)/8; g[6]-=(T) (shove_bmass)/8;
                                                lattices[iM]->get(iXb1,iYb1,iZb1).setPopulations(g);
                                                // Update storage lattices if the selected neighbor location is outside of subdomain boundaries
                                                if ( push_dir > 0 ) {
                                                    plint iXb2 = iX0+vec_offset[iM+numbM].x, iYb2 = iY0+vec_offset[iM+numbM].y, iZb2 = iZ0+vec_offset[iM+numbM].z;
                                                    g[0]=(T) (shove_bmass-1)/4; g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T) (shove_bmass-1)/8;
                                                    lattices[iM+numbM]->get(iXb2,iYb2,iZb2).setPopulations(g);
                                                    lattices[iM+numbM]->get(iXb2,iYb2,iZb2).getDynamics().setOmega(push_dir);
                                                }
                                                // Directly update biomass lattices if the selected neighbor location is inside of subdomain boundaries
                                                else {
                                                    lattices[iM]->get(iXb1+delx,iYb1+dely,iZb1+delz).getPopulations(g);
                                                    g[0]+=(T) (shove_bmass)/4; g[1]+=(T) (shove_bmass)/8; g[2]+=(T) (shove_bmass)/8; g[3]+=(T) (shove_bmass)/8; g[4]+=(T) (shove_bmass)/8; g[5]+=(T) (shove_bmass)/8; g[6]+=(T) (shove_bmass)/8;
                                                    lattices[iM]->get(iXb1+delx,iYb1+dely,iZb1+delz).setPopulations(g);
                                                }
                                            }
                                        }
                                        bMt -= partial_bMassT;
                                    }
                                    // break out of the for loop if the excess biomass has been redistributed
                                    if (chk == 1) { break; }
                                }
                            }
                            else {
                                // [CA-ROBUST] no non-wall neighbour (cell fully encased by wall/solid, e.g. sealed by
                                // precipitated FeS): the excess cannot be redistributed. Leave this cell's biomass in
                                // place at the cap and move on to the next voxel instead of aborting the whole run.
                            }

                            // redistribute the remaining biomass (this is the most time consuming part)
                            plint push_dir = 0, delx = 0, dely = 0, delz = 0;
                            if (chk == 0) {
                                /* use the distance lattice */
                                std::vector<T> tmp1Vector;
                                plint tmp1Len = 0;
                                plint iXd = iX0+vec_offset[distLloc].x, iYd = iY0+vec_offset[distLloc].y, iZd = iZ0+vec_offset[distLloc].z;
                                plint id0 = util::roundToInt( lattices[distLloc]->get(iXd,iYd,iZd).computeDensity() );
                                for (plint iT=0; iT<nbrTlen; ++iT) {
                                    delx = nbrsLocMask[iT*4]; dely = nbrsLocMask[iT*4+1]; delz = nbrsLocMask[iT*4+2];
                                    plint id1 = util::roundToInt( lattices[distLloc]->get(iXd+delx,iYd+dely,iZd+delz).computeDensity() );
                                    if ( id0 >= id1 ) { tmp1Vector.push_back(delx); tmp1Vector.push_back(dely); tmp1Vector.push_back(delz); ++tmp1Len; }
                                }
                                if ( tmp1Len > 1) {
                                    plint randLoc = std::rand() % tmp1Len;
                                    delx = tmp1Vector[randLoc*3]; dely = tmp1Vector[randLoc*3+1]; delz = tmp1Vector[randLoc*3+2];
                                }
                                else if (tmp1Len == 1) { delx = tmp1Vector[0]; dely = tmp1Vector[1]; delz = tmp1Vector[2]; }
                                else { // tmp1Len == 0, purely random redistribution
                                    plint randLoc=0;
                                    if (nbrTlen > 1) { randLoc = std::rand() % nbrTlen; }
                                    delx = nbrsLocMask[randLoc*4]; dely = nbrsLocMask[randLoc*4+1]; delz = nbrsLocMask[randLoc*4+2];
                                }
                                if (bdry_dir > 0) {
                                    if ( delx == 1 && bdry_dir == 1 ) { push_dir = 1; }
                                    else if ( delx == -1 && bdry_dir == 2 ) { push_dir = 2; }
                                    else if ( dely == 1 && bdry_dir == 3 ) { push_dir = 3; }
                                    else if ( dely == -1 && bdry_dir == 4 ) { push_dir = 4; }
                                    else if ( delz == 1 && bdry_dir == 5 ) { push_dir = 5; }
                                    else if ( delz == -1 && bdry_dir == 6 ) { push_dir = 6; }
                                    else { push_dir = 0; }
                                }
                            }
                            for (plint iM=0; iM<numbM; ++iM) {
                                Array<T,7> g;
                                plint iXb1 = iX0+vec_offset[iM].x, iYb1 = iY0+vec_offset[iM].y, iZb1 = iZ0+vec_offset[iM].z; // original biomass lattice
                                T shove_bmass = ( lattices[iM]->get(iXb1,iYb1,iZb1).computeDensity() ) / bMt * bMd;
                                lattices[iM]->get(iXb1,iYb1,iZb1).getPopulations(g);
                                g[0]-=(T) (shove_bmass)/4; g[1]-=(T) (shove_bmass)/8; g[2]-=(T) (shove_bmass)/8; g[3]-=(T) (shove_bmass)/8; g[4]-=(T) (shove_bmass)/8; g[5]-=(T) (shove_bmass)/8; g[6]-=(T) (shove_bmass)/8;
                                lattices[iM]->get(iXb1,iYb1,iZb1).setPopulations(g);
                                // Update storage lattices if the selected neighbor location is outside of subdomain boundaries
                                if ( push_dir > 0 ) {
                                    plint iXb2 = iX0+vec_offset[iM+numbM].x, iYb2 = iY0+vec_offset[iM+numbM].y, iZb2 = iZ0+vec_offset[iM+numbM].z;
                                    g[0]=(T) (shove_bmass-1)/4; g[1]=g[2]=g[3]=g[4]=g[5]=g[6]=(T) (shove_bmass-1)/8;
                                    lattices[iM+numbM]->get(iXb2,iYb2,iZb2).setPopulations(g);
                                    lattices[iM+numbM]->get(iXb2,iYb2,iZb2).getDynamics().setOmega(push_dir);
                                }
                                // Directly update biomass lattices if the selected neighbor location is inside of subdomain boundaries
                                else {
                                    lattices[iM]->get(iXb1+delx,iYb1+dely,iZb1+delz).getPopulations(g);
                                    g[0]+=(T) (shove_bmass)/4; g[1]+=(T) (shove_bmass)/8; g[2]+=(T) (shove_bmass)/8; g[3]+=(T) (shove_bmass)/8; g[4]+=(T) (shove_bmass)/8; g[5]+=(T) (shove_bmass)/8; g[6]+=(T) (shove_bmass)/8;
                                    lattices[iM]->get(iXb1+delx,iYb1+dely,iZb1+delz).setPopulations(g);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Don't apply to envelope, because nearest neighbors need to be accessed.
        return BlockDomain::bulk;
    }
    virtual halfPushExcessBiomass3D<T,Descriptor>* clone() const {
        return new halfPushExcessBiomass3D<T,Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        plint numbM = (length-3)/2;
        for (plint iB=0; iB<numbM; ++iB) {
            modified[iB] = modif::staticVariables;
            modified[iB+numbM] = modif::allVariables;
        }
        modified[length-3] = modif::nothing;
        modified[length-2] = modif::nothing;
        modified[length-1] = modif::nothing;
    }
private:
    T Bmax;
    plint nx, ny, nz, bdryGap, length, solid, bb;
    std::vector<plint> pore;
};

// redistribute the excessive biomass (pull)
template<typename T, template<typename U> class Descriptor>
class pullExcessBiomass3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    pullExcessBiomass3D(plint nx_, plint ny_, plint nz_, plint bdryGap_, plint length_)
    : nx(nx_), ny(ny_), nz(nz_), bdryGap(bdryGap_), length(length_)
    {}
    // lattices[0~(#ofbM-1)] = original biomass lattices
    // lattices[#ofbM~(len-3)] = copy biomass lattices
    // lattices[len-3] = total biomass lattice
    // lattices[len-2] = mask lattice
    // lattices[len-1] = dist lattice
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        std::vector<Dot3D> vec_offset;
        plint numbM = (length-2)/2;
        for (plint iL=0; iL<length; ++iL) {
            vec_offset.push_back(computeRelativeDisplacement(*lattices[0],*lattices[iL]));
        }
        Dot3D absoluteOffset = lattices[0]->getLocation();
        
        // +x boundary (pull from -x neighbor)
        plint iX0 = domain.x0, absX = iX0 + absoluteOffset.x;
        if ( absX >= bdryGap ) {
            for (plint iY0=domain.y0; iY0<=domain.y1; ++iY0) {
                for (plint iZ0=domain.z0; iZ0<=domain.z1; ++iZ0) {
                    for (plint iM=0; iM<numbM; ++iM) {
                        plint iXc = iX0 + vec_offset[iM+numbM].x, iYc = iY0 + vec_offset[iM+numbM].y, iZc = iZ0 + vec_offset[iM+numbM].z;
                        plint dir_id = util::roundToInt( lattices[iM+numbM]->get((iXc-1),iYc,iZc).getDynamics().getOmega() );
                        if (dir_id == 1) {
                            T nbrbM = lattices[iM+numbM]->get((iXc-1),iYc,iZc).computeDensity();
                            if (nbrbM > thrd) {
                                Array<T,7> g;
                                plint iXb = iX0 + vec_offset[iM].x, iYb = iY0 + vec_offset[iM].y, iZb = iZ0 + vec_offset[iM].z;
                                lattices[iM]->get(iXb,iYb,iZb).getPopulations(g);
                                g[0]+=(T) (nbrbM)/4; g[1]+=(T) (nbrbM)/8; g[2]+=(T) (nbrbM)/8; g[3]+=(T) (nbrbM)/8; g[4]+=(T) (nbrbM)/8; g[5]+=(T) (nbrbM)/8; g[6]+=(T) (nbrbM)/8;
                                lattices[iM]->get(iXb,iYb,iZb).setPopulations(g);
                            }
                        }
                    }
                }
            }
        }
        
        // -x boundary (pull from +x neighbor)
        iX0 = domain.x1; absX = iX0 + absoluteOffset.x;
        if ( absX <= (nx-1-bdryGap) ) {
            for (plint iY0=domain.y0; iY0<=domain.y1; ++iY0) {
                for (plint iZ0=domain.z0; iZ0<=domain.z1; ++iZ0) {
                    for (plint iM=0; iM<numbM; ++iM) {
                        plint iXc = iX0 + vec_offset[iM+numbM].x, iYc = iY0 + vec_offset[iM+numbM].y, iZc = iZ0 + vec_offset[iM+numbM].z;
                        plint dir_id = util::roundToInt( lattices[iM+numbM]->get((iXc+1),iYc,iZc).getDynamics().getOmega() );
                        if (dir_id == 2) {
                            T nbrbM = lattices[iM+numbM]->get((iXc+1),iYc,iZc).computeDensity();
                            if (nbrbM > thrd) {
                                Array<T,7> g;
                                plint iXb = iX0 + vec_offset[iM].x, iYb = iY0 + vec_offset[iM].y, iZb = iZ0 + vec_offset[iM].z;
                                lattices[iM]->get(iXb,iYb,iZb).getPopulations(g);
                                g[0]+=(T) (nbrbM)/4; g[1]+=(T) (nbrbM)/8; g[2]+=(T) (nbrbM)/8; g[3]+=(T) (nbrbM)/8; g[4]+=(T) (nbrbM)/8; g[5]+=(T) (nbrbM)/8; g[6]+=(T) (nbrbM)/8;
                                lattices[iM]->get(iXb,iYb,iZb).setPopulations(g);
                            }
                        }
                    }
                }
            }
        }
        
        // +y boundary (pull from -y neighbor)
        plint iY0 = domain.y0, absY = iY0 + absoluteOffset.y;
        if ( absY > 0 ) {
            for (iX0=domain.x0; iX0<=domain.x1; ++iX0) {
                for (plint iZ0=domain.z0; iZ0<=domain.z1; ++iZ0) {
                    for (plint iM=0; iM<numbM; ++iM) {
                        plint iXc = iX0 + vec_offset[iM+numbM].x, iYc = iY0 + vec_offset[iM+numbM].y, iZc = iZ0 + vec_offset[iM+numbM].z;
                        plint dir_id = util::roundToInt( lattices[iM+numbM]->get(iXc,(iYc-1),iZc).getDynamics().getOmega() );
                        if (dir_id == 3) {
                            T nbrbM = lattices[iM+numbM]->get(iXc,(iYc-1),iZc).computeDensity();
                            if (nbrbM > thrd) {
                                Array<T,7> g;
                                plint iXb = iX0 + vec_offset[iM].x, iYb = iY0 + vec_offset[iM].y, iZb = iZ0 + vec_offset[iM].z;
                                lattices[iM]->get(iXb,iYb,iZb).getPopulations(g);
                                g[0]+=(T) (nbrbM)/4; g[1]+=(T) (nbrbM)/8; g[2]+=(T) (nbrbM)/8; g[3]+=(T) (nbrbM)/8; g[4]+=(T) (nbrbM)/8; g[5]+=(T) (nbrbM)/8; g[6]+=(T) (nbrbM)/8;
                                lattices[iM]->get(iXb,iYb,iZb).setPopulations(g);
                            }
                        }
                    }
                }
            }
        }

        // -y boundary (pull from +y neighbor)
        iY0 = domain.y1; absY = iY0 + absoluteOffset.y;
        if ( absY < (ny-1) ) {
            for (iX0=domain.x0; iX0<=domain.x1; ++iX0) {
                for (plint iZ0=domain.z0; iZ0<=domain.z1; ++iZ0) {
                    for (plint iM=0; iM<numbM; ++iM) {
                        plint iXc = iX0 + vec_offset[iM+numbM].x, iYc = iY0 + vec_offset[iM+numbM].y, iZc = iZ0 + vec_offset[iM+numbM].z;
                        plint dir_id = util::roundToInt( lattices[iM+numbM]->get(iXc,(iYc+1),iZc).getDynamics().getOmega() );
                        if (dir_id == 4) {
                            T nbrbM = lattices[iM+numbM]->get(iXc,(iYc+1),iZc).computeDensity();
                            if (nbrbM > thrd) {
                                Array<T,7> g;
                                plint iXb = iX0 + vec_offset[iM].x, iYb = iY0 + vec_offset[iM].y, iZb = iZ0 + vec_offset[iM].z;
                                lattices[iM]->get(iXb,iYb,iZb).getPopulations(g);
                                g[0]+=(T) (nbrbM)/4; g[1]+=(T) (nbrbM)/8; g[2]+=(T) (nbrbM)/8; g[3]+=(T) (nbrbM)/8; g[4]+=(T) (nbrbM)/8; g[5]+=(T) (nbrbM)/8; g[6]+=(T) (nbrbM)/8;
                                lattices[iM]->get(iXb,iYb,iZb).setPopulations(g);
                            }
                        }
                    }
                }
            }
        }
        
        // +z boundary (pull from -z neighbor)
        plint iZ0 = domain.z0, absZ = iZ0 + absoluteOffset.z;
        if ( absZ > 0 ) {
            for (iX0=domain.x0; iX0<=domain.x1; ++iX0) {
                for (iY0=domain.y0; iY0<=domain.y1; ++iY0) {
                    for (plint iM=0; iM<numbM; ++iM) {
                        plint iXc = iX0 + vec_offset[iM+numbM].x, iYc = iY0 + vec_offset[iM+numbM].y, iZc = iZ0 + vec_offset[iM+numbM].z;
                        plint dir_id = util::roundToInt( lattices[iM+numbM]->get(iXc,iYc,(iZc-1)).getDynamics().getOmega() );
                        if (dir_id == 5) {
                            T nbrbM = lattices[iM+numbM]->get(iXc,iYc,(iZc-1)).computeDensity();
                            if (nbrbM > thrd) {
                                Array<T,7> g;
                                plint iXb = iX0 + vec_offset[iM].x, iYb = iY0 + vec_offset[iM].y, iZb = iZ0 + vec_offset[iM].z;
                                lattices[iM]->get(iXb,iYb,iZb).getPopulations(g);
                                g[0]+=(T) (nbrbM)/4; g[1]+=(T) (nbrbM)/8; g[2]+=(T) (nbrbM)/8; g[3]+=(T) (nbrbM)/8; g[4]+=(T) (nbrbM)/8; g[5]+=(T) (nbrbM)/8; g[6]+=(T) (nbrbM)/8;
                                lattices[iM]->get(iXb,iYb,iZb).setPopulations(g);
                            }
                        }
                    }
                }
            }
        }

        // -z boundary (pull from +z neighbor)
        iZ0 = domain.z1; absZ = iZ0 + absoluteOffset.z;
        if ( absZ < (nz-1) ) {
            for (iX0=domain.x0; iX0<=domain.x1; ++iX0) {
                for (iY0=domain.y0; iY0<=domain.y1; ++iY0) {
                    for (plint iM=0; iM<numbM; ++iM) {
                        plint iXc = iX0 + vec_offset[iM+numbM].x, iYc = iY0 + vec_offset[iM+numbM].y, iZc = iZ0 + vec_offset[iM+numbM].z;
                        plint dir_id = util::roundToInt( lattices[iM+numbM]->get(iXc,iYc,(iZc+1)).getDynamics().getOmega() );
                        if (dir_id == 6) {
                            T nbrbM = lattices[iM+numbM]->get(iXc,iYc,(iZc+1)).computeDensity();
                            if (nbrbM > thrd) {
                                Array<T,7> g;
                                plint iXb = iX0 + vec_offset[iM].x, iYb = iY0 + vec_offset[iM].y, iZb = iZ0 + vec_offset[iM].z;
                                lattices[iM]->get(iXb,iYb,iZb).getPopulations(g);
                                g[0]+=(T) (nbrbM)/4; g[1]+=(T) (nbrbM)/8; g[2]+=(T) (nbrbM)/8; g[3]+=(T) (nbrbM)/8; g[4]+=(T) (nbrbM)/8; g[5]+=(T) (nbrbM)/8; g[6]+=(T) (nbrbM)/8;
                                lattices[iM]->get(iXb,iYb,iZb).setPopulations(g);
                            }
                        }
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Don't apply to envelope, because nearest neighbors need to be accessed.
        return BlockDomain::bulk;
    }
    virtual pullExcessBiomass3D<T,Descriptor>* clone() const {
        return new pullExcessBiomass3D<T,Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        plint numbM = (length-3)/2;
        for (plint iB=0; iB<numbM; ++iB) {
            modified[iB] = modif::staticVariables;
            modified[iB+numbM] = modif::nothing;
        }
        modified[length-3] = modif::nothing;
        modified[length-2] = modif::nothing;
        modified[length-1] = modif::nothing;
    }
private:
    plint nx, ny, nz, bdryGap, length;
};

#endif // COMPLAB3D_PROCESSORS_PART1_HH
