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
 *  complab3d_processors_surrogate.hh
 *  --------------------------------------------------------------------------
 *  THE SURROGATE MODEL PROCESSOR
 *
 *  A surrogate model is a cheap stand-in for a full flux balance calculation.
 *  Instead of solving a linear program at every voxel, you train something fast
 *  -- a neural network, a polynomial fit, a lookup table -- on the answers the
 *  linear program gives, and then evaluate that instead.  It is typically two
 *  to three orders of magnitude faster, at the cost of only being trustworthy
 *  inside the range of conditions it was trained on.
 *
 *  This processor prepares the inputs, calls the user's model in
 *  ../surrogateModel.hh, and turns the answer into concentration and biomass
 *  increments.  The model itself is a user-editable file, like
 *  defineKinetics.hh, and changing it requires recompiling.
 *
 *  --------------------------------------------------------------------------
 *  THE INDEXING CONTRACT -- this is the bug fix that matters most
 *
 *  CompLaB 2D built its input arrays over the microbes PRESENT IN THIS VOXEL,
 *  compacted, and the shipped example network then read Fin[0][0], Fin[0][1]
 *  and wrote bioR[0].  So "microbe 0" meant "whichever organism happened to be
 *  listed first in this particular voxel" -- which changes from voxel to voxel.
 *  With more than one surrogate organism the network was applied to an
 *  arbitrary one and every other organism silently fell through to pure decay.
 *
 *  Here the arrays are always FULL LENGTH, indexed by the microbe's number in
 *  CompLaB.xml.  Fin[3] is always microbe3, in every voxel, on every rank.
 *  Rows for organisms that are absent from this voxel are left at zero, and the
 *  user's model is called once per PRESENT organism with its own global id.
 *
 *  --------------------------------------------------------------------------
 *  UNITS
 *      Fin, Fout   mmol / gDW / h    positive = consumed by the organism
 *      bioR        1 / h             specific growth rate
 *      C, B        mol / L           CompLB3D's convention throughout
 *
 *  The conversion between the two lives in MetabolicConfig::fluxToDeltaC() and
 *  growthToDeltaB(), and nowhere else.
 *
 *  --------------------------------------------------------------------------
 *  LATTICE VECTOR LAYOUT   (identical to the FBA processors)
 *
 *      [0            .. subsNum-1]        substrate concentrations   C
 *      [subsNum      .. subsNum+bioNum-1] biomass of the surrogate microbes
 *      [dCloc        .. dCloc+subsNum-1]  substrate increments       dC
 *      [dBloc        .. dBloc+bioNum-1]   biomass increments         dB
 *      [maskLloc]                         mask lattice
 * ============================================================================
 */

#ifndef COMPLAB3D_PROCESSORS_SURROGATE_HH
#define COMPLAB3D_PROCESSORS_SURROGATE_HH

#include "complab3d_metabolic.hh"
#include "../surrogateModel.hh"

template<typename T, template<typename U> class Descriptor>
class run_surrogate3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    /* totalMicrobes is the FULL microbe count from CompLaB.xml, not the number
     * of surrogate microbes -- it sizes the arrays the user's model indexes. */
    run_surrogate3D (plint nx_, plint subsNum_, plint bioNum_, plint totalMicrobes_, T dt_,
                     const std::vector<plint> &globalId_,
                     const MetabolicConfig    *cfg_,
                     const std::vector< std::vector<T> > &vec2_Kc_,
                     const std::vector<T>                &vec1_mu_,
                     plint solid_, plint bb_)
        : nx(nx_), subsNum(subsNum_), bioNum(bioNum_), totalMicrobes(totalMicrobes_), dt(dt_),
          globalId(globalId_), cfg(cfg_), vec2_Kc(vec2_Kc_), vec1_mu(vec1_mu_),
          solid(solid_), bb(bb_),
          dCloc(subsNum_ + bioNum_),
          dBloc(2*subsNum_ + bioNum_),
          maskLloc(2*(subsNum_ + bioNum_))
    {}

    virtual void process (Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices)
    {
        Dot3D absoluteOffset = lattices[0]->getLocation();
        Dot3D maskOffset     = computeRelativeDisplacement(*lattices[0], *lattices[maskLloc]);

        /* loop-invariant: computed once, not once per voxel */
        std::vector<Dot3D> off;
        off.reserve(maskLloc);
        for (plint iT = 0; iT < maskLloc; ++iT)
            off.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT]));

        std::vector<T> conc (subsNum, T());
        std::vector<T> avail(subsNum, T());
        std::vector<T> dC   (subsNum, T());
        std::vector<T> bmass(bioNum,  T());
        std::vector<T> dB   (bioNum,  T());
        std::vector<plint> bLoc;

        /* full-length, global-id-indexed -- see the contract note above */
        std::vector< std::vector<double> > Fin (totalMicrobes, std::vector<double>(subsNum, 0.0));
        std::vector< std::vector<double> > Fout(totalMicrobes, std::vector<double>(subsNum, 0.0));
        std::vector<double> bioR(totalMicrobes, 0.0);

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            const plint absX = iX + absoluteOffset.x;
            if (absX <= 0 || absX >= nx-1) continue;

            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {

                    const plint mask = util::roundToInt(
                        lattices[maskLloc]->get(iX+maskOffset.x, iY+maskOffset.y, iZ+maskOffset.z).computeDensity());
                    if (mask == solid || mask == bb) continue;

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

                    for (plint iS = 0; iS < subsNum; ++iS) {
                        const T c = lattices[iS]->get(iX+off[iS].x, iY+off[iS].y, iZ+off[iS].z).computeDensity();
                        conc[iS] = (c > thrd) ? c : T();
                        dC[iS]   = T();
                        avail[iS] = T();
                    }
                    for (plint iS = 0; iS < subsNum; ++iS)
                        avail[iS] = cfg->useTotals ? cfg->totals.total(conc, iS) : conc[iS];

                    /* ---- build the Michaelis-Menten uptake estimates ----- *
                     * Only the rows of organisms actually present are filled; the rest
                     * stay at zero so a user model that inspects them sees "absent". */
                    for (plint gM = 0; gM < totalMicrobes; ++gM) {
                        std::fill(Fin [gM].begin(), Fin [gM].end(), 0.0);
                        std::fill(Fout[gM].begin(), Fout[gM].end(), 0.0);
                        bioR[gM] = 0.0;
                    }

                    for (size_t k = 0; k < bLoc.size(); ++k) {
                        const plint iB = bLoc[k];
                        const plint gM = globalId[iB];
                        for (plint iS = 0; iS < subsNum; ++iS) {
                            const T mm = (vec2_Kc[gM][iS] + avail[iS] > T())
                                       ? avail[iS] / (vec2_Kc[gM][iS] + avail[iS]) : T();
                            /* [FIX-3D] Vmax comes from <maximum_uptake_flux>, its own field.
                             * The first version read cfg->maxUptake, which holds
                             * <substrate_lower_bounds> -- a hard floor, not a rate -- so the
                             * documented tag was ignored and an unrelated one was used in its
                             * place.  Same predicate as both FBA paths. */
                            T v;
                            if (cfg->vmax[gM][iS] > thrd) {
                                v = cfg->vmax[gM][iS] * mm;          // enzyme-limited, positive = consumed
                            } else {
                                /* no Vmax: the flux that would exhaust the local supply.
                                 * exhaustionFlux is negative (uptake); the surrogate
                                 * convention is positive = consumed. */
                                v = -cfg->exhaustionFlux(avail[iS], bmass[iB], gM, dt) * mm;
                            }

                            /* [FIX-3D] The depletion limiter used to read
                             *     if (maxDraw < 0 && v > -maxDraw) v = -maxDraw;
                             * with v already equal to -maxDraw*mm and mm in [0,1], so
                             * v <= -maxDraw always and the test could never be true.  It
                             * was a no-op in exactly the branch that needs it -- swept over
                             * thousands of states it fired zero times.  Compare against the
                             * unattenuated capacity instead. */
                            const T cap = -cfg->exhaustionFlux(avail[iS], bmass[iB], gM, dt); // >= 0
                            if (cap > T() && v > cap) v = cap;

                            /* [FIX-3D] Honour the same two clamps the FBA paths apply, so a
                             * user who switches a microbe between surrogate and glpk does not
                             * silently get different limits.  maxUptake is stored as a negative
                             * floor; the surrogate convention is positive = consumed. */
                            const T floorMag = -cfg->maxUptake[gM][iS];
                            if (cfg->maxUptake[gM][iS] > (T) -1e29 && floorMag >= T() && v > floorMag)
                                v = floorMag;
                            if (cfg->fix_lower_bounds[iS]) {
                                const T avCap = -cfg->exhaustionFlux(avail[iS], bmass[iB], gM, dt);
                                if (avCap >= T() && v > avCap) v = avCap;
                            }

                            Fin[gM][iS] = (double) v;
                        }
                    }

                    /* Fout starts as a copy of Fin: a user model that only predicts
                     * growth and never touches Fout gets the plain Monod consumption. */
                    Fout = Fin;

                    /* ---- call the user's model, once per present organism -- */
                    for (size_t k = 0; k < bLoc.size(); ++k) {
                        defineSurrogateModel(globalId[bLoc[k]], Fin, Fout, bioR, mask);
                    }

                    /* ---- turn the answer into increments ------------------ */
                    for (size_t k = 0; k < bLoc.size(); ++k) {
                        const plint iB = bLoc[k];
                        const plint gM = globalId[iB];
                        const T mu = (T) bioR[gM];
                        dB[iB] = (std::fabs(mu) < thrd)
                               ? cfg->decayToDeltaB(vec1_mu[gM], bmass[iB], dt)   // [FIX-3D] per second
                               : cfg->growthToDeltaB(mu, bmass[iB], dt);
                        if (dB[iB] + bmass[iB] < T()) dB[iB] = -bmass[iB];
                    }

                    /* [FIX-3D] Substrate budget summed over EVERY organism present, then
                     * clamped once.  The first version clamped inside the per-organism
                     * loop against an `avail` that was never decremented, so two surrogate
                     * organisms in one voxel each took the whole supply and the lattice
                     * went negative on the first step.  Matches the FBA paths now. */
                    for (plint iS = 0; iS < subsNum; ++iS) {
                        if (cfg->fix_concentration[iS]) continue;
                        T draw = T();
                        for (size_t k = 0; k < bLoc.size(); ++k) {
                            const plint iB = bLoc[k];
                            const plint gM = globalId[iB];
                            /* positive Fout = consumed, so the increment is negative */
                            draw += -cfg->fluxToDeltaC((T) Fout[gM][iS], bmass[iB], gM, dt);
                        }
                        if (std::fabs(draw) <= thrd) continue;
                        /* [FIX-3D] Clamp against what is ACTUALLY left, not the pristine
                         * `avail`.  See the identical note in the FBA processors: with
                         * fba_concentration_basis = total, a carrier species that is itself
                         * a transported substrate was debited twice. */
                        const T headroom = avail[iS] + (cfg->useTotals ? T() : dC[iS]);
                        if (draw + headroom < T()) draw = -headroom;
                        if (draw > T()) draw = T();          // consumption only
                        if (cfg->useTotals) cfg->totals.drawDown(conc, iS, -draw, dC);
                        else                dC[iS] += draw;
                    }

                    /* ---- deposit ----------------------------------------- */
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

    virtual run_surrogate3D<T,Descriptor>* clone() const {
        return new run_surrogate3D<T,Descriptor>(*this);
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
    std::vector< std::vector<T> > vec2_Kc;
    std::vector<T> vec1_mu;
    plint solid, bb, dCloc, dBloc, maskLloc;
};

#endif  // COMPLAB3D_PROCESSORS_SURROGATE_HH
