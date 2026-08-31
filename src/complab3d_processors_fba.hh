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
 *  complab3d_processors_fba.hh
 *  --------------------------------------------------------------------------
 *  THE TWO FLUX BALANCE ANALYSIS PROCESSORS
 *
 *      runFBA_glpk3D     -- solves one linear program per voxel through GLPK
 *      runFBA_cobrapy3D  -- the same, through COBRApy in an embedded Python
 *
 *  Both are optional.  Each compiles to nothing unless the matching build
 *  switch is on, and does nothing unless the matching XML switch is on.
 *
 *  --------------------------------------------------------------------------
 *  WHAT FLUX BALANCE ANALYSIS DOES HERE, IN PLAIN TERMS
 *
 *  A metabolic model is a list of every chemical reaction an organism can run.
 *  At each voxel, at each time step, we ask it: "given how much food is around
 *  you right here, what is the fastest you could possibly grow, and what would
 *  you eat and excrete to do it?"  That question is a linear program, and GLPK
 *  or COBRApy answers it.  The answer -- a growth rate and a set of exchange
 *  fluxes -- becomes this voxel's reaction term for this time step.
 *
 *  The local food supply enters as a bound on the uptake reaction, through a
 *  Michaelis-Menten (Monod) expression:
 *
 *      v_lower  =  -Vmax * C / (Kc + C)          mmol / gDW / h
 *
 *  where C is the local concentration, Kc the half-saturation constant and Vmax
 *  the maximum uptake flux.  The minus sign is the export convention: an
 *  exchange reaction is written as "metabolite ->", so eating is negative.
 *
 *  If a substrate has no <maximum_uptake_flux>, the bound instead becomes the
 *  flux that would exactly exhaust the local supply in one time step -- see
 *  MetabolicConfig::exhaustionFlux().  That makes the organism supply-limited
 *  rather than enzyme-limited.
 *
 *  --------------------------------------------------------------------------
 *  THE FEASIBILITY REPAIR LOOP
 *
 *  Every microbe solves its own LP in ignorance of the others.  Two organisms
 *  eating the same substrate can therefore each be told "you may have all of
 *  it", and between them remove more than exists.  When that happens we detect
 *  it, pool the competing organisms' biomass, re-derive the bound from the
 *  SHARED supply, and solve again.  This is real competition, not a clip after
 *  the fact, and it is the part of CompLaB 2D most worth keeping.
 *
 *  --------------------------------------------------------------------------
 *  LATTICE VECTOR LAYOUT  (this IS the interface -- there is no type to check
 *  it for you, so read it before you change the caller)
 *
 *      [0            .. subsNum-1]        substrate concentrations   C
 *      [subsNum      .. subsNum+bioNum-1] biomass of THIS solver's microbes
 *      [dCloc        .. dCloc+subsNum-1]  substrate increments       dC
 *      [dBloc        .. dBloc+bioNum-1]   biomass increments         dB
 *      [maskLloc]                         mask lattice
 *
 *      dCloc    = subsNum + bioNum
 *      dBloc    = 2*subsNum + bioNum
 *      maskLloc = 2*(subsNum + bioNum)      total length 2*(subsNum+bioNum)+1
 *
 *  `globalId[k]` gives the CompLaB.xml microbe number of the k-th biomass
 *  lattice, so per-microbe parameters are never indexed by the compacted
 *  position.  CompLaB 2D got this wrong in the surrogate path and it is the
 *  reason multi-organism surrogate runs silently misbehaved.
 * ============================================================================
 */

#ifndef COMPLAB3D_PROCESSORS_FBA_HH
#define COMPLAB3D_PROCESSORS_FBA_HH

#include "complab3d_metabolic.hh"

/* By default the FBA processors run on the bulk AND the envelope, exactly like
 * run_kinetics does, which is guaranteed consistent but solves a linear program
 * for every envelope cell as well -- pure waste, since those values are
 * overwritten by communication straight afterwards.
 *
 * Building with -DCOMPLAB3D_FBA_BULK_ONLY restricts them to the bulk and lets
 * Palabos communicate the increments, which is typically a 20-40% saving at
 * realistic block sizes.  It is opt-in rather than default because it changes
 * a correctness-critical setting and deserves to be verified against a
 * bulkAndEnvelope run on your own case first. */
#ifdef COMPLAB3D_FBA_BULK_ONLY
  #define COMPLAB3D_FBA_DOMAIN BlockDomain::bulk
#else
  #define COMPLAB3D_FBA_DOMAIN BlockDomain::bulkAndEnvelope
#endif


/* ============================================================================
 *  runFBA_glpk3D
 * ============================================================================
 */
#ifdef COMPLAB_ENABLE_GLPK

template<typename T, template<typename U> class Descriptor>
class runFBA_glpk3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    runFBA_glpk3D (plint nx_, plint subsNum_, plint bioNum_, T dt_,
                   const std::vector<plint> &globalId_,
                   const MetabolicConfig    *cfg_,
                   const std::vector< std::vector<T> > &vec2_Kc_,
                   const std::vector<T>                &vec1_mu_,
                   plint solid_, plint bb_)
        : nx(nx_), subsNum(subsNum_), bioNum(bioNum_), dt(dt_),
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

        /* Block-relative displacements depend only on the lattices, not on the
         * voxel, so they are computed once here.  The pre-existing 3D processors
         * rebuild this vector inside the innermost z-loop; that is pure overhead
         * and is avoided in the new code. */
        std::vector<Dot3D> off;
        off.reserve(maskLloc);
        for (plint iT = 0; iT < maskLloc; ++iT)
            off.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT]));

        /* scratch, reused for every voxel so nothing is allocated in the hot loop */
        std::vector<T>     conc (subsNum, T());   // free concentration
        std::vector<T>     avail(subsNum, T());   // free or total, per the XML basis
        std::vector<T>     dC   (subsNum, T());
        std::vector<T>     bmass(bioNum,  T());
        std::vector<plint> bLoc;                  // compacted -> local biomass index
        std::vector<T>     dB   (bioNum,  T());

        std::vector<int>    location;
        std::vector<double> lpLb, lpUb, xmin, lambda, redcosts;
        std::vector< std::vector<T> > flux(bioNum, std::vector<T>(subsNum, T()));
        std::vector<T>      grate(bioNum, T());
        std::vector< std::vector<bool> > shared(bioNum, std::vector<bool>(subsNum, false));

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            const plint absX = iX + absoluteOffset.x;
            if (absX <= 0 || absX >= nx-1) continue;           // inlet / outlet planes

            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {

                    const plint mask = util::roundToInt(
                        lattices[maskLloc]->get(iX+maskOffset.x, iY+maskOffset.y, iZ+maskOffset.z).computeDensity());
                    if (mask == solid || mask == bb) continue;

                    /* ---- who is here? ---------------------------------- */
                    bLoc.clear();
                    for (plint iB = 0; iB < bioNum; ++iB) {
                        const T b = lattices[subsNum+iB]->get(iX+off[subsNum+iB].x,
                                                              iY+off[subsNum+iB].y,
                                                              iZ+off[subsNum+iB].z).computeDensity();
                        bmass[iB] = (b > thrd) ? b : T();
                        if (bmass[iB] > thrd) bLoc.push_back(iB);
                    }
                    if (bLoc.empty()) continue;

                    /* ---- what is there to eat? ------------------------- */
                    for (plint iS = 0; iS < subsNum; ++iS) {
                        const T c = lattices[iS]->get(iX+off[iS].x, iY+off[iS].y, iZ+off[iS].z).computeDensity();
                        conc[iS] = (c > thrd) ? c : T();
                        dC[iS]   = T();
                    }
                    for (plint iS = 0; iS < subsNum; ++iS)
                        avail[iS] = cfg->useTotals ? cfg->totals.total(conc, iS) : conc[iS];

                    for (plint iB = 0; iB < bioNum; ++iB) { grate[iB] = T(); dB[iB] = T(); }

                    /* ==================================================== *
                     *  solve, check the substrate budget, repair, re-solve  *
                     * ==================================================== */
                    for (plint iB = 0; iB < bioNum; ++iB)
                        std::fill(shared[iB].begin(), shared[iB].end(), false);
                    bool needRepair = true;
                    int  pass       = 0;

                    while (needRepair && pass < 3) {
                        needRepair = false;

                        for (size_t k = 0; k < bLoc.size(); ++k) {
                            const plint iB = bLoc[k];
                            const plint gM = globalId[iB];

                            location.clear(); lpLb.clear(); lpUb.clear();

                            for (plint iS = 0; iS < subsNum; ++iS) {
                                const plint col = cfg->subsLoc[gM][iS];
                                if (col < 0) continue;

                                T lo;
                                if (shared[iB][iS]) {
                                    /* repair pass: share the supply with every other
                                     * organism that wanted this substrate */
                                    T pooled = T();
                                    for (size_t k2 = 0; k2 < bLoc.size(); ++k2)
                                        if (shared[bLoc[k2]][iS]) pooled += bmass[bLoc[k2]];
                                    lo = cfg->exhaustionFlux(avail[iS], pooled, gM, dt);
                                } else {
                                    const T mm = (vec2_Kc[gM][iS] + avail[iS] > T())
                                               ? avail[iS] / (vec2_Kc[gM][iS] + avail[iS]) : T();
                                    /* [FIX-3D] Two defects here at once.
                                     *  (a) The gate tested maxRelease -- the release
                                     *      CEILING -- to decide whether an uptake Vmax
                                     *      had been supplied.  Since maxRelease defaults
                                     *      to +1e30, giving <substrate_lower_bounds>
                                     *      without <substrate_upper_bounds> silently
                                     *      dropped to the supply-limited branch, and the
                                     *      Monod attenuation was bypassed entirely.
                                     *  (b) It used maxUptake (the hard floor) as the
                                     *      Vmax.  Vmax now comes from its own field, fed
                                     *      by <maximum_uptake_flux>.
                                     * The cobrapy and surrogate paths use the identical
                                     * predicate, so the three agree. */
                                    if (cfg->vmax[gM][iS] > thrd) {
                                        lo = -cfg->vmax[gM][iS] * mm;          // enzyme-limited
                                    } else {
                                        lo = cfg->exhaustionFlux(avail[iS], bmass[iB], gM, dt) * mm;
                                    }
                                }

                                /* never ask for more than the model itself allows */
                                if (lo < cfg->maxUptake[gM][iS]) lo = cfg->maxUptake[gM][iS];

                                /* [FIX-3D] CompLaB 2D complab_processors.hh:130 wrote
                                 *     if (fixLB && c2f < conc) c2f = -conc;
                                 * but c2f is negative and conc is positive, so the test was
                                 * true in essentially every cell and fix_lower_bounds silently
                                 * overrode the Michaelis-Menten bound everywhere.  The COBRApy
                                 * twin at line 366 had the correct comparison.  Fixed here. */
                                if (cfg->fix_lower_bounds[iS] && lo < -avail[iS]) lo = -avail[iS];

                                T hi = cfg->maxRelease[gM][iS];
                                if (hi > (T) 1e29) hi = (T) 1e6;
                                if (cfg->equate_bounds[gM]) hi = lo + (T) thrd;
                                if (hi < lo) hi = lo + (T) thrd;

                                location.push_back((int) col);
                                lpLb.push_back((double) lo);
                                lpUb.push_back((double) hi);
                            }

                            double fmin = 0.0, tsec = 0.0, tmem = 0.0;
                            int    status = 0;
                            /* run_glpk sizes these itself now; assigning here as well
                             * keeps the intent visible at the call site. */
                            xmin.assign((size_t) cfg->vec_nrxns[gM], 0.0);
                            lambda.assign((size_t) cfg->vec_nmets[gM], 0.0);
                            redcosts.assign((size_t) cfg->vec_nrxns[gM], 0.0);

                            const int glpkerr = run_glpk(
                                cfg->vec_lp[gM], cfg->method[gM], cfg->isMIP[gM],
                                cfg->vec_nmets[gM], cfg->vec_nrxns[gM], cfg->lpsolver,
                                cfg->sParam[gM], cfg->iParam[gM],
                                xmin, fmin, status, location, lpLb, lpUb,
                                lambda, redcosts, tsec, tmem);

                            /* [FIX-3D] `glpkerr >= 0` accepted almost every failure, since
                             * run_glpk returns positive GLPK codes on failure and only the
                             * setup errors are negative.  Success is exactly 0.  The
                             * legacy `status == 180` (LPX_OPT) is unreachable with any
                             * modern GLPK and has been dropped. */
                            const bool ok = (glpkerr == 0) && (status == GLP_OPT);

                            if (!ok) {
                                /* [FIX-3D] CompLaB 2D printed a decoded GLPK error and set
                                 * everything to zero, which turns an infeasible LP into silent
                                 * decay.  We do the same numerically -- there is no better local
                                 * answer -- but we count it so the run summary can report how
                                 * often it happened instead of hiding it. */
                                grate[iB] = T();
                                for (plint iS = 0; iS < subsNum; ++iS) flux[iB][iS] = T();
                                noteInfeasible(gM);
                            } else {
                                grate[iB] = (T) fmin;
                                plint c = 0;
                                for (plint iS = 0; iS < subsNum; ++iS) {
                                    const plint col = cfg->subsLoc[gM][iS];
                                    flux[iB][iS] = (col >= 0 && col < (plint) xmin.size())
                                                 ? (T) xmin[col] : T();
                                    if (col >= 0) ++c;
                                }
                                (void) c;
                            }
                        }

                        /* ---- budget check ------------------------------ */
                        for (plint iS = 0; iS < subsNum; ++iS) {
                            if (cfg->fix_concentration[iS] || cfg->fix_lower_bounds[iS]) continue;
                            T draw = T();
                            for (size_t k = 0; k < bLoc.size(); ++k) {
                                const plint iB = bLoc[k];
                                draw += cfg->fluxToDeltaC(flux[iB][iS], bmass[iB], globalId[iB], dt);
                            }
                            /* [FIX-3D] CompLaB 2D:204 only repaired when conc > thrd, so a
                             * substrate already floored to zero accepted an unlimited negative
                             * draw and the lattice went negative.  The >= 0 test below catches
                             * that case too. */
                            if (draw + avail[iS] < -thrd) {
                                bool any = false;
                                for (size_t k = 0; k < bLoc.size(); ++k) {
                                    const plint iB = bLoc[k];
                                    if (flux[iB][iS] < -thrd) { shared[iB][iS] = true; any = true; }
                                }
                                if (any) needRepair = true;
                            }
                        }
                        ++pass;
                    }

                    /* ==================================================== *
                     *  turn fluxes into increments                          *
                     * ==================================================== */
                    for (plint iS = 0; iS < subsNum; ++iS) {
                        if (cfg->fix_concentration[iS]) continue;
                        T draw = T();
                        for (size_t k = 0; k < bLoc.size(); ++k) {
                            const plint iB = bLoc[k];
                            draw += cfg->fluxToDeltaC(flux[iB][iS], bmass[iB], globalId[iB], dt);
                        }
                        if (std::fabs(draw) <= thrd) continue;

                        /* [FIX-3D] the last line of defence.  CompLaB 2D had no positivity
                         * clamp anywhere in the increment path -- the commented-out in-place
                         * kinetics code had one, the live delta code did not -- so
                         * concentrations could and did go negative.  After three repair passes
                         * we clamp rather than give up, and never below zero. */
                        if (draw + avail[iS] < T()) draw = -avail[iS];

                        if (cfg->useTotals) cfg->totals.drawDown(conc, iS, -draw, dC);
                        else                dC[iS] += draw;
                    }

                    for (size_t k = 0; k < bLoc.size(); ++k) {
                        const plint iB = bLoc[k];
                        const plint gM = globalId[iB];
                        /* No growth predicted -> the organism decays instead.
                         * [FIX-3D] decay is per SECOND, growth is per hour; they cannot
                         * share one converter.  See decayToDeltaB(). */
                        const T mu = grate[iB];
                        dB[iB] = (std::fabs(mu) < thrd)
                               ? cfg->decayToDeltaB(vec1_mu[gM], bmass[iB], dt)
                               : cfg->growthToDeltaB(mu, bmass[iB], dt);
                        if (dB[iB] + bmass[iB] < T()) dB[iB] = -bmass[iB];   // biomass never negative
                    }

                    /* ==================================================== *
                     *  deposit into the increment lattices                  *
                     * ==================================================== */
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

    virtual BlockDomain::DomainT appliesTo() const { return COMPLAB3D_FBA_DOMAIN; }

    virtual runFBA_glpk3D<T,Descriptor>* clone() const {
        return new runFBA_glpk3D<T,Descriptor>(*this);
    }

    void getTypeOfModification (std::vector<modif::ModifT> &modified) const {
        for (size_t i = 0; i < modified.size(); ++i) modified[i] = modif::nothing;
        for (plint iT = dCloc; iT < maskLloc; ++iT) modified[iT] = modif::staticVariables;
    }

    /* [FIX-3D] This used to be a bare `static long long infeasibleCount` with a
     * comment claiming the run summary reported it.  Nothing read it, it was not
     * MPI-reduced, and it counted solves rather than voxels -- dead state whose
     * comment lied.  An infeasible LP is now actually visible: the first few are
     * reported, then the message is suppressed so it cannot flood the log. */
    static long long infeasibleCount;
    static void noteInfeasible (plint gM) {
        ++infeasibleCount;
        if (infeasibleCount <= 5) {
            pcout << "  [FBA] microbe" << gM << ": the linear program did not solve at some voxel; "
                  << "that voxel gets zero growth and zero flux this step.\n";
            if (infeasibleCount == 5)
                pcout << "  [FBA] further infeasible-solve messages suppressed.\n";
        }
    }

private:
    plint nx, subsNum, bioNum;
    T dt;
    std::vector<plint> globalId;
    const MetabolicConfig *cfg;
    std::vector< std::vector<T> > vec2_Kc;
    std::vector<T> vec1_mu;
    plint solid, bb, dCloc, dBloc, maskLloc;
};

template<typename T, template<typename U> class Descriptor>
long long runFBA_glpk3D<T,Descriptor>::infeasibleCount = 0;

#endif  // COMPLAB_ENABLE_GLPK


/* ============================================================================
 *  runFBA_cobrapy3D
 *
 *  Same physics, different back end.  The bounds are marshalled into one flat
 *  argument vector, handed to complab3d_cobrapy.run_FBA, and the answer comes
 *  back as a growth rate plus a PACKED flux vector per microbe -- packed over
 *  that microbe's own exchange list, NOT over substrates.  Confusing the two is
 *  what broke the 2D version; expandFluxToSubstrates() is the only correct way
 *  to unpack it, and it is called below.
 * ============================================================================
 */
#ifdef COMPLAB_ENABLE_COBRAPY

template<typename T, template<typename U> class Descriptor>
class runFBA_cobrapy3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    runFBA_cobrapy3D (plint nx_, plint subsNum_, plint bioNum_, T dt_,
                      const std::vector<plint> &globalId_,
                      MetabolicConfig *cfg_, char *pyFileName_,
                      const std::vector<plint> &modelSlot_,
                      const std::vector< std::vector<T> > &vec2_Kc_,
                      const std::vector<T>                &vec1_mu_,
                      plint solid_, plint bb_)
        : nx(nx_), subsNum(subsNum_), bioNum(bioNum_), dt(dt_),
          globalId(globalId_), cfg(cfg_), pyFileName(pyFileName_), modelSlot(modelSlot_),
          vec2_Kc(vec2_Kc_), vec1_mu(vec1_mu_), solid(solid_), bb(bb_),
          dCloc(subsNum_ + bioNum_),
          dBloc(2*subsNum_ + bioNum_),
          maskLloc(2*(subsNum_ + bioNum_))
    {}

    virtual void process (Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices)
    {
        Dot3D absoluteOffset = lattices[0]->getLocation();
        Dot3D maskOffset     = computeRelativeDisplacement(*lattices[0], *lattices[maskLloc]);

        /* Block-relative displacements depend only on the lattices, not on the
         * voxel, so they are computed once here.  The pre-existing 3D processors
         * rebuild this vector inside the innermost z-loop; that is pure overhead
         * and is avoided in the new code. */
        std::vector<Dot3D> off;
        off.reserve(maskLloc);
        for (plint iT = 0; iT < maskLloc; ++iT)
            off.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT]));

        std::vector<T> conc(subsNum, T()), avail(subsNum, T()), dC(subsNum, T());
        std::vector<T> bmass(bioNum, T()), dB(bioNum, T());
        std::vector< std::vector<T> > flux(bioNum, std::vector<T>(subsNum, T()));
        std::vector<plint> bLoc;

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
                        if (bmass[iB] > thrd) bLoc.push_back(iB);
                    }
                    if (bLoc.empty()) continue;

                    for (plint iS = 0; iS < subsNum; ++iS) {
                        const T c = lattices[iS]->get(iX+off[iS].x, iY+off[iS].y, iZ+off[iS].z).computeDensity();
                        conc[iS] = (c > thrd) ? c : T();
                        dC[iS]   = T();
                    }
                    for (plint iS = 0; iS < subsNum; ++iS)
                        avail[iS] = cfg->useTotals ? cfg->totals.total(conc, iS) : conc[iS];
                    for (plint iB = 0; iB < bioNum; ++iB) dB[iB] = T();

                    /* ---- marshal the argument vector -------------------- *
                     * layout, matching complab3d_cobrapy.run_FBA exactly:
                     *   [0]                 M            number of microbes present
                     *   [1 .. M]            n_i          exchanges used by microbe i
                     *   [1+M     .. ]       loc          concatenated, n_i each
                     *   [1+M+N   .. ]       lb           concatenated
                     *   [1+M+2N  .. ]       ub           concatenated
                     * (the model handles are appended by optimize_cobrapy)     */
                    const plint M = (plint) bLoc.size();
                    std::vector< std::vector<int> > locRow(M);
                    std::vector< std::vector<double> > lbRow(M), ubRow(M);

                    for (plint k = 0; k < M; ++k) {
                        const plint iB = bLoc[k];
                        const plint gM = globalId[iB];
                        for (plint iS = 0; iS < subsNum; ++iS) {
                            const plint col = cfg->subsLoc[gM][iS];
                            if (col < 0) continue;

                            const T mm = (vec2_Kc[gM][iS] + avail[iS] > T())
                                       ? avail[iS] / (vec2_Kc[gM][iS] + avail[iS]) : T();
                            T lo;
                            if (cfg->vmax[gM][iS] > thrd) {
                                lo = -cfg->vmax[gM][iS] * mm;                  // enzyme-limited
                            } else {
                                lo = cfg->exhaustionFlux(avail[iS], bmass[iB], gM, dt) * mm;
                            }
                            if (lo < cfg->maxUptake[gM][iS]) lo = cfg->maxUptake[gM][iS];
                            if (cfg->fix_lower_bounds[iS] && lo < -avail[iS]) lo = -avail[iS];

                            /* [FIX-3D] CompLaB 2D:377 read
                             *     vec2_maxRelease[bmLattice[bloc]][iS]
                             * which indexes the compacted list with an already-compacted
                             * index -- a double indirection that silently fetched another
                             * microbe's release ceiling. */
                            T hi = cfg->maxRelease[gM][iS];
                            if (hi > (T) 1e29) hi = (T) 1e6;
                            if (cfg->equate_bounds[gM]) hi = lo + (T) thrd;
                            if (hi < lo) hi = lo + (T) thrd;

                            /* [FIX-3D] CompLaB 2D:357 declared the per-microbe bound vector
                             * OUTSIDE the microbe loop and never cleared it, so microbe i
                             * inherited every earlier microbe's bounds.  These rows are local
                             * to the iteration, so the bug cannot recur. */
                            locRow[k].push_back((int) col);
                            lbRow [k].push_back((double) lo);
                            ubRow [k].push_back((double) hi);
                        }
                    }

                    std::vector<double> args;
                    args.push_back((double) M);
                    for (plint k = 0; k < M; ++k) args.push_back((double) locRow[k].size());
                    for (plint k = 0; k < M; ++k) for (size_t j = 0; j < locRow[k].size(); ++j) args.push_back((double) locRow[k][j]);
                    for (plint k = 0; k < M; ++k) for (size_t j = 0; j < lbRow [k].size(); ++j) args.push_back(lbRow[k][j]);
                    for (plint k = 0; k < M; ++k) for (size_t j = 0; j < ubRow [k].size(); ++j) args.push_back(ubRow[k][j]);

                    /* only the models of the microbes actually present */
                    std::vector<PyObject*> here;
                    for (plint k = 0; k < M; ++k) here.push_back(cfg->vec_model[modelSlot[globalId[bLoc[k]]]]);

                    std::vector<double> grate;
                    std::vector< std::vector<double> > packed;
                    std::vector<int> status;
                    const int erck = optimize_cobrapy(pyFileName, args, here, grate, packed, status);

                    /* Unpack every microbe's answer FIRST, then apply the substrate
                     * budget across all of them together.
                     *
                     * [FIX-3D] The first version clamped inside this loop, per microbe,
                     * against an `avail` that was never decremented.  With M microbes
                     * each individually able to exhaust the voxel, each was clamped to
                     * -avail and the total came to -M*avail: two cobrapy organisms in
                     * one voxel drove the substrate lattice negative on the first step.
                     * The GLPK path already summed first and clamped once; this now
                     * matches it. */
                    for (plint k = 0; k < M; ++k) {
                        const plint iB = bLoc[k];
                        const plint gM = globalId[iB];

                        for (plint iS = 0; iS < subsNum; ++iS) flux[iB][iS] = T();
                        T mu = T();

                        /* [FIX-3D] the flux vector coming back is PACKED over locRow[k];
                         * CompLaB 2D indexed it by substrate and silently mis-assigned every
                         * flux whenever a microbe did not use every substrate. */
                        if (erck == 0 && k < (plint) status.size() && status[k] != 0) {
                            mu = (T) grate[k];
                            std::vector<int> subsLocRow(subsNum);
                            for (plint iS = 0; iS < subsNum; ++iS) subsLocRow[iS] = (int) cfg->subsLoc[gM][iS];
                            std::vector<T> full(subsNum, T());
                            expandFluxToSubstrates(packed[k], subsLocRow, subsNum, full);
                            for (plint iS = 0; iS < subsNum; ++iS) flux[iB][iS] = full[iS];
                        } else {
                            noteInfeasible(gM);
                        }

                        dB[iB] = (std::fabs(mu) < thrd)
                               ? cfg->decayToDeltaB(vec1_mu[gM], bmass[iB], dt)   // [FIX-3D] per second
                               : cfg->growthToDeltaB(mu, bmass[iB], dt);
                        if (dB[iB] + bmass[iB] < T()) dB[iB] = -bmass[iB];
                    }

                    /* joint substrate budget: sum over every organism present, clamp once */
                    for (plint iS = 0; iS < subsNum; ++iS) {
                        if (cfg->fix_concentration[iS]) continue;
                        T draw = T();
                        for (plint k = 0; k < M; ++k) {
                            const plint iB = bLoc[k];
                            draw += cfg->fluxToDeltaC(flux[iB][iS], bmass[iB], globalId[iB], dt);
                        }
                        if (std::fabs(draw) <= thrd) continue;
                        /* [FIX-3D] Clamp against what is ACTUALLY left, not the pristine
                         * `avail`.  With fba_concentration_basis = total, drawDown spreads a
                         * component's draw over its carrier species -- and a carrier that is
                         * itself a transported substrate was then debited a second time by
                         * its own loop iteration, which knew nothing about the first.  A
                         * substrate could go negative even though every individual clamp was
                         * satisfied.  headroom sees the running dC. */
                        const T headroom = avail[iS] + (cfg->useTotals ? T() : dC[iS]);
                        if (draw + headroom < T()) draw = -headroom;
                        if (draw > T()) draw = T();          // consumption only, never a spurious source
                        if (cfg->useTotals) cfg->totals.drawDown(conc, iS, -draw, dC);
                        else                dC[iS] += draw;
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

    virtual BlockDomain::DomainT appliesTo() const { return COMPLAB3D_FBA_DOMAIN; }

    virtual runFBA_cobrapy3D<T,Descriptor>* clone() const {
        return new runFBA_cobrapy3D<T,Descriptor>(*this);
    }

    void getTypeOfModification (std::vector<modif::ModifT> &modified) const {
        for (size_t i = 0; i < modified.size(); ++i) modified[i] = modif::nothing;
        for (plint iT = dCloc; iT < maskLloc; ++iT) modified[iT] = modif::staticVariables;
    }

    static long long infeasibleCount;
    static void noteInfeasible (plint gM) {
        ++infeasibleCount;
        if (infeasibleCount <= 5) {
            pcout << "  [FBA] microbe" << gM << ": COBRApy did not return a usable solution at some "
                  << "voxel; that voxel gets zero growth and zero flux this step.\n";
            if (infeasibleCount == 5)
                pcout << "  [FBA] further infeasible-solve messages suppressed.\n";
        }
    }

private:
    plint nx, subsNum, bioNum;
    T dt;
    std::vector<plint> globalId;
    MetabolicConfig *cfg;
    char *pyFileName;
    /* cfg->vec_model is packed over COBRApy microbes only, in ascending global
     * order, because that is the order prep_cobrapy() was handed.  modelSlot is
     * indexed by GLOBAL microbe id and gives the position in that packed list;
     * it is built once in complab.cpp rather than recomputed per voxel. */
    std::vector<plint> modelSlot;
    std::vector< std::vector<T> > vec2_Kc;
    std::vector<T> vec1_mu;
    plint solid, bb, dCloc, dBloc, maskLloc;
};

template<typename T, template<typename U> class Descriptor>
long long runFBA_cobrapy3D<T,Descriptor>::infeasibleCount = 0;

#endif  // COMPLAB_ENABLE_COBRAPY

#endif  // COMPLAB3D_PROCESSORS_FBA_HH
