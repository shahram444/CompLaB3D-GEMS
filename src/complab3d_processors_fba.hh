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
#include "complab3d_thermo.hh"      /* the F_T gate; a no-op unless <thermodynamics> is on */

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

        std::vector<double> tconc;                 /* the gate's view of the same chemistry */
        const bool gated = complab_thermo::enabled();

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

                            /* ------------------------------------------------------------ *
                             *  THE SOLVE.
                             *
                             *  Three arrangements share this one call site, and which one runs
                             *  is decided by configuration alone:
                             *
                             *    plain          no <multi_step>, no <cybernetic>.  One LP.
                             *    multi-step     <multi_step> only.  run_glpk_lex walks the
                             *                   lexicographic chain and returns the growth rate
                             *                   from the first stage and the fluxes from the
                             *                   last.  See complab3d_lexicographic.hh.
                             *    cybernetic     <cybernetic>.  The block below solves the whole
                             *                   chain once per carbon source and blends the
                             *                   answers; each of those solves comes back through
                             *                   here.  See complab3d_cybernetic.hh.
                             *
                             *  Keeping them on one call site is deliberate: the bound
                             *  construction above, the budget repair below and the infeasibility
                             *  accounting are shared, so the three arrangements cannot drift
                             *  apart in the parts that are not about the objective.
                             * ------------------------------------------------------------ */
                            /* ---- how many growth options are there, and how is each weighted?
                             *
                             *  Without <cybernetic> there is exactly one option, weighted 1, and
                             *  everything below collapses to a single solve with the bounds built
                             *  above -- the same call, the same arguments, the same answer.
                             *
                             *  With <cybernetic> there is one option per carbon source.  Each is
                             *  solved with the OTHER carbon sources shut, so the model has to
                             *  make a living on that source alone, and the results are averaged
                             *  with the cybernetic weights.  Shutting a source means forbidding
                             *  its UPTAKE, not its release: an organism growing on lactate must
                             *  still be free to excrete pyruvate, and that excretion is the whole
                             *  reason the switch has something to switch to later.           */
                            const complab_cyb::CyberneticPlan &cyb = cfg->cybPlan[gM];
                            const bool cybOn = cyb.active();

                            std::vector<double> uw;
                            if (cybOn && !complab_cyb::weights(cyb, avail, uw)) {
                                /* no carbon at all in this voxel: nothing to solve, nothing to
                                 * blend.  Zero rates are the right answer, not an infeasibility. */
                                grate[iB] = T();
                                for (plint iS = 0; iS < subsNum; ++iS) flux[iB][iS] = T();
                                continue;
                            }

                            const size_t nOpt = cybOn ? cyb.size() : 1;
                            T   accGrate = T();
                            bool anyOk   = false, anyFail = false;
                            std::vector<T> accFlux(subsNum, T());

                            /* the bounds as built above, kept intact so each option starts from
                             * the same place rather than from the previous option's closures */
                            const std::vector<double> baseLb = lpLb;

                            for (size_t opt = 0; opt < nOpt; ++opt) {

                                const double w = cybOn ? uw[opt] : 1.0;
                                if (cybOn && w <= cyb.weightFloor) continue;   /* see the header:
                                                                                * a near-zero
                                                                                * weight is not
                                                                                * worth an LP */

                                lpLb = baseLb;
                                if (cybOn) {
                                    /* shut every other carbon source's uptake for this option */
                                    const std::vector<int> &excl = cyb.sources[opt].exclusive;
                                    for (size_t e = 0; e < excl.size(); ++e) {
                                        const plint col = cfg->subsLoc[gM][excl[e]];
                                        if (col < 0) continue;
                                        for (size_t sl = 0; sl < location.size(); ++sl)
                                            if (location[sl] == (int) col) {
                                                if (lpLb[sl] < 0.0) lpLb[sl] = 0.0;
                                                break;
                                            }
                                    }
                                }

                                double fOpt = 0.0, tOpt = 0.0, mOpt = 0.0;
                                int    sOpt = 0;

                                const int glpkerr = run_glpk_lex(
                                    cfg->vec_lp[gM], cfg->method[gM], cfg->isMIP[gM],
                                    cfg->vec_nmets[gM], cfg->vec_nrxns[gM], cfg->lpsolver,
                                    cfg->sParam[gM], cfg->iParam[gM],
                                    xmin, fOpt, sOpt, location, lpLb, lpUb,
                                    lambda, redcosts, tOpt, mOpt,
                                    cfg->lexPlan[gM], cfg->vec_c[gM],
                                    cfg->vec_lb[gM], cfg->vec_ub[gM],
                                    (int) cfg->vec_objLoc[gM]);
                                tsec += tOpt; tmem = mOpt; status = sOpt; fmin = fOpt;

                                /* [FIX-3D] `glpkerr >= 0` accepted almost every failure, since
                                 * run_glpk returns positive GLPK codes on failure and only the
                                 * setup errors are negative.  Success is exactly 0.  The
                                 * legacy `status == 180` (LPX_OPT) is unreachable with any
                                 * modern GLPK and has been dropped. */
                                if (glpkerr != 0 || sOpt != GLP_OPT) { anyFail = true; continue; }

                                anyOk = true;
                                accGrate += (T) (w * fOpt);
                                for (plint iS = 0; iS < subsNum; ++iS) {
                                    const plint col = cfg->subsLoc[gM][iS];
                                    if (col >= 0 && col < (plint) xmin.size())
                                        accFlux[iS] += (T) (w * xmin[col]);
                                }
                            }

                            lpLb = baseLb;      /* leave the vectors as the repair pass expects */

                            if (!anyOk) {
                                /* [FIX-3D] CompLaB 2D printed a decoded GLPK error and set
                                 * everything to zero, which turns an infeasible LP into silent
                                 * decay.  We do the same numerically -- there is no better local
                                 * answer -- but we count it so the run summary can report how
                                 * often it happened instead of hiding it. */
                                grate[iB] = T();
                                for (plint iS = 0; iS < subsNum; ++iS) flux[iB][iS] = T();
                                noteInfeasible(gM);
                            } else {
                                /* One option out of several failing is not a failed voxel: the
                                 * organism simply cannot live on that source here, which is a
                                 * physical statement, not a solver problem.  It is still counted,
                                 * because a run in which it happens everywhere is telling you the
                                 * medium is wrong. */
                                if (anyFail) noteInfeasible(gM);
                                grate[iB] = accGrate;
                                for (plint iS = 0; iS < subsNum; ++iS) flux[iB][iS] = accFlux[iS];
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

                    /* ---- the thermodynamic gate -----------------------------------------
                     * Applied to the SOLVED fluxes, not to the uptake bounds.  A linear program
                     * asked for a smaller bound would redistribute its whole flux distribution
                     * and could hand back a different byproduct pattern, which is a modelling
                     * change; scaling the solution keeps the metabolic answer and asks only what
                     * fraction of it the local energy balance permits.  The budget above was
                     * computed on the ungated fluxes, so scaling here can only reduce the draw
                     * and cannot break it.  A no-op when <thermodynamics> is off. */
                    if (gated) {
                        complab_thermo::fillConc(avail, (int) subsNum, tconc);
                        for (size_t k = 0; k < bLoc.size(); ++k) {
                            const plint iB = bLoc[k];
                            const double ft = complab_thermo::gateFor((int) globalId[iB], tconc);
                            if (ft >= 1.0) continue;
                            grate[iB] *= (T) ft;
                            for (plint iS = 0; iS < subsNum; ++iS) flux[iB][iS] *= (T) ft;
                        }
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
        /* One row per organism, one column per substrate: which substrates this organism has been
         * told to share on the repair pass. The GLPK path has carried this since the port; the
         * COBRApy path now runs the same loop and needs the same marking. */
        std::vector< std::vector<bool> > shared(bioNum, std::vector<bool>(subsNum, false));

        std::vector<double> tconc;                 /* the gate's view of the same chemistry */
        const bool gated = complab_thermo::enabled();

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

                    /* The gate reads the same vector the bounds were built from. */
                    if (gated) complab_thermo::fillConc(avail, (int) subsNum, tconc);
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

                    /* [FIX] THE REPAIR LOOP, which this path did not have.
                     *
                     * Every organism present is given an uptake bound built from the FULL local
                     * supply, because each is solved without knowing the others are there. Two
                     * organisms competing for one substrate therefore each claim all of it, and
                     * the joint budget below then scales both back. That keeps the concentration
                     * positive -- nothing here could ever go negative -- but it splits the
                     * substrate by clamping rather than by re-optimising, and a linear program
                     * given a smaller bound does not return a scaled version of its previous
                     * answer: it can pick a different flux distribution altogether.
                     *
                     * The GLPK path has solved this since the port: mark the substrates that were
                     * over-drawn, rebuild those bounds against the POOLED biomass, and solve
                     * again, up to three passes. This is that loop, over the batched call.
                     *
                     * The cost is one extra Python round trip in the voxels where it fires, and
                     * none anywhere else. With a single organism -- which is every shipped COBRApy
                     * case -- it never fires at all, because one organism cannot over-draw a
                     * budget built from what it can actually reach. */
                    for (plint k = 0; k < M; ++k)
                        std::fill(shared[bLoc[k]].begin(), shared[bLoc[k]].end(), false);
                    bool needRepair = true;
                    int  pass = 0;
                    std::vector<double> grate;
                    std::vector< std::vector<double> > packed;
                    std::vector<int> status;
                    int erck = 0;

                    while (needRepair && pass < 3) {
                    needRepair = false;
                    for (plint k = 0; k < M; ++k) { locRow[k].clear(); lbRow[k].clear(); ubRow[k].clear(); }

                    for (plint k = 0; k < M; ++k) {
                        const plint iB = bLoc[k];
                        const plint gM = globalId[iB];
                        for (plint iS = 0; iS < subsNum; ++iS) {
                            const plint col = cfg->subsLoc[gM][iS];
                            if (col < 0) continue;

                            const T mm = (vec2_Kc[gM][iS] + avail[iS] > T())
                                       ? avail[iS] / (vec2_Kc[gM][iS] + avail[iS]) : T();
                            T lo;
                            if (shared[iB][iS]) {
                                /* repair pass: this substrate was over-drawn, so share the
                                 * supply with every other organism that wanted it. Same rule
                                 * as the GLPK path, same function. */
                                T pooled = T();
                                for (plint k2 = 0; k2 < M; ++k2)
                                    if (shared[bLoc[k2]][iS]) pooled += bmass[bLoc[k2]];
                                lo = cfg->exhaustionFlux(avail[iS], pooled, gM, dt);
                            } else if (cfg->vmax[gM][iS] > thrd) {
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

                    /* only the models of the microbes actually present */
                    std::vector<PyObject*> here;
                    for (plint k = 0; k < M; ++k) here.push_back(cfg->vec_model[modelSlot[globalId[bLoc[k]]]]);

                    /* ---- growth options, and the weight of each --------------------------
                     *
                     *  Mirror of the GLPK path, with one structural difference forced by the
                     *  bridge: COBRApy is called once for ALL microbes present in the voxel, not
                     *  once per microbe, because the round trip through the interpreter costs far
                     *  more than the solve. So the loop is over OPTION INDEX, and every microbe
                     *  contributes its own option of that index -- or nothing, if it has fewer
                     *  options than the microbe with the most. One Python call per option index,
                     *  not one per (microbe, option) pair.
                     *
                     *  Without <cybernetic> anywhere in the voxel this runs exactly once with
                     *  weight 1 and the bounds built above, which is the original single call. */
                    std::vector< std::vector<double> > uwM((size_t) M);
                    size_t nOptMax = 1;
                    for (plint k = 0; k < M; ++k) {
                        const plint gM = globalId[bLoc[k]];
                        const complab_cyb::CyberneticPlan &cyb = cfg->cybPlan[gM];
                        if (!cyb.active()) continue;
                        if (!complab_cyb::weights(cyb, avail, uwM[(size_t) k])) {
                            uwM[(size_t) k].assign(cyb.size(), 0.0);   /* no carbon: contributes
                                                                        * nothing, and must not
                                                                        * fall back to weight 1 */
                        }
                        if (cyb.size() > nOptMax) nOptMax = cyb.size();
                    }

                    std::vector<T>                gAcc((size_t) M, T());
                    std::vector< std::vector<T> > pAcc((size_t) M);
                    std::vector<int>              sAcc((size_t) M, 0);
                    bool anyCall = false;

                    for (size_t opt = 0; opt < nOptMax; ++opt) {

                        /* which microbes take part in this option index, and with what weight */
                        std::vector<double> wRow((size_t) M, 0.0);
                        bool anyThisOpt = false;
                        for (plint k = 0; k < M; ++k) {
                            const plint gM = globalId[bLoc[k]];
                            const complab_cyb::CyberneticPlan &cyb = cfg->cybPlan[gM];
                            if (cyb.active()) {
                                if (opt < uwM[(size_t) k].size() &&
                                    uwM[(size_t) k][opt] > cyb.weightFloor) {
                                    wRow[(size_t) k] = uwM[(size_t) k][opt];
                                    anyThisOpt = true;
                                }
                            } else if (opt == 0) {
                                wRow[(size_t) k] = 1.0;
                                anyThisOpt = true;
                            }
                        }
                        if (!anyThisOpt) continue;

                        /* Bounds for this option: the row built above, with the OTHER carbon
                         * sources' uptake shut for any microbe that is switching. Shutting means
                         * lb -> 0: no uptake, but release is untouched, because an organism
                         * growing on lactate must still be free to excrete the pyruvate that a
                         * later option will consume. */
                        std::vector< std::vector<double> > lbOpt = lbRow;
                        for (plint k = 0; k < M; ++k) {
                            if (wRow[(size_t) k] <= 0.0) continue;
                            const plint gM = globalId[bLoc[k]];
                            const complab_cyb::CyberneticPlan &cyb = cfg->cybPlan[gM];
                            if (!cyb.active() || opt >= cyb.size()) continue;
                            const std::vector<int> &excl = cyb.sources[opt].exclusive;
                            for (size_t e = 0; e < excl.size(); ++e) {
                                const plint col = cfg->subsLoc[gM][excl[e]];
                                if (col < 0) continue;
                                for (size_t sl = 0; sl < locRow[k].size(); ++sl)
                                    if (locRow[k][sl] == (int) col) {
                                        if (lbOpt[k][sl] < 0.0) lbOpt[k][sl] = 0.0;
                                        break;
                                    }
                            }
                        }

                        std::vector<double> args;
                        args.push_back((double) M);
                        for (plint k = 0; k < M; ++k) args.push_back((double) locRow[k].size());
                        for (plint k = 0; k < M; ++k) for (size_t j = 0; j < locRow[k].size(); ++j) args.push_back((double) locRow[k][j]);
                        for (plint k = 0; k < M; ++k) for (size_t j = 0; j < lbOpt [k].size(); ++j) args.push_back(lbOpt[k][j]);
                        for (plint k = 0; k < M; ++k) for (size_t j = 0; j < ubRow [k].size(); ++j) args.push_back(ubRow[k][j]);

                        grate.clear(); packed.clear(); status.clear();
                        const int e1 = optimize_cobrapy(pyFileName, args, here, grate, packed, status);
                        if (!anyCall) { erck = e1; anyCall = true; }
                        if (e1 != 0) { erck = e1; continue; }

                        for (plint k = 0; k < M; ++k) {
                            const double w = wRow[(size_t) k];
                            if (w <= 0.0) continue;
                            if (k >= (plint) status.size() || status[k] == 0) continue;   /* this
                                                                * organism could not live on this
                                                                * source here -- a physical
                                                                * statement, not a solver failure */
                            sAcc[(size_t) k] = status[k];
                            if (k < (plint) grate.size()) gAcc[(size_t) k] += (T) (w * grate[k]);
                            if (k < (plint) packed.size()) {
                                if (pAcc[(size_t) k].size() < packed[k].size())
                                    pAcc[(size_t) k].resize(packed[k].size(), T());
                                for (size_t j = 0; j < packed[k].size(); ++j)
                                    pAcc[(size_t) k][j] += (T) (w * packed[k][j]);
                            }
                        }
                    }

                    /* Hand the blended answer back under the names the rest of this function
                     * already uses, so the budget check, the repair pass and the unpacking below
                     * are untouched by any of the above. */
                    grate.assign(gAcc.begin(), gAcc.end());
                    packed = pAcc;
                    status = sAcc;

                    /* ---- budget check, and mark what to repair ----------------------------
                     * Identical in rule to the GLPK path: work out what every organism together
                     * would draw, and where that exceeds the supply, mark the substrate for every
                     * organism that asked for it and solve again against the pooled biomass.
                     *
                     * The fluxes are unpacked twice as a result -- once here to test the budget,
                     * once below to apply it -- which is a few multiplications against an 800 us
                     * Python round trip, and keeps the repair decision beside the loop it
                     * controls rather than threaded through the code that follows. */
                    if (M > 1) {
                        for (plint iS = 0; iS < subsNum; ++iS) {
                            if (cfg->fix_concentration[iS] || cfg->fix_lower_bounds[iS]) continue;
                            T draw = T();
                            for (plint k = 0; k < M; ++k) {
                                const plint iB = bLoc[k];
                                const plint gM = globalId[iB];
                                if (!(erck == 0 && k < (plint) status.size() && status[k] != 0)) continue;
                                std::vector<int> row(subsNum);
                                for (plint j = 0; j < subsNum; ++j) row[j] = (int) cfg->subsLoc[gM][j];
                                std::vector<T> full(subsNum, T());
                                expandFluxToSubstrates(packed[k], row, subsNum, full);
                                draw += cfg->fluxToDeltaC(full[iS], bmass[iB], gM, dt);
                            }
                            if (draw + avail[iS] < -thrd) {
                                bool any = false;
                                for (plint k = 0; k < M; ++k) {
                                    const plint iB = bLoc[k];
                                    const plint gM = globalId[iB];
                                    if (!(erck == 0 && k < (plint) status.size() && status[k] != 0)) continue;
                                    std::vector<int> row(subsNum);
                                    for (plint j = 0; j < subsNum; ++j) row[j] = (int) cfg->subsLoc[gM][j];
                                    std::vector<T> full(subsNum, T());
                                    expandFluxToSubstrates(packed[k], row, subsNum, full);
                                    if (full[iS] < -thrd) { shared[iB][iS] = true; any = true; }
                                }
                                if (any) needRepair = true;
                            }
                        }
                    }
                    ++pass;
                    }   /* while (needRepair && pass < 3) */

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

                        /* ---- the thermodynamic gate -------------------------------------
                         * Applied to the SOLVED fluxes, not to the uptake bounds.  A linear
                         * program asked for a smaller bound would redistribute its whole flux
                         * distribution and could hand back a different byproduct pattern, which
                         * is a modelling change; scaling the solution keeps the metabolic answer
                         * and asks only what fraction of it the local energy balance permits.
                         * The joint budget below therefore sees the gated fluxes.  A no-op when
                         * <thermodynamics> is off. */
                        if (gated) {
                            const double ft = complab_thermo::gateFor((int) gM, tconc);
                            if (ft < 1.0) {
                                mu *= (T) ft;
                                for (plint iS = 0; iS < subsNum; ++iS) flux[iB][iS] *= (T) ft;
                            }
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
                        /* [v1.3] The ternary was inverted: it added the running dC in the `free`
                         * branch, where dC[iS] is provably still zero (dC is cleared per voxel and
                         * each iS is visited once, so this loop is its only writer), and added
                         * nothing in the `total` branch, which is the one the comment above is
                         * about -- drawDown writes dC at CARRIER indices of an earlier component,
                         * which is exactly the second debit described. As written the guard
                         * reduced to the pristine avail[iS] it says it must not use. */
                        const T headroom = avail[iS] + (cfg->useTotals ? dC[iS] : T());
                        if (draw + headroom < T()) draw = -headroom;
                        /* [FIX] `if (draw > T()) draw = T();` used to stand here, and it was the
                         * one substantive difference between this path and the GLPK one.
                         *
                         * It discards every net RELEASE, so a COBRApy organism could consume but
                         * never excrete: an acetate-secreting E. coli produced no acetate, an AOM
                         * organism produced no sulfide, and the two back ends disagreed on any
                         * model that excretes -- while agreeing on uptake and growth, which is
                         * exactly the pattern that makes a difference look like a rounding
                         * question rather than a missing product.
                         *
                         * The GLPK path has no such line and never did. The line was a defence
                         * against a spurious source that the positivity clamp above already
                         * covers: a draw is limited by what the voxel holds, and a release is a
                         * flux the linear program returned, from mass the model balanced. */
                        if (cfg->useTotals && draw < T()) cfg->totals.drawDown(conc, iS, -draw, dC);
                        else                             dC[iS] += draw;
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
