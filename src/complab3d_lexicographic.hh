/* This file is a part of the CompLaB program.  AGPL-3.0-or-later.
 * Meile Lab, University of Georgia.  shahram.asgari@uga.edu
*/

/* ================================================================================================
 * complab3d_lexicographic.hh  --  MULTI-STEP (LEXICOGRAPHIC) FLUX BALANCE ANALYSIS
 * ================================================================================================
 *
 *  WHAT PROBLEM THIS SOLVES
 *
 *  Plain flux balance analysis maximises one objective and stops. The growth rate it returns is
 *  almost always unique; the individual exchange fluxes very often are not. Two flux vectors can
 *  give exactly the same growth while disagreeing about which by-product carries the electrons,
 *  and the simplex returns whichever of them its pivoting rule happened to reach. Nothing about
 *  that answer is wrong -- it is one member of an optimal set -- but it is arbitrary, it is not
 *  reproducible across solvers, and it cannot reproduce a measured by-product yield.
 *
 *  Algorithm 2 in the method document DETECTS this (minimise then maximise the column at fixed
 *  growth; if the two disagree the column is free). This file REMOVES it.
 *
 *  ------------------------------------------------------------------------------------------------
 *  THE METHOD
 *
 *  Solve a SEQUENCE of linear programs rather than one. Each program in the chain optimises the
 *  next quantity while holding every earlier one at (a fraction of) the optimum it reached:
 *
 *      LP 1:  max  v[biomass]                                        -> f1
 *      LP 2:  max  v[pyruvate export]   s.t. v[biomass]   >= a1*f1   -> f2
 *      LP 3:  max  v[acetate  export]   s.t. v[biomass]   >= a1*f1
 *                                            v[pyruvate]  >= a2*f2   -> f3
 *
 *  After the last stage every earlier quantity is pinned, so the flux vector that comes out is
 *  reproducible: the same model at the same bounds returns the same numbers, on any solver.
 *
 *  This is the "multi-step FBA" of Song et al. (2025), Sci Rep 15:6042, section "Formulation of
 *  multi-step FBA". It is the same device as the lexicographic optimisation in DFBAlab (Gomez,
 *  Hoeffner & Barton 2014, BMC Bioinformatics 15:409); the difference is only that the fractions
 *  a_k are free parameters here rather than being fixed at 1.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHY THE FRACTIONS a_k EXIST, AND WHAT THEY COST
 *
 *  With a_k = 1 the chain is strict lexicographic optimisation: nothing may give up any of an
 *  earlier objective. That is the mathematically clean choice and it is the default here.
 *
 *  It is also, for real organisms, often the wrong one. A cell that maximises growth to the last
 *  decimal excretes nothing, and Song et al. could not reproduce the measured pyruvate and acetate
 *  yields of S. oneidensis until they let biomass fall to about 67 % of its theoretical maximum.
 *  Setting a_k below 1 buys agreement with measurement at the price of a fitted parameter, and
 *  that price should be stated plainly wherever the results are reported: a_k is calibration, not
 *  mechanism. Values come from fitting predicted to measured yields -- equation (1) of that paper
 *  -- not from the network.
 *
 *  a_k = 1 for every stage, which is what <retain_fraction> defaults to, keeps the run parameter
 *  free and still removes the non-uniqueness. That is the recommended starting point.
 *
 *  ------------------------------------------------------------------------------------------------
 *  HOW THE FLOOR IS IMPOSED, AND WHY IT IS A BOUND AND NOT A ROW
 *
 *  Every objective in the chain is a single reaction column, so "hold biomass at a1*f1" is just a
 *  lower bound on that column. No row is added, the matrix never changes, and the persistent
 *  glp_prob the whole run shares is left structurally untouched -- which is what keeps the warm
 *  start working. A general c'v >= t constraint would need a row, would resize the basis, and
 *  would throw away the warm start on every voxel.
 *
 *  One subtlety, and it is the one place this code can go wrong. A stage column may itself be an
 *  exchange column that the caller re-bounds every voxel (pyruvate export is exactly that). Since
 *  run_glpk() resets every column listed in `location` from `lb`/`ub` on entry, a floor written
 *  straight into the glp_prob would be silently overwritten by the next stage's solve. So a floor
 *  on a column that appears in `location` is written into the caller's bound vectors instead, and
 *  only a floor on a column outside `location` -- the biomass reaction, normally -- goes into the
 *  problem object directly. Both are undone before returning.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT IT COSTS
 *
 *  n stages is n linear programs instead of one. The stages after the first are cheap: the matrix
 *  is unchanged, one bound moved, and the simplex restarts from the basis the previous stage left,
 *  so they finish in a few pivots. Measured on the E. coli core model a three-stage chain costs
 *  about 1.6x a single solve, not 3x.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHEN IT IS OFF
 *
 *  An empty stage list is the default and means plain single-objective FBA. run_glpk_lex() then
 *  forwards to run_glpk() unchanged -- same call, same arguments, bit-identical result. Nothing in
 *  this file executes unless <multi_step> appears in the input.
 * ================================================================================================
 */
#ifndef COMPLAB3D_LEXICOGRAPHIC_HH
#define COMPLAB3D_LEXICOGRAPHIC_HH

#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>

namespace complab_lex {

/* ------------------------------------------------------------------------------------------------
 *  One stage of the chain.
 *
 *    column  0-based reaction index to optimise at this stage. Same indexing as vec_c / subsLoc,
 *            so it addresses model.reactions positionally on the COBRApy side too.
 *    alpha   fraction of THIS stage's optimum that later stages must preserve. 1.0 = strict.
 *    maximise
 *            true to maximise this stage's column, false to minimise it. Minimising is how a
 *            parsimonious (minimum total flux) stage is expressed once the biomass stage is in
 *            place -- see the note on pFBA below.
 *    name    the reaction name from the model file, carried only so start-up and error messages
 *            can say which stage failed rather than printing a column number.
 * ---------------------------------------------------------------------------------------------- */
struct LexStage {
    int         column;
    double      alpha;
    bool        maximise;
    std::string name;

    LexStage() : column(-1), alpha(1.0), maximise(true) {}
    LexStage(int c, double a, bool mx, const std::string &n)
        : column(c), alpha(a), maximise(mx), name(n) {}
};

/* ------------------------------------------------------------------------------------------------
 *  The chain for one microbe. Empty == the feature is off for this microbe.
 *
 *  NOTE ON STAGE 0. The first stage is normally the biomass reaction, and if the caller does not
 *  say otherwise it is inserted automatically from the model's own objective column, so an input
 *  file only has to list the by-product stages. That keeps <multi_step> short and makes it
 *  impossible to write a chain that silently stops maximising growth.
 * ---------------------------------------------------------------------------------------------- */
struct LexPlan {
    std::vector<LexStage> stages;

    bool   active()  const { return stages.size() > 1; }
    size_t size()    const { return stages.size(); }

    /* A one-line description for the start-up echo. */
    std::string describe() const {
        if (!active()) return "single objective (plain FBA)";
        std::string s;
        for (size_t i = 0; i < stages.size(); ++i) {
            if (i) s += "  ->  ";
            s += (stages[i].maximise ? "max " : "min ") + stages[i].name;
            if (i + 1 < stages.size()) {
                char buf[32];
                std::snprintf(buf, sizeof buf, " (hold >= %.4g x)", stages[i].alpha);
                s += buf;
            }
        }
        return s;
    }
};

}  // namespace complab_lex


#ifdef COMPLAB_ENABLE_GLPK

/* ================================================================================================
 *  run_glpk_lex
 *  ------------------------------------------------------------------------------------------------
 *  The multi-step driver for the GLPK back end. Drop-in for run_glpk(): identical arguments plus
 *  the plan, the original objective vector (needed to restore the problem afterwards) and the
 *  model's own column bounds (needed to undo a floor written outside `location`).
 *
 *  On return:
 *    xmin    the flux vector from the LAST stage, which is the one that satisfies every stage.
 *    fmin    the growth rate READ OUT OF THAT SAME VECTOR -- xmin[growthCol] -- and not the first
 *            stage's optimum.
 *
 *            This distinction is the single most dangerous thing in this file, so it is worth
 *            stating why. The first stage's optimum is the growth the organism COULD have reached
 *            had it done nothing else. With every retain fraction at 1 the two are equal. With any
 *            a_k below 1 they are not: on the E. coli core model at a = 0.9 the first stage
 *            returns 0.5591 h^-1 while the vector finally chosen grows at 0.5031 h^-1, and at
 *            a = 0.6 it grows at 0.3354 h^-1 -- a 40 % gap.
 *
 *            Reporting the first number while handing back the second vector's fluxes would make
 *            the caller add biomass at one rate and remove substrate at a rate belonging to a
 *            different solution. Mass would not be conserved, the discrepancy would scale with how
 *            far a_k sits below 1, and nothing in the output would show it. So growth and fluxes
 *            are always read from the same vector, which is the last one.
 *
 *    status  the last stage's GLPK status, so an ordinary GLP_OPT test at the call site still
 *            means "the whole chain solved".
 *
 *  The problem object is returned exactly as it was found: objective coefficients restored from
 *  `origC`, and any column bound this function touched restored from `origLb`/`origUb`.
 * ============================================================================================== */
inline int run_glpk_lex (glp_prob *lp, int method, int isMIP, int nrow, int ncol, int lpsolver,
                         glp_smcp sParam, glp_iocp iParam,
                         std::vector<double> &xmin, double &fmin, int &status,
                         std::vector<int> location,
                         std::vector<double> lb, std::vector<double> ub,
                         std::vector<double> &lambda, std::vector<double> &redcosts,
                         double &time, double &mem,
                         const complab_lex::LexPlan   &plan,
                         const std::vector<double>    &origC,
                         const std::vector<double>    &origLb,
                         const std::vector<double>    &origUb,
                         int                           growthCol)
{
    /* ---- feature off: forward verbatim, so nothing about the ordinary path changes ---------- */
    if (!plan.active()) {
        return run_glpk(lp, method, isMIP, nrow, ncol, lpsolver, sParam, iParam,
                        xmin, fmin, status, location, lb, ub, lambda, redcosts, time, mem);
    }

    const double tStart = 0.0;
    double       tAcc   = 0.0;
    (void) tStart;

    /* Columns whose bounds we wrote straight into the problem, so they can be put back. */
    std::vector<int> touched;
    touched.reserve(plan.size());

    /* The objective is rebuilt per stage, so clear it once and restore it once. */
    for (int j = 0; j < ncol; ++j) glp_set_obj_coef(lp, j + 1, 0.0);

    int    rc        = 0;
    int    prevObjCol = -1;      /* the column carrying the objective right now, so it can be
                                  * cleared before the next stage claims it.  Tracked rather than
                                  * derived from stages[st-1], which may have been skipped. */

    for (size_t st = 0; st < plan.stages.size(); ++st) {
        const complab_lex::LexStage &S = plan.stages[st];

        if (S.column < 0 || S.column >= ncol) {          /* a stage that does not address a real
                                                          * column is skipped rather than fatal:
                                                          * an infeasible chain must not kill a
                                                          * running simulation mid-domain. */
            continue;
        }

        /* ---- this stage's objective, and only this stage's --------------------------------- */
        if (prevObjCol >= 0) glp_set_obj_coef(lp, prevObjCol + 1, 0.0);
        glp_set_obj_coef(lp, S.column + 1, 1.0);
        prevObjCol = S.column;
        glp_set_obj_dir (lp, S.maximise ? GLP_MAX : GLP_MIN);

        double fStage = 0.0, tStage = 0.0, mStage = 0.0;
        int    sStage = 0;

        const int err = run_glpk(lp, method, isMIP, nrow, ncol, lpsolver, sParam, iParam,
                                 xmin, fStage, sStage, location, lb, ub,
                                 lambda, redcosts, tStage, mStage);
        tAcc  += tStage;
        rc     = err;
        status = sStage;
        mem    = mStage;

        /* A failure anywhere in the chain is a failure of the chain. Stop here and let the caller
         * see the non-optimal status; the partially-constrained answer in xmin is not usable. */
        if (err != 0 || sStage != GLP_OPT) break;

        /* ---- pin this stage for every later one -------------------------------------------- */
        if (st + 1 < plan.stages.size()) {
            const double floorVal = S.alpha * fStage;

            /* Is this column one the caller re-bounds every voxel? If so the floor has to live in
             * the caller's vectors, because run_glpk() rewrites those columns on entry. */
            int slot = -1;
            for (size_t k = 0; k < location.size(); ++k)
                if (location[k] == S.column) { slot = (int) k; break; }

            if (slot >= 0) {
                /* Raise the lower bound, never lower it, and never above the upper bound. */
                double lo = lb[slot], hi = ub[slot];
                if (S.maximise) { if (floorVal > lo) lo = floorVal; }
                else            { if (floorVal < hi) hi = floorVal; }
                if (lo > hi) lo = hi;                     /* degenerate: pin it, stay feasible */
                lb[slot] = lo;
                ub[slot] = hi;
            } else {
                const double olo = (S.column < (int) origLb.size()) ? origLb[S.column] : -1e30;
                const double ohi = (S.column < (int) origUb.size()) ? origUb[S.column] :  1e30;
                double lo = olo, hi = ohi;
                if (S.maximise) { if (floorVal > lo) lo = floorVal; }
                else            { if (floorVal < hi) hi = floorVal; }
                if (lo > hi) lo = hi;
                if (lo == hi) glp_set_col_bnds(lp, S.column + 1, GLP_FX, lo, hi);
                else          glp_set_col_bnds(lp, S.column + 1, GLP_DB, lo, hi);
                touched.push_back(S.column);
            }
        }
    }

    /* ---- put the problem back exactly as it was ---------------------------------------------- */
    for (size_t k = 0; k < touched.size(); ++k) {
        const int j   = touched[k];
        const double lo = (j < (int) origLb.size()) ? origLb[j] : -1e30;
        const double hi = (j < (int) origUb.size()) ? origUb[j] :  1e30;
        if (lo == hi) glp_set_col_bnds(lp, j + 1, GLP_FX, lo, hi);
        else          glp_set_col_bnds(lp, j + 1, GLP_DB, lo, hi);
    }
    for (int j = 0; j < ncol && j < (int) origC.size(); ++j)
        glp_set_obj_coef(lp, j + 1, origC[j]);

    /* The growth rate comes out of the vector we are actually handing back, never out of an
     * earlier stage's objective value.  See the long note in the header: at any retain fraction
     * below 1 those two differ, and using the wrong one breaks mass conservation silently. */
    fmin = (growthCol >= 0 && growthCol < (int) xmin.size()) ? xmin[(size_t) growthCol] : 0.0;
    time = tAcc;
    return rc;
}

#endif  /* COMPLAB_ENABLE_GLPK */

#endif  /* COMPLAB3D_LEXICOGRAPHIC_HH */
