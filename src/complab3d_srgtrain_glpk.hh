/* This file is a part of the CompLaB program.  AGPL-3.0-or-later.
 * Meile Lab, University of Georgia.  shahram.asgari@uga.edu
*/

/* ================================================================================================
 * complab3d_srgtrain_glpk.hh  --  THE CALLBACK THAT LETS SURROGATE TRAINING REACH THE LP SOLVER
 * ================================================================================================
 *
 *  complab3d_surrogate.hh knows how to sweep a function and fit a network to it, but deliberately
 *  knows nothing about GLPK: it takes the solver as a callback, which is what lets it be tested
 *  against an analytic function with no LP solver present.
 *
 *  This is that callback, for the GLPK path. Twenty lines of substance: set the uptake bounds on
 *  the microbe's already-built glp_prob, solve, return the objective.
 *
 *  ------------------------------------------------------------------------------------------------
 *  IT USES THE SAME PROBLEM OBJECT THE SIMULATION WILL USE
 *
 *  Not a copy, and not a freshly built one. cfg.vec_lp[microbe] is the persistent glp_prob that
 *  setup_metabolic_solvers() created, with the model's own bounds already applied. Training through
 *  it means the network is fitted to exactly the linear program the run will solve -- including any
 *  <constraint_indices> overrides, which a separately built problem would miss.
 *
 *  The bounds are restored afterwards, so training leaves the problem exactly as it found it.
 *
 *  ------------------------------------------------------------------------------------------------
 *  WHAT "UPTAKE" MEANS HERE
 *
 *  The trainer passes POSITIVE uptake rates, one per training input, in mmol/gDW/h -- the same
 *  convention as Fin everywhere else in CompLaB. Exchange reactions are written "metabolite ->",
 *  so consuming is negative flux, and the lower bound is set to -uptake. Getting that sign wrong
 *  would train the network on an organism that excretes everything, so it is done in exactly one
 *  place: here.
 *
 *  Substrates that are not training inputs keep whatever bound the model gave them. That is
 *  deliberate: a genome-scale model's medium is part of the model, and silently closing it would
 *  make every solve infeasible.
 * ================================================================================================
 */
#ifndef COMPLAB3D_SRGTRAIN_GLPK_HH
#define COMPLAB3D_SRGTRAIN_GLPK_HH

#ifdef COMPLAB_ENABLE_GLPK

#include <vector>

#include "complab3d_metabolic.hh"
#include "complab3d_glpkcpp.hh"
#include "complab3d_surrogate.hh"

namespace complab_srgtrain {

/* Everything the callback needs, bundled so it can travel through a void*. */
struct GlpkContext {
    MetabolicConfig *cfg;
    plint microbe;                    // global microbe index
    std::vector<plint> inputCol;      // LP column for each training input, -1 if unmapped
    long calls;
    long infeasible;

    GlpkContext() : cfg(0), microbe(0), calls(0), infeasible(0) {}
};

/* Map each training input, named by substrate, onto that microbe's exchange column.
 *
 * Returns "" on success. A name that does not match any substrate is fatal rather than ignored:
 * training on a bound that is never applied would fit a network to a constant, which looks like a
 * beautifully converged result. */
inline std::string bindInputs(GlpkContext &ctx,
                              const std::vector<complab_srg::InputSpec> &inputs,
                              const std::vector<std::string> &substrateNames)
{
    ctx.inputCol.assign(inputs.size(), -1);
    for (size_t k = 0; k < inputs.size(); ++k) {
        plint sIdx = -1;
        for (size_t s = 0; s < substrateNames.size(); ++s)
            if (substrateNames[s] == inputs[k].name) { sIdx = (plint) s; break; }
        if (sIdx < 0)
            return "training input '" + inputs[k].name + "' is not one of the substrates in "
                   "CompLaB.xml. It must name a <name_of_substrates>.";

        const plint col = ctx.cfg->subsLoc[(size_t) ctx.microbe][(size_t) sIdx];
        if (col < 0)
            return "training input '" + inputs[k].name + "' has no exchange reaction for this "
                   "microbe (its <exchange_reaction_indices> entry is -1), so its uptake bound "
                   "would never be applied.";
        ctx.inputCol[k] = col;
    }
    return "";
}

/* The callback itself. Signature fixed by complab_srg::FbaFn. */
inline double solveWithUptake(const std::vector<double> &uptake, void *vctx)
{
    GlpkContext *ctx = (GlpkContext *) vctx;
    if (!ctx || !ctx->cfg) return 0.0;
    MetabolicConfig &cfg = *ctx->cfg;
    const size_t gM = (size_t) ctx->microbe;
    glp_prob *lp = cfg.vec_lp[gM];
    if (!lp) return 0.0;

    ++ctx->calls;

    const int ncol = cfg.vec_nrxns[gM];
    const int nrow = cfg.vec_nmets[gM];

    /* Remember what we are about to change, so the problem is left as we found it. */
    std::vector<int> touched;
    std::vector<double> savedLo, savedHi;
    touched.reserve(ctx->inputCol.size());

    for (size_t k = 0; k < ctx->inputCol.size() && k < uptake.size(); ++k) {
        const int col = (int) ctx->inputCol[k];
        if (col < 0 || col >= ncol) continue;
        touched.push_back(col);
        savedLo.push_back(glp_get_col_lb(lp, col + 1));      // GLPK columns are 1-based
        savedHi.push_back(glp_get_col_ub(lp, col + 1));

        /* POSITIVE uptake in, NEGATIVE lower bound out. The upper bound is left alone so the
         * solver may still excrete if the network wants to. */
        const double lo = -uptake[k];
        const double hi = savedHi.back();
        glp_set_col_bnds(lp, col + 1, (lo <= hi) ? GLP_DB : GLP_FX, lo, (lo <= hi) ? hi : lo);
    }

    std::vector<double> xmin((size_t) ncol, 0.0), lambda((size_t) nrow, 0.0),
                        redcosts((size_t) ncol, 0.0);
    std::vector<int> location;
    std::vector<double> emptyLb, emptyUb;
    double fmin = 0.0, tsec = 0.0, tmem = 0.0;
    int status = 0;

    const int err = run_glpk(lp, cfg.method[gM], cfg.isMIP[gM], nrow, ncol, cfg.lpsolver,
                             cfg.sParam[gM], cfg.iParam[gM],
                             xmin, fmin, status, location, emptyLb, emptyUb,
                             lambda, redcosts, tsec, tmem);

    for (size_t k = 0; k < touched.size(); ++k) {
        const double lo = savedLo[k], hi = savedHi[k];
        glp_set_col_bnds(lp, touched[k] + 1, (lo < hi) ? GLP_DB : GLP_FX, lo, hi);
    }

    if (err != 0 || status != GLP_OPT) {
        /* Infeasible means the organism cannot balance at that supply. Zero growth is the
         * physically right label, and keeping the sample teaches the network where the feasible
         * region ends instead of leaving a hole in its domain. */
        ++ctx->infeasible;
        return 0.0;
    }
    return (fmin > 0.0) ? fmin : 0.0;
}

/* ================================================================================================
 *  A LINEAR PROGRAM BUILT JUST FOR TRAINING
 *
 *  The common case is a run with ONE organism whose <reaction_type> is `surrogate`. There is then
 *  no GLPK microbe anywhere in the run, so there is no glp_prob to sweep -- and requiring the user
 *  to add a throwaway GLPK microbe just to be allowed to train would be an absurd thing to ask.
 *
 *  So build one here, for training only, from the surrogate microbe's own <model_filename> (or the
 *  file <model_source> resolved) and its own <exchange_reaction_names>. It is released as soon as
 *  training finishes; the simulation never sees it and never pays for it.
 *
 *  This deliberately reuses load_metabolic_models3D() and setup_metabolic_solvers() rather than
 *  building the problem by hand: the network must be fitted to the SAME linear program the GLPK
 *  path would have solved, down to the bound overrides, or the surrogate is a fit to something
 *  that was never going to be run.
 * ================================================================================================ */
struct StandaloneLp {
    MetabolicConfig cfg;
    bool built;

    StandaloneLp() : built(false) {}
    ~StandaloneLp() { if (built) release_metabolic_solvers(cfg); }

private:
    StandaloneLp(const StandaloneLp &);              // MetabolicConfig owns raw handles
    StandaloneLp &operator=(const StandaloneLp &);
};

/* Copy one microbe's metabolic settings out of the run's config, load its model, and build the LP.
 * Returns "" on success. `src` is the run's own MetabolicConfig, `iM` the microbe to copy. */
inline std::string buildTrainingLp(StandaloneLp &S, const MetabolicConfig &src, plint iM,
                                   const std::string &inputDir, plint numSubstrates)
{
    MetabolicConfig &c = S.cfg;

    c.enable_fba_glpk = true;
    c.glpk_count = 1;
    c.mm_count   = 1;
    c.lpsolver   = src.lpsolver;
    c.save_pb    = 0;

    const size_t m = (size_t) iM;
    c.model_filename.assign(1, src.model_filename[m]);
    c.subsLoc.assign(1, src.subsLoc[m]);
    c.exchange_names.assign(1, src.exchange_names[m]);
    c.constraint_loc.assign(1, src.constraint_loc[m]);
    c.constraint_lb .assign(1, src.constraint_lb[m]);
    c.constraint_ub .assign(1, src.constraint_ub[m]);
    c.equate_bounds .assign(1, src.equate_bounds[m]);
    c.objDir        .assign(1, src.objDir[m]);
    c.vmax          .assign(1, src.vmax[m]);
    c.maxUptake     .assign(1, src.maxUptake[m]);
    c.maxRelease    .assign(1, src.maxRelease[m]);
    c.biomass_molar_mass.assign(1, src.biomass_molar_mass[m]);
    c.fix_concentration.assign((size_t) numSubstrates, false);
    c.fix_lower_bounds .assign((size_t) numSubstrates, false);

    c.S3.assign(1, std::vector< std::vector<T> >());
    c.S1.assign(1, std::vector<T>());
    c.vec_b .assign(1, std::vector<T>());
    c.vec_c .assign(1, std::vector<T>());
    c.vec_lb.assign(1, std::vector<T>());
    c.vec_ub.assign(1, std::vector<T>());
    c.vec_objLoc.assign(1, 0);
    c.vec_nmets .assign(1, 0);
    c.vec_nrxns .assign(1, 0);

    if (c.model_filename[0].empty())
        return "microbe" + std::to_string((long long) iM) + " has no <model_filename>, and "
               "<model_source> supplied none either, so there is no metabolic model to sweep.";

    std::vector<plint> rt(1, (plint) rxntype::GLPK);
    if (load_metabolic_models3D(c, inputDir, rt, 1) != 0)
        return "the metabolic model could not be read (see above).";

    char pyNone[] = "none";
    char srcNone[] = ".";
    if (setup_metabolic_solvers(c, rt, 1, pyNone, srcNone) != 0)
        return "the linear program could not be built (see above).";

    S.built = true;
    return "";
}

}  // namespace complab_srgtrain

#endif  // COMPLAB_ENABLE_GLPK
#endif  // COMPLAB3D_SRGTRAIN_GLPK_HH
