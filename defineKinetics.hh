/* =============================================================================
 * defineKinetics.hh   --   BIOTIC RATE LAWS.  This file is yours to edit.
 *
 * CompLB3D compiles this in. It is a TEMPLATE: every capability is shown, live,
 * with generic placeholder values. Nothing here is commented out and nothing is
 * switched off, so you can see the whole shape of what a rate law may do and
 * then replace the numbers with yours.
 *
 * The substrate and microbe ORDER is the order they are declared in
 * CompLaB.xml. C[0] is <substrate0>, B[0] is <microbe0>, and so on. Getting
 * that order wrong is the one mistake this file cannot catch for you.
 *
 * UNITS
 *     C[i]      mol/L          substrate concentration
 *     B[i]      mol/L          biomass
 *     subsR[i]  mol/L/s        rate of change of substrate i
 *     bioR[i]   1/s            SPECIFIC growth rate of microbe i (per unit biomass)
 *     rxnR[k]   mol/L/s        optional: one entry per reaction, for <save_reaction_rates>
 *
 * WHEN THIS FILE IS CALLED
 *     Only for a microbe whose <reaction_type> is kinetics, or one of the
 *     _and_kinetics combinations. A microbe on glpk, cobrapy, surrogate,
 *     symbolic or graphnet does not come through here.
 * ============================================================================= */
#ifndef DEFINE_KINETICS_HH
#define DEFINE_KINETICS_HH

#include <vector>
#include <cmath>
#include <algorithm>

/* -----------------------------------------------------------------------------
 * REACTION RATE OUTPUT  -  the opt-in contract from src/complab3d_rates.hh
 *
 * With <save_reaction_rates>true</> in CompLaB.xml, the solver writes an extra
 * .vti per snapshot holding how fast each reaction ran in every voxel.  Per
 * chemical and per microbe channels appear on their own, because those are
 * sized from CompLaB.xml.  PER REACTION channels appear only if this header
 * says how many reactions it has and what they are called -- which is what the
 * three things below do.  Name them in the order you fill rxnR.
 * --------------------------------------------------------------------------- */
#define COMPLAB_HAS_RXN_RATES 1
#include <string>
namespace RxnRates {
    inline const std::vector<std::string>& names() {
        static const std::vector<std::string> n = {
            "growth"          /* rxnR[0] : the single reaction this template defines */
        };
        return n;                      /* as many as you have, in your own order */
    }
}

namespace KineticParams {
    /* -------- generic placeholders: replace with your own ------------------ */
    const double mu_max   = 2.78e-5;   /* 1/s      max specific growth rate (= 0.1 /h)      */
    const double K_donor  = 1.0e-4;    /* mol/L    half saturation, electron donor          */
    const double K_accept = 5.0e-5;    /* mol/L    half saturation, electron acceptor       */
    const double Y_donor  = 0.42;      /* mol biomass per mol donor    (yield)              */
    const double s_accept = 4.8;       /* mol acceptor per mol donor   (stoichiometry)      */
    const double s_prod   = 1.845;     /* mol product  per mol donor   (stoichiometry)      */
    const double k_decay  = 6.94e-8;   /* 1/s      first order decay   (= 0.006 /day)       */
    const double MAX_FRAC = 0.5;       /* never consume more than half a pool in one step   */
    const double dt_kin   = 1.0e-2;    /* s        the reaction step; match <ade_dt>        */
    const double eps      = 1.0e-30;   /* guards a division, nothing more                   */
}

/* Index names, so the body below reads as chemistry rather than as arithmetic.
 * Rename these to your species and the rate law still says what it means. */
namespace KineticIdx {
    const int DONOR = 0;    /* C[0]  e.g. acetate  */
    const int ACCEP = 1;    /* C[1]  e.g. Fe(III)  */
    const int PROD  = 2;    /* C[2]  e.g. HCO3-    */
    const int MIC   = 0;    /* B[0]  the organism  */
}

void defineRxnKinetics(
    std::vector<double> B,          /* IN : biomass here, one entry per microbe            */
    std::vector<double> C,          /* IN : concentrations here, one entry per substrate   */
    std::vector<double>& subsR,     /* OUT: mol/L/s per substrate                          */
    std::vector<double>& bioR,      /* OUT: 1/s per microbe                                */
    plb::plint mask,                /* IN : voxel label; the processor decides where       */
    std::vector<double>* rxnR = 0   /* OUT: optional, one entry per reaction               */
) {
    using namespace KineticParams;
    using namespace KineticIdx;

    for (size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    for (size_t i = 0; i < bioR.size();  ++i) bioR[i]  = 0.0;
    if (rxnR) for (size_t i = 0; i < rxnR->size(); ++i) (*rxnR)[i] = 0.0;

    (void) mask;   /* available if you want a different law inside biofilm than in pore */

    if ((int) C.size() <= PROD || (int) B.size() <= MIC) return;

    const double donor = std::max(C[DONOR], 0.0);
    const double accep = std::max(C[ACCEP], 0.0);
    const double bio   = std::max(B[MIC],   0.0);

    /* ---- dual Monod growth -------------------------------------------------
     * Both substrates limit: run out of either and growth stops. */
    const double mu = mu_max
                    * donor / (K_donor  + donor + eps)
                    * accep / (K_accept + accep + eps);

    /* ---- the reaction rate, per volume ------------------------------------
     *     donor + s_accept * acceptor  ->  s_prod * product
     * Consumption follows growth through the yield, so the carbon balance is
     * a property of the numbers above rather than of three separate rate laws. */
    double R = (mu * bio) / (Y_donor + eps);        /* mol donor / L / s */

    /* ---- never take more than a pool holds in one step --------------------
     * Scale the whole reaction, so the stoichiometry survives the limiter. */
    const double capD = donor * MAX_FRAC / dt_kin;
    const double capA = accep * MAX_FRAC / (dt_kin * s_accept);
    R = std::min(R, std::min(capD, capA));
    if (R < 0.0) R = 0.0;

    subsR[DONOR] = -R;
    subsR[ACCEP] = -R * s_accept;
    subsR[PROD]  = +R * s_prod;

    /* ---- growth, and decay when there is nothing to eat --------------------
     * bioR is a SPECIFIC rate: the solver multiplies it by the local biomass. */
    bioR[MIC] = (R > 0.0) ? (R * Y_donor / (bio + eps)) : (-k_decay);

    /* Report the per-reaction rate, in the order RxnRates::names() declares.
     *
     * A call with rxnR set is a REPORTING call made by the sampler, not a step
     * of the simulation.  If this header ever keeps running totals, guard them
     * with  if (!rxnR)  so a reporting call cannot inflate them. */
    if (rxnR && rxnR->size() > 0) (*rxnR)[0] = R;
}

#endif
