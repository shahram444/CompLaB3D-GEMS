/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* =============================================================================
 * defineAbioticKinetics.hh  --  A + B -> C,  SECOND ORDER,  VOLUMETRIC
 *
 * The whole chemistry of the ABC sweep. One irreversible bimolecular reaction in
 * every pore voxel, no microbes, no minerals, no equilibrium:
 *
 *     A + B -> C          R = k [A][B]        mol/L/s
 *
 *     dA/dt = -R          dB/dt = -R          dC/dt = +R
 *
 *     C[0] = A   C[1] = B   C[2] = C
 *
 * THE GEOMETRY NEVER CHANGES. Nothing here precipitates or dissolves, and the
 * case does not enable <precipitation> or <dissolution>, so the pore space at
 * the last time step is the pore space that was read from geometry.dat. That is
 * deliberate: PRT-DeepONet-3D learns an operator from a FIXED pore space to a
 * concentration field, so the geometry has to be an input, not a variable.
 *
 * -----------------------------------------------------------------------------
 * WHY SECOND ORDER AND NOT FIRST
 *
 * A first-order sink R = -k[A] is separable: the answer is an exponential along
 * a streamline, the second species is decoration, and a network can fit it from
 * the Peclet number alone without ever looking at the pore space. There is no
 * transfer to learn.
 *
 * R = k[A][B] cannot be solved species by species. Where the two plumes have not
 * met, R is zero however much of either is present. The reaction lives on the
 * surface where they overlap, and WHERE THAT SURFACE SITS is set by the pore
 * structure and the flow together. That is the thing worth learning.
 *
 * -----------------------------------------------------------------------------
 * NOTHING STARTS IN THE DOMAIN
 *
 * Every <initial_concentration> in the case XML is 0.0. A enters at x = 0 and B
 * enters at x = nx-1, both held by Dirichlet planes, and they travel toward each
 * other. A is carried by the flow; B has to work upstream against it, so its
 * reach from the outlet is roughly L/Pe.
 *
 * The consequence for this file is that the rate is exactly zero everywhere at
 * t = 0, and stays zero in most of the domain for most of the run. Guard for
 * that rather than paying for it: the two early returns below skip the multiply
 * in every voxel a plume has not reached yet, which at the top of the Peclet
 * ladder is most of
 * them.
 *
 * -----------------------------------------------------------------------------
 * PARAMETERS FROM THE ENVIRONMENT
 *
 * One binary serves the whole sweep. The rate constant is read once from the
 * environment and cached, so a Damkohler series needs no rebuild:
 *
 *     PRT_KABIO     second-order rate constant k          L/(mol s)
 *     PRT_DT        must equal the printed [ADE] dt       s
 *     PRT_MAXFRAC   per-step reactant cap fraction        dimensionless
 *
 * tools/make_campaign.py writes all three into each case's env.sh, and
 * tools/check_campaign.py proves the binary actually reads them by running it
 * once under a sentinel value and looking for that value in its own log. A build
 * made against config/kinetics/defineAbioticKinetics.default.hh ignores the
 * environment completely: every case then runs identical chemistry while its
 * params.json records a different Damkohler, every run succeeds, the dataset
 * builds, and it teaches the network that Damkohler does nothing. That failure
 * is silent in every other way, which is why there is a check for it.
 *
 * The [KIN] line this file prints on its first call is the other half of that
 * check. Read it on the first run of any campaign.
 *
 * -----------------------------------------------------------------------------
 * ONE SCALING PATH, NOT TWO  [v1.3.2]
 *
 * v1.3.2 added <abiotic_rate_scale>, an XML multiplier applied to whatever this
 * function returns. It exists for cases that are NOT driven by a campaign, where
 * editing a header and rebuilding for each rate is the only alternative.
 *
 * It and PRT_KABIO are two ways to do one thing, and setting both multiplies the
 * rate twice. The shipped XML pins <abiotic_rate_scale> to 1.0 and the campaign
 * varies PRT_KABIO only. check_campaign.py refuses a campaign where both move.
 *
 * -----------------------------------------------------------------------------
 * WHAT DAMKOHLER MEANS FOR THIS LAW
 *
 * The startup banner reports Da_r = k_r L^2 / D with k_r = r(Cref)/Cref, the
 * pseudo-first-order constant. k L^2 / D is dimensionless only for a first-order
 * constant; a second-order k has units that do not cancel, so it has to be
 * collapsed against a reference composition first. Cref here is the FEED, A0,
 * because nothing starts in the domain:
 *
 *     r(A0, A0) = k A0^2        k_r = r/A0 = k A0        Da = k A0 L^2 / D
 *
 * so a campaign that wants a given Da sets
 *
 *     k = Da * D / (A0 * L^2)
 *
 * which is exactly what make_campaign.py computes. Stock v1.3.2 builds Cref from
 * <initial_concentration>, which is zero in every case here, so it would report
 * Da_r = 0 for the whole sweep. patch/apply_damkohler_feed_reference.py fixes that
 * to fall back to the Dirichlet feed value. Apply it, or read the banner knowing
 * it is reporting the wrong reference.
 * ============================================================================= */
#ifndef DEFINE_ABIOTIC_KINETICS_HH
#define DEFINE_ABIOTIC_KINETICS_HH

#include <vector>
#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <iostream>

namespace AbioticParams {

    inline double envd(const char* key, double dflt) {
        const char* v = std::getenv(key);
        if (v == 0 || *v == '\0') return dflt;
        try { return std::stod(std::string(v)); } catch (...) { return dflt; }
    }

    struct Params { double k_abio, dt_kinetics, MAX_RATE_FRACTION, MIN_CONC; };

    /* Read once, on the first voxel of the first step, and kept. getenv() per
     * voxel per step would be a measurable fraction of the run. */
    inline const Params& get() {
        static const Params p = [](){
            Params q;
            q.k_abio            = envd("PRT_KABIO",   1.0);
            q.dt_kinetics       = envd("PRT_DT",      1.0e-2);
            q.MAX_RATE_FRACTION = envd("PRT_MAXFRAC", 0.25);
            q.MIN_CONC          = 1.0e-20;
            std::cout << "[KIN] abiotic A+B->C  k_abio=" << q.k_abio
                      << " L/(mol s)  dt_kinetics=" << q.dt_kinetics
                      << " s  maxfrac=" << q.MAX_RATE_FRACTION << std::endl;
            return q;
        }();
        return p;
    }
}

namespace AbioticKineticsStats {
    static double iter_total_reaction = 0.0;
    static long   iter_cells_reacting = 0;
    inline void resetIteration() { iter_total_reaction = 0.0; iter_cells_reacting = 0; }
    inline void accumulate(double rate) {
        if (std::fabs(rate) > 1e-20) { iter_cells_reacting++; iter_total_reaction += rate; }
    }
}

/* -----------------------------------------------------------------------------
 * THE RATE LAW
 *
 * Returns rates per second; run_abiotic_kinetics multiplies by dt and hands the
 * increment to update_abiotic_rxnLattices, which clamps it against what the
 * voxel actually holds. The cap below is not that clamp and does not replace it:
 * the cap keeps the step from overshooting so far that the clamp has to fire,
 * because a clamped step is no longer the reaction the case was configured for.
 * The [CLAMP] line at the end of the run says how often that happened.
 * ----------------------------------------------------------------------------- */
void defineAbioticRxnKinetics(
    std::vector<double> C,          // C[0]=A  C[1]=B  C[2]=C
    std::vector<double>& subsR,     // OUTPUT substrate rates, mol/L/s
    plb::plint mask
) {
    (void) mask;                    // volumetric: every pore voxel, no surface gate
    const AbioticParams::Params& K = AbioticParams::get();

    for (size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;

    /* SIZE GUARD, on the output as well as the input. A guard left over from a
     * longer species list disables this file completely and silently. */
    if (C.size() < 3 || subsR.size() < 3) return;

    const double A = std::max(C[0], 0.0);
    const double B = std::max(C[1], 0.0);

    /* Where the two plumes have not met yet, there is nothing to do. This is the
     * common case for most of the domain for most of the run. */
    if (A <= K.MIN_CONC || B <= K.MIN_CONC) return;

    double R = K.k_abio * A * B;

    /* THE PER-STEP CAP.
     *
     * A and B are consumed one for one, so the binding reactant is whichever is
     * scarcer and ONE cap covers both. Capping them separately would be the same
     * number here, but it would stop being the same number the moment the
     * stoichiometry changes, and a cap that silently breaks stoichiometry is
     * worse than no cap. */
    if (K.dt_kinetics > 0.0) {
        const double scarcest = (A < B) ? A : B;
        const double capR = scarcest * K.MAX_RATE_FRACTION / K.dt_kinetics;
        if (R > capR) R = capR;
    }

    subsR[0] = -R;      // A consumed
    subsR[1] = -R;      // B consumed
    subsR[2] = +R;      // C produced

    AbioticKineticsStats::accumulate(R);

    for (size_t i = 0; i < subsR.size(); ++i)
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
}

/* -----------------------------------------------------------------------------
 * THE DISSOLUTION HOOK -- REQUIRED, AND INERT HERE
 *
 * v1.3 calls this from dissolutionVOP.hh for every voxel carrying a declared
 * solid phase, and the call is compiled unconditionally, so the symbol must
 * exist in any defineAbioticKinetics.hh even for a case with no mineral. This
 * sweep has none: the geometry is frozen on purpose, so nothing dissolves.
 *
 * It sets mineralR and touches nothing else. It deliberately does NOT zero
 * subsR, although the shipped default does: the caller already cleared it
 * (dissolutionVOP.hh:471), so clearing it again is redundant, and a hook that
 * wipes an accumulator it did not fill is the wrong shape for one that might
 * later be called after another contributor.
 * ----------------------------------------------------------------------------- */
inline void defineDissolutionRate(plb::plint phaseId, const std::vector<double>& C,
                                  double mineral, std::vector<double>& subsR,
                                  double& mineralR, plb::plint mask)
{
    (void) phaseId; (void) C; (void) mineral; (void) subsR; (void) mask;
    mineralR = 0.0;
}

#endif // DEFINE_ABIOTIC_KINETICS_HH
