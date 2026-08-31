/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* =============================================================================
 * defineAbioticKinetics.hh  --  EXAMPLE 15: precipitation AND dissolution
 *
 * Both hooks are implemented, on the same mineral, so a voxel can seal and
 * later reopen:
 *
 *     precipitation   Fe2+ + HS-  ->  FeS(s)          rate = kp * [Fe2+] * [HS-]
 *     dissolution     FeS(s) + H+ ->  Fe2+ + HS-      rate = kd * [H+] * mineral
 *
 * Substrate order, matching CompLaB.xml:
 *     C[0] = Fe2   C[1] = HS   C[2] = H   C[3] = FeS  (immobile, the mineral)
 *
 * WHY <reopen_fraction> MATTERS HERE AND NOWHERE ELSE
 *   With both directions active a voxel can sit near the sealing threshold. A
 *   single threshold would make it flip open and shut every update interval,
 *   and every flip re-solves the flow. reopen_fraction 0.9 puts a gap between
 *   sealing at full_density and reopening below 0.9 x full_density, so it
 *   settles instead of chattering. Watch the porosity in the log: it should
 *   move in steps and then hold, not oscillate.
 *
 * THE TWO HOOKS ARE NOT INVERSES. Precipitation is second order in the
 * dissolved ions; dissolution is first order in acid and in remaining mineral.
 * They are separate rate laws that happen to move the same mineral, which is
 * how the code is meant to be used: the equilibrium solver does not compute
 * solubility products, so there is no saturation state to drive a single
 * reversible rate.
 * ============================================================================= */
#ifndef DEFINE_ABIOTIC_KINETICS_HH
#define DEFINE_ABIOTIC_KINETICS_HH

#include <vector>
#include <cstddef>
#include <cmath>
#include <algorithm>

namespace ExampleBoth {
/* WHY THIS RATE CONSTANT AND NOT A LARGER ONE
 *
 *   The solver applies these rates explicitly:  C += subsR * dt.  Nothing between
 *   this function and the concentration field stops a rate from consuming more
 *   reactant than the voxel holds, so the rate law itself must be small enough
 *   that it cannot.  The condition is
 *
 *       R * dt  <<  min(reactant concentrations)
 *
 *   With the second-order form R = k * A * B, the worst case is both reactants at
 *   their inlet value C_in, which gives  k * C_in^2 * dt  <<  C_in, i.e.
 *
 *       k  <<  1 / (C_in * dt)
 *
 *   This case has C_in = 5 mol/L and dt near 1.2e-2 s, so k must stay well below
 *   about 16 L/mol/s.  The value below leaves an order of magnitude of margin.
 *
 *   Getting this wrong does not crash: concentrations simply go negative, the run
 *   continues, and the answer is nonsense.  <diagnostics> reports it -- that is
 *   what the "went negative" line is for.
 */
    const double kp_FeS = 1.0;       // L/mol/s     precipitation (see the note above)
    const double kd_FeS = 5.0e-2;    // L/mol/s     dissolution, per mol/L left
}

void defineAbioticRxnKinetics(std::vector<double> C, std::vector<double>& subsR,
                              plb::plint mask)
{
    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    if (C.size() < 4 || subsR.size() < 4) return;
    if (mask < 2) return;

    const double Fe2 = std::max(C[0], 0.0);
    const double HS  = std::max(C[1], 0.0);

    const double R = ExampleBoth::kp_FeS * Fe2 * HS;

    subsR[0] = -R;
    subsR[1] = -R;
    subsR[3] = +R;

    for (std::size_t i = 0; i < subsR.size(); ++i)
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
}

inline void defineDissolutionRate(plb::plint phaseId, const std::vector<double>& C,
                                  double mineral, std::vector<double>& subsR,
                                  double& mineralR, plb::plint mask)
{
    (void) mask;

    mineralR = 0.0;
    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    if (C.size() < 4 || subsR.size() < 4) return;

    if (phaseId == 1) {                        // the FeS precipitate
        const double H = std::max(C[2], 0.0);
        const double rate = ExampleBoth::kd_FeS * H * std::max(mineral, 0.0);

        mineralR = -rate;
        subsR[0] = +rate;                      // Fe2+ back into the water
        subsR[1] = +rate;                      // HS- back into the water
        subsR[2] = -rate;                      // proton consumed
    }

    if (std::isnan(mineralR) || std::isinf(mineralR)) mineralR = 0.0;
    for (std::size_t i = 0; i < subsR.size(); ++i)
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
}

#endif
