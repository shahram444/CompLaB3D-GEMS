/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* =============================================================================
 * defineAbioticKinetics.hh  --  EXAMPLE 14: calcite dissolution
 *
 *     CaCO3(s) + H+  ->  Ca2+ + HCO3-
 *
 *     rate = k * [H+] * mineral       first order in acid and in how much
 *                                     mineral is left, the usual crude stand
 *                                     in for reactive surface area
 *
 * Substrate order, matching CompLaB.xml:
 *     C[0] = H        C[1] = Ca       C[2] = calcite  (immobile, the mineral)
 *
 * NOTHING PRECIPITATES IN THIS EXAMPLE. defineAbioticRxnKinetics below is
 * empty on purpose, so the only thing moving mineral is the dissolution hook.
 * That makes the mass balance unambiguous: every mole of calcite that leaves a
 * grain must appear as a mole of Ca2+ in the water.
 *
 * WHAT TO EXPECT
 *   Acid enters from the left, eats into the upstream faces of the four
 *   grains, and Ca2+ leaves on the right. As each grain thins, its outer
 *   voxels drop below <reopen_fraction> of full and reopen to flow, so the
 *   porosity in the log rises step by step. That is the mirror image of what
 *   example 13 does.
 *
 * THE ARGUMENTS, briefly (the full contract is in the shipped file):
 *   phaseId  1-based, in <phaseN> order.  C  the water TOUCHING the mineral,
 *   averaged over the open neighbours.  mineral  mol/L left in this voxel.
 *   Return mineralR NEGATIVE, and the released solutes in subsR, positive.
 *   Do not write the mineral's own subsR entry; the solver does that from
 *   mineralR. Do not add your own step limiter; the solver caps removal at
 *   half the mineral present and scales everything you return by the same
 *   factor, so the element balance survives the cap exactly.
 * ============================================================================= */
#ifndef DEFINE_ABIOTIC_KINETICS_HH
#define DEFINE_ABIOTIC_KINETICS_HH

#include <vector>
#include <cstddef>
#include <cmath>
#include <algorithm>

namespace ExampleDissol {
    const double k_calcite = 1.0e-2;   // L/mol/s per mol/L of remaining mineral
}

void defineAbioticRxnKinetics(std::vector<double> C, std::vector<double>& subsR,
                              plb::plint mask)
{
    (void) C; (void) mask;
    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;   // nothing precipitates here
}

inline void defineDissolutionRate(plb::plint phaseId, const std::vector<double>& C,
                                  double mineral, std::vector<double>& subsR,
                                  double& mineralR, plb::plint mask)
{
    (void) mask;

    mineralR = 0.0;
    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    if (C.size() < 3 || subsR.size() < 3) return;

    if (phaseId == 1) {                       // calcite: the only declared phase
        const double H = std::max(C[0], 0.0);
        const double rate = ExampleDissol::k_calcite * H * std::max(mineral, 0.0);  // mol/L/s

        mineralR = -rate;                     // the mineral disappears
        subsR[1] = +rate;                     // one Ca2+ per CaCO3
        subsR[0] = -rate;                     // one proton consumed
        // subsR[2] is the mineral itself: left alone on purpose
    }

    if (std::isnan(mineralR) || std::isinf(mineralR)) mineralR = 0.0;
    for (std::size_t i = 0; i < subsR.size(); ++i)
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
}

#endif
