/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* =============================================================================
 * defineAbioticKinetics.hh  --  EXAMPLE 03: A + B  ->  C
 *
 *     rate = k * [A] * [B]        second order, irreversible
 *
 * Substrate order, matching CompLaB.xml:
 *     C[0] = A      C[1] = B      C[2] = C
 *
 * WHAT TO EXPECT
 *   A enters from the left boundary, B from the right, and C forms in the
 *   middle where the two overlap. With no flow the reaction front sits at the
 *   centre; with flow it is pushed downstream. That is the whole point of the
 *   case: the peak in C tells you the coupling between transport and reaction
 *   is working.
 *
 * MASS BALANCE
 *   One A plus one B makes exactly one C, so at every step the reaction
 *   INCREMENTS satisfy
 *       dA + dC = 0   and   dB + dC = 0
 *   Any drift in that sum is a bug, not chemistry.
 *
 *   [v1.3] Those sums are over the increments, meaning the amounts this law
 *   adds to and removes from each field in one time step, and NOT over the
 *   reported totals. Earlier text here said the totals must satisfy them.
 *   They cannot: CompLaB.xml starts A and B at 0 and holds each at 1.0 on one
 *   face with a Dirichlet condition, that is, a boundary that fixes the
 *   concentration, so both are supplied from outside the domain and all three
 *   totals rise. The checkable quantities are that C rises, C being the only
 *   closed species, and that the increments balance. Naming A or B in a
 *   <conserve> sum would guarantee a mass-balance FAIL that means nothing.
 *   postprocess.py prints the change in each field after a run, which is where
 *   you read it off.
 * ============================================================================= */
#ifndef DEFINE_ABIOTIC_KINETICS_HH
#define DEFINE_ABIOTIC_KINETICS_HH

#include <vector>
#include <cstddef>
#include <cmath>
#include <algorithm>

namespace ExampleAbiotic {
    const double k_AB = 5.0e-2;      // L/mol/s
}

void defineAbioticRxnKinetics(std::vector<double> C, std::vector<double>& subsR,
                              plb::plint mask)
{
    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    if (C.size() < 3 || subsR.size() < 3) return;
    if (mask < 2) return;                       // nothing happens in solid or wall voxels

    const double A = std::max(C[0], 0.0);       // floor tiny negatives from the LB solver
    const double B = std::max(C[1], 0.0);

    const double R = ExampleAbiotic::k_AB * A * B;   // mol/L/s

    subsR[0] = -R;
    subsR[1] = -R;
    subsR[2] = +R;

    for (std::size_t i = 0; i < subsR.size(); ++i)
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
}

inline void defineDissolutionRate(plb::plint phaseId, const std::vector<double>& C,
                                  double mineral, std::vector<double>& subsR,
                                  double& mineralR, plb::plint mask)
{
    (void) phaseId; (void) C; (void) mineral; (void) mask;
    mineralR = 0.0;
    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
}

#endif
