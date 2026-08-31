/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* =============================================================================
 * defineAbioticKinetics.hh  --  EXAMPLE: no abiotic reactions
 *
 * Every rate is zero, and defineDissolutionRate does nothing. Transport only.
 * ============================================================================= */
#ifndef DEFINE_ABIOTIC_KINETICS_HH
#define DEFINE_ABIOTIC_KINETICS_HH

#include <vector>
#include <cstddef>
#include <cmath>

void defineAbioticRxnKinetics(std::vector<double> C, std::vector<double>& subsR,
                              plb::plint mask)
{
    (void) C; (void) mask;
    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
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
