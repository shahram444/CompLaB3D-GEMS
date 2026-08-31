/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* =============================================================================
 * defineKinetics.hh  --  EXAMPLE: no biological reactions
 *
 * Every rate is zero. Transport only. The file exists so that this example
 * cannot accidentally compile against the shipped uranium network, which
 * indexes C[0] to C[94] and would read past the end of a short array.
 * ============================================================================= */
#ifndef DEFINE_KINETICS_HH
#define DEFINE_KINETICS_HH

#include <vector>
#include <cstddef>

void defineRxnKinetics(std::vector<double> B, std::vector<double> C,
                       std::vector<double>& subsR, std::vector<double>& bioR,
                       plb::plint mask)
{
    (void) B; (void) C; (void) mask;
    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    for (std::size_t i = 0; i < bioR.size();  ++i) bioR[i]  = 0.0;
}

#endif
