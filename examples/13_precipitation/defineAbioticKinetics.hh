/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* =============================================================================
 * defineAbioticKinetics.hh  --  EXAMPLE 13: FeS precipitation
 *
 *     Fe2+ + HS-  ->  FeS(s)        rate = k * [Fe2+] * [HS-]
 *
 * Substrate order, matching CompLaB.xml:
 *     C[0] = Fe2      C[1] = HS      C[2] = FeS   (immobile, the mineral)
 *
 * FeS is declared <immobile>true</immobile> in the XML, so it never advects or
 * diffuses. It accumulates exactly where it forms, which is what lets
 * <precipitation> convert a voxel to solid once it reaches max_precipRho.
 *
 * WHAT TO EXPECT
 *   Fe2+ enters from the left, HS- from the right. FeS builds up where they
 *   meet, the middle voxels seal, and the flow is forced around them. Watch
 *   the porosity line in the log fall step by step.
 *
 * MASS BALANCE: one Fe2+ plus one HS- makes exactly one FeS.
 * ============================================================================= */
#ifndef DEFINE_ABIOTIC_KINETICS_HH
#define DEFINE_ABIOTIC_KINETICS_HH

#include <vector>
#include <cstddef>
#include <cmath>
#include <algorithm>

namespace ExamplePrecip {
    const double k_FeS = 2.0e+1;     // L/mol/s
}

void defineAbioticRxnKinetics(std::vector<double> C, std::vector<double>& subsR,
                              plb::plint mask)
{
    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    if (C.size() < 3 || subsR.size() < 3) return;
    if (mask < 2) return;

    const double Fe2 = std::max(C[0], 0.0);
    const double HS  = std::max(C[1], 0.0);

    const double R = ExamplePrecip::k_FeS * Fe2 * HS;    // mol/L/s

    subsR[0] = -R;
    subsR[1] = -R;
    subsR[2] = +R;                                       // the mineral, immobile

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
