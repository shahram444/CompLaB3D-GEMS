/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* =============================================================================
 * defineKinetics.hh  --  EXAMPLES 05 to 08: Monod growth
 *
 *     growth   = mu_max * S/(Ks + S) * B      gDW/L/s
 *     donor    consumed at growth / Y
 *     product  made at growth * fP
 *     decay    kd * B, first order, always on
 *
 * Substrate order, matching CompLaB.xml:
 *     C[0] = donor      C[1] = product
 *
 * This file is written for ANY number of microbes: it loops over B.size() and
 * takes each organism's parameters from the tables below. Examples 05, 06 and
 * 07 have one microbe and use index 0; example 08 has two and uses both.
 *
 * WHY THE PARAMETERS LIVE HERE AND NOT IN THE XML
 *   <half_saturation_constants> and <maximum_uptake_flux> in CompLaB.xml feed
 *   the FBA and surrogate paths, which build their own uptake bounds. A
 *   hand-written kinetics file is not obliged to read them, and this one does
 *   not: everything it needs is in the table below, where you can see it next
 *   to the equation that uses it. Change a number here and rebuild.
 *
 * WHAT TO EXPECT
 *   Biomass grows where the donor reaches it, the donor is drawn down around
 *   the colony, and the product plume spreads from it. Because decay is first
 *   order and growth saturates, biomass approaches a steady value rather than
 *   growing without bound.
 * ============================================================================= */
#ifndef DEFINE_KINETICS_HH
#define DEFINE_KINETICS_HH

#include <vector>
#include <cstddef>
#include <cmath>
#include <algorithm>

namespace ExampleMonod {
    // one entry per microbe, in CompLaB.xml order
    const double mu_max[2] = { 2.0e-4, 1.0e-4 };   // 1/s     max specific growth rate
    const double Ks    [2] = { 5.0e-2, 5.0e-3 };   // mol/L   half saturation on the donor
    const double Y     [2] = { 0.4,    0.4    };   // gDW per mol of donor
    const double fP    [2] = { 2.0,    2.0    };   // mol of product per gDW made
    const double kd    [2] = { 1.0e-6, 1.0e-6 };   // 1/s     first order decay
    const int    NMICROBE   = 2;
}

inline double monod(double S, double Ks) { return S / (Ks + S); }

void defineRxnKinetics(std::vector<double> B, std::vector<double> C,
                       std::vector<double>& subsR, std::vector<double>& bioR,
                       plb::plint mask)
{
    using namespace ExampleMonod;

    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    for (std::size_t i = 0; i < bioR.size();  ++i) bioR[i]  = 0.0;
    if (C.size() < 2 || subsR.size() < 2) return;
    if (mask < 2) return;                       // no biology in solid or wall voxels

    const double S = std::max(C[0], 0.0);

    for (std::size_t m = 0; m < B.size() && m < bioR.size(); ++m) {
        const int p = (m < (std::size_t) NMICROBE) ? (int) m : NMICROBE - 1;
        const double Bm = std::max(B[m], 0.0);
        if (Bm <= 0.0) continue;

        const double growth = mu_max[p] * monod(S, Ks[p]) * Bm;   // gDW/L/s

        bioR[m]  += growth - kd[p] * Bm;
        subsR[0] -= growth / Y[p];                                // donor consumed
        subsR[1] += growth * fP[p];                               // product released
    }

    for (std::size_t i = 0; i < subsR.size(); ++i)
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
    for (std::size_t i = 0; i < bioR.size(); ++i)
        if (std::isnan(bioR[i]) || std::isinf(bioR[i])) bioR[i] = 0.0;
}

#endif
