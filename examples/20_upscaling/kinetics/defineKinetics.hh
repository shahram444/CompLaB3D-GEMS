/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* =============================================================================
 * defineKinetics.hh  --  EXAMPLE 19: anaerobic oxidation of methane
 *
 *     CH4 + SO4(2-)  ->  HS(-) + HCO3(-) + H2O
 *
 *     growth = mu_max * CH4/(K_CH4 + CH4) * SO4/(K_SO4 + SO4) * B    gDW/L/s
 *     CH4  consumed at growth / Y
 *     SO4  consumed at growth / Y            one sulfate per methane
 *     HS   produced at growth / Y
 *     HCO3 produced at growth / Y
 *
 * Substrate order, matching CompLaB.xml:
 *     C[0] = CH4   C[1] = SO4   C[2] = HS   C[3] = HCO3
 *
 * THIS FILE KNOWS NOTHING ABOUT THERMODYNAMICS, AND THAT IS THE POINT.  It is
 * an ordinary dual-Monod law: it reads the two substrates, ignores the two
 * products, and returns a positive rate wherever there is any methane and any
 * sulfate at all.  Left to itself it keeps the reaction running in the middle
 * of the aggregate, where sulfide and bicarbonate have built up to the point
 * that no organism could conserve energy from it.
 *
 * The gate that fixes that is in input/aom.thm, applied by the solver AFTER
 * this function returns.  Nothing here has to change to switch it on or off.
 * Compare the two runs and the difference is entirely the gate.
 *
 * ON THE RATE CONSTANT.  mu_max below is about six orders of magnitude faster
 * than a real ANME-SRB consortium, which doubles on a timescale of months
 * (Orphan et al. 2009).  A case at the measured rate would need a simulated
 * year to build the product gradient this one builds in a minute, and the
 * gradient is what the example is about.
 *
 * WHAT IS ACTUALLY BEING MATCHED is not the rate constant but the Thiele
 * modulus, r L^2 / (D C): the dimensionless group that decides whether a
 * reaction inside a diffusive aggregate is limited by supply, and therefore
 * whether an interior gradient forms at all.  In a real aggregate over months
 * that group is large.  Raising mu_max is how a short run reproduces the same
 * regime.  Every other number here is the measured one.  Do not quote this
 * rate constant.
 * ============================================================================= */
#ifndef DEFINE_KINETICS_HH
#define DEFINE_KINETICS_HH

#include <vector>
#include <cstddef>
#include <cmath>
#include <algorithm>

/* -----------------------------------------------------------------------------
 * KineticsStats -- REQUIRED BY THE SOLVER, NOT OPTIONAL
 *
 * complab.cpp calls KineticsStats::getStats() and KineticsStats::resetIteration()
 * at every output interval.  Those calls are compiled unconditionally, so ANY
 * defineKinetics.hh must declare this namespace or the program fails to link.
 * ----------------------------------------------------------------------------- */
namespace KineticsStats {
    static const double MIN_BIOMASS = 1e-12;
    static double iter_sum_dB = 0.0, iter_max_biomass = 0.0, iter_max_dB = 0.0, iter_min_DOC = 1e30;
    static long   iter_cells_with_biomass = 0, iter_cells_with_growth = 0;

    inline void resetIteration() {
        iter_sum_dB = 0; iter_max_biomass = 0; iter_max_dB = 0; iter_min_DOC = 1e30;
        iter_cells_with_biomass = 0; iter_cells_with_growth = 0;
    }
    inline void accumulate(double biomass, double donor, double dB) {
        if (biomass > MIN_BIOMASS) {
            iter_cells_with_biomass++; iter_sum_dB += dB;
            if (biomass > iter_max_biomass) iter_max_biomass = biomass;
            if (dB > iter_max_dB) iter_max_dB = dB;
            if (donor < iter_min_DOC && donor > 0) iter_min_DOC = donor;
            if (dB > 0) iter_cells_with_growth++;
        }
    }
    inline void getStats(long& cb, long& cg, double& s, double& mB, double& mdB, double& mD) {
        cb = iter_cells_with_biomass; cg = iter_cells_with_growth;
        s  = iter_sum_dB; mB = iter_max_biomass; mdB = iter_max_dB;
        mD = (iter_min_DOC < 1e20) ? iter_min_DOC : 0.0;
    }
}

namespace AOM {
    const double mu_max = 1.0e-1;   // 1/s     see the note above: NOT the measured rate
    const double K_CH4  = 1.0e-3;   // mol/L   Nauhaus et al. (2002)
    const double K_SO4  = 5.0e-4;   // mol/L   Nauhaus et al. (2002)
    const double Y      = 0.25;     // gDW per mol CH4  == 1% of the methane carbon in biomass
    const double kd     = 0.0;      // 1/s     decay off, so the gate is the only thing acting
}

inline double monod(double S, double Ks) { return (S > 0.0) ? S / (Ks + S) : 0.0; }

void defineRxnKinetics(std::vector<double> B, std::vector<double> C,
                       std::vector<double>& subsR, std::vector<double>& bioR,
                       plb::plint mask)
{
    using namespace AOM;

    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    for (std::size_t i = 0; i < bioR.size();  ++i) bioR[i]  = 0.0;
    if (C.size() < 4 || subsR.size() < 4) return;
    if (mask < 2) return;                       // no biology in solid or wall voxels

    const double CH4 = std::max(C[0], 0.0);
    const double SO4 = std::max(C[1], 0.0);

    for (std::size_t m = 0; m < B.size() && m < bioR.size(); ++m) {
        const double Bm = std::max(B[m], 0.0);
        if (Bm <= 0.0) continue;

        const double growth = mu_max * monod(CH4, K_CH4) * monod(SO4, K_SO4) * Bm;   // gDW/L/s
        const double r      = growth / Y;                                            // mol/L/s

        bioR[m]  += growth - kd * Bm;
        subsR[0] -= r;          // CH4  consumed
        subsR[1] -= r;          // SO4  consumed
        subsR[2] += r;          // HS   produced
        subsR[3] += r;          // HCO3 produced

        KineticsStats::accumulate(Bm, CH4, growth - kd * Bm);
    }

    for (std::size_t i = 0; i < subsR.size(); ++i)
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
    for (std::size_t i = 0; i < bioR.size(); ++i)
        if (std::isnan(bioR[i]) || std::isinf(bioR[i])) bioR[i] = 0.0;
}

#endif
