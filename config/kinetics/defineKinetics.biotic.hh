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

/* -----------------------------------------------------------------------------
 * KineticsStats -- REQUIRED BY THE SOLVER, NOT OPTIONAL
 *
 * complab.cpp calls KineticsStats::getStats() and KineticsStats::resetIteration()
 * at every output interval, to print the per-interval kinetics summary. Those
 * calls are compiled unconditionally, so ANY defineKinetics.hh must declare this
 * namespace or the whole program fails to link -- even a case with no kinetics
 * at all, where the calls are guarded at run time by `if (kns_count > 0)`.
 *
 * The shipped defineKinetics.hh in the project root has it. Generated example
 * headers must have it too; it is reproduced here rather than #included so that
 * each example header stays self-contained and can be dropped into the source
 * tree on its own.
 *
 * A rate law that never calls accumulate() simply reports zeros, which is
 * correct for a case that does no biotic kinetics.
 * ----------------------------------------------------------------------------- */
namespace KineticsStats {
    static const double MIN_BIOMASS = 1e-12;    // below this a voxel is counted as empty
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

        KineticsStats::accumulate(Bm, S, growth - kd[p] * Bm);    // feeds the per-interval report
    }

    for (std::size_t i = 0; i < subsR.size(); ++i)
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
    for (std::size_t i = 0; i < bioR.size(); ++i)
        if (std::isnan(bioR[i]) || std::isinf(bioR[i])) bioR[i] = 0.0;
}

#endif
