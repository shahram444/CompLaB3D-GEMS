/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

/* =============================================================================
 * defineKinetics.hh  --  EXAMPLE 12: maintenance, running alongside FBA
 *
 * Example 12 gives microbe0 <reaction_type>glpk_and_kinetics</reaction_type>.
 * Both paths run and their rates ADD. That is the point of the combined types:
 * flux balance analysis supplies the growth, and this file supplies whatever
 * the metabolic model does not describe.
 *
 * Here that is maintenance: a constant per-biomass drain on the donor that
 * yields no growth. Real cells pay it; a plain FBA objective does not include
 * it, so the model over-predicts growth at low substrate without it.
 *
 * Substrate order, matching CompLaB.xml:
 *     C[0] = S (donor)   C[1] = O (acceptor)   C[2] = P (product)
 *
 * NOTE THE SIGN DISCIPLINE. This file must not write bioR: growth belongs to
 * the FBA path, and writing it here would double count. It only draws down the
 * donor. Getting this wrong is the classic mistake with combined reaction
 * types, so it is worth being explicit about.
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

namespace ExampleMaintenance {
    const double m_S = 2.0e-6;      // mol of donor per gDW per second, no growth
}

void defineRxnKinetics(std::vector<double> B, std::vector<double> C,
                       std::vector<double>& subsR, std::vector<double>& bioR,
                       plb::plint mask)
{
    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    for (std::size_t i = 0; i < bioR.size();  ++i) bioR[i]  = 0.0;
    if (C.size() < 3 || subsR.size() < 3 || B.empty()) return;
    if (mask < 2) return;

    const double S = std::max(C[0], 0.0);
    if (S <= 0.0) return;                        // no donor, no maintenance draw

    double drain = 0.0;
    for (std::size_t m = 0; m < B.size(); ++m) drain += ExampleMaintenance::m_S * std::max(B[m], 0.0);

    // Never draw more than is present in one step. dt is not visible here, so
    // this is a coarse guard; the solver clamps properly downstream.
    subsR[0] = -std::min(drain, S);

    // bioR deliberately left at zero: growth is the FBA path's job.

    for (std::size_t i = 0; i < subsR.size(); ++i)
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
}

#endif
