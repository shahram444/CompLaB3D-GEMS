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

void defineRxnKinetics(std::vector<double> B, std::vector<double> C,
                       std::vector<double>& subsR, std::vector<double>& bioR,
                       plb::plint mask)
{
    (void) B; (void) C; (void) mask;
    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    for (std::size_t i = 0; i < bioR.size();  ++i) bioR[i]  = 0.0;
}

#endif
