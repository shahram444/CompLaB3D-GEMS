"""The per-example kinetics headers.

CompLaB compiles its rate laws in: defineKinetics.hh and
defineAbioticKinetics.hh are #included, not read at run time. Two cases with
different chemistry therefore need different headers, which is why every
example folder carries its own copy and why runAllExamples.sh rebuilds between
cases. The 2D suite worked the same way; each example/scalability case shipped
its own defineKinetics.hh.

The shipped headers in CompLB3D implement the 95-substrate uranium network from
Zhao et al. Those would read past the end of a 2-substrate array, so every
example here replaces them with something the size of its own chemistry.

Units, throughout, matching the solver:
    concentrations   the units of <initial_concentration>, per litre
    biomass          gDW/L
    every rate       PER SECOND
A rate constant published per day must be divided by 86400. That one mistake
causes more dead runs than anything else in this code.
"""

LICENCE = """/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/
"""

# ---------------------------------------------------------------------------
# Inert: no reactions at all. Used by the transport-only cases so that the
# shipped uranium network cannot fire by accident.
# ---------------------------------------------------------------------------
KINETICS_INERT = LICENCE + """
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
"""

ABIOTIC_INERT = LICENCE + """
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
"""

# ---------------------------------------------------------------------------
# Example 03: one irreversible second-order abiotic reaction, A + B -> C
# ---------------------------------------------------------------------------
ABIOTIC_ABC = LICENCE + """
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
 *   One A plus one B makes exactly one C, so at every step
 *       dA + dC = 0   and   dB + dC = 0
 *   Any drift in that sum is a bug, not chemistry. runAllExamples.sh checks it.
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
"""

# ---------------------------------------------------------------------------
# Examples 05 to 08 and 12: Monod growth on one donor, one waste product
# ---------------------------------------------------------------------------
KINETICS_MONOD = LICENCE + """
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
"""

# ---------------------------------------------------------------------------
# Example 12: kinetics that run ALONGSIDE flux balance analysis
# ---------------------------------------------------------------------------
KINETICS_MAINTENANCE = LICENCE + """
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
"""

# ---------------------------------------------------------------------------
# Example 13: mineral precipitation, Fe2+ + HS- -> FeS(s)
# ---------------------------------------------------------------------------
ABIOTIC_PRECIP = LICENCE + """
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
"""

# ---------------------------------------------------------------------------
# Example 14: mineral dissolution, calcite grains in acid
# ---------------------------------------------------------------------------
ABIOTIC_DISSOL = LICENCE + """
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
"""

# ---------------------------------------------------------------------------
# Example 15: both directions at once
# ---------------------------------------------------------------------------
ABIOTIC_BOTH = LICENCE + """
/* =============================================================================
 * defineAbioticKinetics.hh  --  EXAMPLE 15: precipitation AND dissolution
 *
 * Both hooks are implemented, on the same mineral, so a voxel can seal and
 * later reopen:
 *
 *     precipitation   Fe2+ + HS-  ->  FeS(s)          rate = kp * [Fe2+] * [HS-]
 *     dissolution     FeS(s) + H+ ->  Fe2+ + HS-      rate = kd * [H+] * mineral
 *
 * Substrate order, matching CompLaB.xml:
 *     C[0] = Fe2   C[1] = HS   C[2] = H   C[3] = FeS  (immobile, the mineral)
 *
 * WHY <reopen_fraction> MATTERS HERE AND NOWHERE ELSE
 *   With both directions active a voxel can sit near the sealing threshold. A
 *   single threshold would make it flip open and shut every update interval,
 *   and every flip re-solves the flow. reopen_fraction 0.9 puts a gap between
 *   sealing at full_density and reopening below 0.9 x full_density, so it
 *   settles instead of chattering. Watch the porosity in the log: it should
 *   move in steps and then hold, not oscillate.
 *
 * THE TWO HOOKS ARE NOT INVERSES. Precipitation is second order in the
 * dissolved ions; dissolution is first order in acid and in remaining mineral.
 * They are separate rate laws that happen to move the same mineral, which is
 * how the code is meant to be used: the equilibrium solver does not compute
 * solubility products, so there is no saturation state to drive a single
 * reversible rate.
 * ============================================================================= */
#ifndef DEFINE_ABIOTIC_KINETICS_HH
#define DEFINE_ABIOTIC_KINETICS_HH

#include <vector>
#include <cstddef>
#include <cmath>
#include <algorithm>

namespace ExampleBoth {
    const double kp_FeS = 2.0e+1;    // L/mol/s     precipitation
    const double kd_FeS = 5.0e-2;    // L/mol/s     dissolution, per mol/L left
}

void defineAbioticRxnKinetics(std::vector<double> C, std::vector<double>& subsR,
                              plb::plint mask)
{
    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    if (C.size() < 4 || subsR.size() < 4) return;
    if (mask < 2) return;

    const double Fe2 = std::max(C[0], 0.0);
    const double HS  = std::max(C[1], 0.0);

    const double R = ExampleBoth::kp_FeS * Fe2 * HS;

    subsR[0] = -R;
    subsR[1] = -R;
    subsR[3] = +R;

    for (std::size_t i = 0; i < subsR.size(); ++i)
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
}

inline void defineDissolutionRate(plb::plint phaseId, const std::vector<double>& C,
                                  double mineral, std::vector<double>& subsR,
                                  double& mineralR, plb::plint mask)
{
    (void) mask;

    mineralR = 0.0;
    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;
    if (C.size() < 4 || subsR.size() < 4) return;

    if (phaseId == 1) {                        // the FeS precipitate
        const double H = std::max(C[2], 0.0);
        const double rate = ExampleBoth::kd_FeS * H * std::max(mineral, 0.0);

        mineralR = -rate;
        subsR[0] = +rate;                      // Fe2+ back into the water
        subsR[1] = +rate;                      // HS- back into the water
        subsR[2] = -rate;                      // proton consumed
    }

    if (std::isnan(mineralR) || std::isinf(mineralR)) mineralR = 0.0;
    for (std::size_t i = 0; i < subsR.size(); ++i)
        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;
}

#endif
"""
