/* This file is a part of the CompLaB program.
 *
 * The CompLaB3D software (3-D pore-scale extension) is developed since 2024
 * by the University of Georgia (United States, Meile Lab, Department of
 * Marine Sciences). The original 2-D CompLaB v1.0 was a collaboration of
 * the University of Georgia and Chungnam National University (South Korea).
 *
 * The CompLaB3D extension (3-D pore-scale model with biotic/abiotic
 * kinetics, Anderson-Accelerated Newton-Raphson PCF equilibrium chemistry,
 * and the CompLaB Studio GUI) is developed by Shahram Asgari and
 * Christof Meile (Meile Lab) at the University of Georgia. The 2-D
 * predecessor (CompLaB v1.0) was developed by Heewon Jung et al.
 *
 * Contact:
 * Shahram Asgari
 * Department of Marine Sciences (Meile Lab)
 * University of Georgia
 * Athens, GA 30602, USA
 * shahram.asgari@uga.edu
 *
 * Christof Meile
 * Department of Marine Sciences
 * University of Georgia
 * Athens, GA 30602, USA
 *
 * The most recent release of CompLaB can be downloaded at
 * https://bitbucket.org/MeileLab/complab/downloads/
 *
 * CompLaB is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * ┌────────────────────────────────────────────────────────────────────────┐
 * │  GUARD: These two lines prevent the compiler from reading this file   │
 * │  more than once. "#ifndef" means "if not already defined," and        │
 * │  "#define" marks it as defined. Think of it as a "read once" stamp.   │
 * └────────────────────────────────────────────────────────────────────────┘
 */
#ifndef COMPLAB3D_PROCESSORS_PART4_HH
#define COMPLAB3D_PROCESSORS_PART4_HH

/*
 * ┌────────────────────────────────────────────────────────────────────────┐
 * │  INCLUDES: These lines import standard C++ math and data-structure    │
 * │  libraries so we can use things like sqrt(), pow(), vectors, etc.     │
 * └────────────────────────────────────────────────────────────────────────┘
 */
#include <cmath>       /* math functions: sqrt, pow, log10, abs, etc.           */
#include <limits>      /* provides numeric_limits (e.g., largest possible double) */
#include <vector>      /* dynamic arrays: std::vector<T> (resizable lists)       */
#include <string>      /* text strings: std::string                              */
#include <algorithm>   /* find, min, max, swap — utility algorithms              */
#include <iomanip>     /* formatting for printed output (setw, setprecision)     */

/*
 * "using namespace plb;" means we can write "plint" instead of "plb::plint".
 * "plb" is the Palabos library namespace — all Palabos types live there.
 *
 * "typedef double T;" creates a shorthand: anywhere we write "T" it means
 * "double" (a 64-bit floating-point number, ~15 digits of precision).
 */
using namespace plb;
typedef double T;


/* ╔════════════════════════════════════════════════════════════════════════════╗
 * ║                                                                          ║
 * ║  FILE OVERVIEW: EQUILIBRIUM CHEMISTRY SOLVER                             ║
 * ║  ═══════════════════════════════════════════════════════════════          ║
 * ║                                                                          ║
 * ║  WHAT THIS FILE DOES:                                                    ║
 * ║  This file solves "fast equilibrium reactions" — chemical reactions that  ║
 * ║  happen nearly instantly compared to flow and biological timescales.      ║
 * ║  Example: the carbonate system (CO2 + H2O ↔ H2CO3 ↔ HCO3- + H+ ↔ ...) ║
 * ║                                                                          ║
 * ║  THE PROBLEM:                                                            ║
 * ║  Given the total amount of each chemical element in a tiny volume (voxel)║
 * ║  of pore space, figure out how that total is distributed among the       ║
 * ║  different chemical forms (species) at equilibrium.                      ║
 * ║                                                                          ║
 * ║  Example: if total carbon = 0.001 mol/L, how much is H2CO3, how much    ║
 * ║  is HCO3-, and how much is CO3-- ?                                       ║
 * ║                                                                          ║
 * ║  THE METHOD:                                                             ║
 * ║  We use the "PCF" method (Positive Continuous Fraction) accelerated by   ║
 * ║  "Anderson Acceleration" (AA). This combination was chosen over          ║
 * ║  Newton-Raphson because it is more robust and converges more reliably.   ║
 * ║  See the reserve file (*_NR_reserve.hh) for the Newton-Raphson version.  ║
 * ║                                                                          ║
 * ║  WHERE IT FITS IN THE SIMULATION:                                        ║
 * ║  At each timestep, CompLaB runs three phases in order:                   ║
 * ║    1. Transport   — move dissolved species through pore space (LBM)      ║
 * ║    2. Kinetics    — run slow biological/chemical reactions               ║
 * ║    3. Equilibrium — THIS solver — instantly balance fast reactions        ║
 * ║                                                                          ║
 * ║  REFERENCES:                                                             ║
 * ║  [1] Anderson (1965) - Original Anderson Acceleration method             ║
 * ║  [2] Walker & Ni (2011) - Convergence analysis, condition monitoring     ║
 * ║  [3] Carrayrou et al. (2002) - PCF method for equilibrium chemistry      ║
 * ║  [4] Awada et al. (2025) - Anderson + PCF combination, 99.5% success    ║
 * ║                                                                          ║
 * ╚════════════════════════════════════════════════════════════════════════════╝
 *
 *
 *  KEY VOCABULARY (for readers not familiar with C++ or Palabos):
 *  ──────────────────────────────────────────────────────────────
 *  "template<typename T>"  = This code works with any number type (float, double).
 *                            We always use T = double in CompLaB.
 *
 *  "class"                 = A blueprint that bundles data and functions together.
 *                            Think of it like a "module" containing everything
 *                            the equilibrium solver needs.
 *
 *  "std::vector<T>"        = A resizable list of numbers. Like a Python list.
 *                            Example: {0.001, 0.005, 1e-7} is a vector of 3 numbers.
 *
 *  "size_t"                = A non-negative integer type used for counting and indexing.
 *
 *  "plint"                 = Palabos integer type (same as long int).
 *
 *  "const"                 = "This value will not be changed by this function."
 *
 *  "mutable"               = Exception to const: this variable CAN be changed
 *                            even inside a function marked const.
 *
 *  "&" (reference)         = "Pass the actual variable, not a copy." Saves memory.
 *
 *  "static constexpr"      = A constant that is set once and never changes.
 *                            The compiler bakes it directly into the program.
 *
 *  "BlockLattice3D"        = A Palabos object representing a 3D grid of numbers.
 *                            Each grid cell holds 7 "populations" (for D3Q7 lattice)
 *                            whose sum equals the concentration at that cell.
 *
 *  "computeDensity()"      = A Palabos function that sums the 7 populations at a
 *                            grid cell to get the concentration value.
 *
 *  "getPopulations(g)"     = Read the 7 population values into array g.
 *
 *  "setPopulations(g)"     = Write the 7 population values from array g back.
 *
 *  "Box3D domain"          = A rectangular 3D region defined by x0..x1, y0..y1, z0..z1.
 *
 *  "Dot3D"                 = A 3D integer coordinate (x, y, z).
 *
 *  "computeRelativeDisplacement" = Palabos function to align two lattices that
 *                            may be stored at different memory offsets.
 */


/* ╔════════════════════════════════════════════════════════════════════════════╗
 * ║                                                                          ║
 * ║  CLASS: EquilibriumChemistry                                             ║
 * ║  ═══════════════════════════                                             ║
 * ║                                                                          ║
 * ║  This class contains everything needed to solve the equilibrium problem: ║
 * ║  - The list of species and components                                    ║
 * ║  - The stoichiometry matrix and equilibrium constants (logK)             ║
 * ║  - The PCF fixed-point iteration                                         ║
 * ║  - The Anderson Acceleration that speeds up PCF                          ║
 * ║  - Helper math functions (vector operations, QR decomposition)           ║
 * ║                                                                          ║
 * ║  TERMINOLOGY:                                                            ║
 * ║  - "Species" = all chemical forms (H2CO3, HCO3-, CO3--, H+, ...)        ║
 * ║  - "Components" = the independent building blocks from which all         ║
 * ║    species are formed (e.g., HCO3- and H+). These are like the          ║
 * ║    "basis vectors" of the chemical system.                               ║
 * ║  - "logK" = log10 of the equilibrium constant for forming each species  ║
 * ║    from the components.                                                  ║
 * ║  - "Stoichiometry matrix S" = a table where S[i][j] tells how many      ║
 * ║    units of component j are needed to form one unit of species i.        ║
 * ║                                                                          ║
 * ╚════════════════════════════════════════════════════════════════════════════╝
 */
template<typename T>
class EquilibriumChemistry
{
public:
    /*
     * ┌──────────────────────────────────────────────────────────────────┐
     * │  CONSTANTS                                                      │
     * │                                                                  │
     * │  These are fixed numbers that set physical and numerical limits. │
     * │  "static constexpr" means they are baked into the code at        │
     * │  compile time and never change during the simulation.            │
     * └──────────────────────────────────────────────────────────────────┘
     */

    /* MIN_CONC: The smallest concentration we allow (mol/L).
     * Anything below 1e-30 is essentially zero — computers can't
     * distinguish it from zero anyway with 64-bit precision.
     */
    static constexpr T MIN_CONC = 1e-30;

    /* MAX_CONC: The largest concentration we allow (10 mol/L).
     * Real geochemical concentrations rarely exceed 1 mol/L,
     * so 10 is a generous upper bound to catch runaway values.
     */
    static constexpr T MAX_CONC = 10.0;

    /* MIN_LOG_C and MAX_LOG_C: log10 versions of the above.
     * log10(1e-30) = -30  and  log10(10) = 1.
     * The solver works in "log-space" (log10 of concentrations)
     * because concentrations span many orders of magnitude
     * (e.g., [H+] at pH 2 is 0.01, at pH 12 is 0.000000000001).
     * Working in log-space makes these just "2" and "12".
     */
    static constexpr T MIN_LOG_C = -30.0;
    static constexpr T MAX_LOG_C = 1.0;

    /* ANDERSON ACCELERATION PARAMETERS:
     *
     * DEFAULT_ANDERSON_DEPTH = 4:
     *   How many past iterations Anderson Acceleration "remembers."
     *   More memory = potentially faster convergence, but also more
     *   risk of numerical instability. Values 3-5 work well for
     *   geochemical systems.
     *
     * DEFAULT_CONDITION_TOL = 1e10:
     *   The maximum allowed "condition number" of the AA system.
     *   The condition number measures how sensitive the solution is
     *   to small perturbations. If it exceeds this threshold, AA
     *   automatically reduces its memory depth to stay stable.
     *   Think of it as a "safety valve."
     *
     * DEFAULT_BETA = 1.0:
     *   A damping/mixing parameter. At 1.0, we take the full Anderson
     *   step. Values < 1.0 would slow convergence but add stability.
     *   1.0 works well in practice.
     */
    static constexpr plint DEFAULT_ANDERSON_DEPTH = 4;
    static constexpr T DEFAULT_CONDITION_TOL = 1e10;
    static constexpr T DEFAULT_BETA = 1.0;

    /*
     * ┌──────────────────────────────────────────────────────────────────┐
     * │  CONSTRUCTOR                                                     │
     * │                                                                  │
     * │  This runs automatically when an EquilibriumChemistry object     │
     * │  is created. It sets all the default values. The actual species, │
     * │  components, and reactions are set later via the "set" functions. │
     * └──────────────────────────────────────────────────────────────────┘
     */
    EquilibriumChemistry()
        : maxIterations(200),       /* Try up to 200 iterations before giving up         */
          tolerance(1e-8),          /* Consider converged if residual < 1e-8              */
          lastConverged(true),      /* Track whether the last solve succeeded             */
          lastIterations(0),        /* Track how many iterations the last solve took      */
          lastResidual(0.0),        /* Track the final residual of the last solve         */
          verbose(false),           /* If true, print detailed iteration-by-iteration info*/
          totalSolves(0),           /* Counter: how many times we've been called          */
          totalConverged(0),        /* Counter: how many times we converged successfully  */
          totalDiverged(0),         /* Counter: how many times we failed to converge      */
          andersonDepth(DEFAULT_ANDERSON_DEPTH),
          conditionTol(DEFAULT_CONDITION_TOL),
          beta(DEFAULT_BETA)
    {}

    /*
     * ┌──────────────────────────────────────────────────────────────────┐
     * │  SETTERS: Functions to configure the solver                      │
     * │                                                                  │
     * │  These are called during setup (before the simulation starts)    │
     * │  to tell the solver what species, components, and reactions      │
     * │  exist in this particular chemical system.                       │
     * └──────────────────────────────────────────────────────────────────┘
     */

    /* setSpeciesNames: Provide the list of all chemical species.
     * Example: {"H2CO3", "HCO3-", "H+", "CO3--"}
     * We skip "H2O" because water is the solvent, not a variable.
     */
    void setSpeciesNames(const std::vector<std::string>& names) {
        species_names.clear();     /* Empty the current list */
        for (const auto& name : names) {   /* Loop through each name */
            if (name != "H2O" && name != "h2o") {
                species_names.push_back(name);   /* Add it to our list */
            }
        }
    }

    /* setComponentNames: Provide the list of "master" components.
     * Example: {"HCO3-", "H+"} for the carbonate system.
     * Components are the independent species from which all others
     * can be expressed via equilibrium reactions.
     */
    void setComponentNames(const std::vector<std::string>& names) {
        component_names = names;
    }

    /* setLogK: Provide the log10 of equilibrium formation constants.
     * One value per species. For components themselves, logK = 0.
     * Example: {6.35, 0.0, 0.0, -10.33} meaning:
     *   H2CO3 has logK = 6.35  (it's 10^6.35 times the product of components)
     *   HCO3- has logK = 0     (it IS a component)
     *   H+    has logK = 0     (it IS a component)
     *   CO3-- has logK = -10.33
     */
    void setLogK(const std::vector<T>& logK) {
        logK_values = logK;
    }

    /* setStoichiometryMatrix: Provide the stoichiometry matrix S.
     * S[i][j] = how many units of component j are in species i.
     *
     * For the carbonate system with components = {HCO3-, H+}:
     *   S[H2CO3]  = [1,  1]   (H2CO3 = K * [HCO3-]^1 * [H+]^1)
     *   S[HCO3-]  = [1,  0]   (HCO3- = [HCO3-]^1, it's a component)
     *   S[H+]     = [0,  1]   (H+ = [H+]^1, it's a component)
     *   S[CO3--]  = [1, -1]   (CO3-- = K * [HCO3-]^1 * [H+]^(-1))
     *               ↑   ↑
     *             HCO3-  H+
     *
     * A negative stoichiometry like -1 for H+ in CO3-- means H+ appears
     * in the denominator: CO3-- = K * [HCO3-] / [H+].
     */
    void setStoichiometryMatrix(const std::vector<std::vector<T>>& S) {
        stoich_matrix = S;
    }

    /* Simple setters for solver parameters. */
    void setMaxIterations(plint maxIter) { maxIterations = maxIter; }
    void setTolerance(T tol) { tolerance = tol; }
    void setVerbose(bool v) { verbose = v; }
    void setAndersonDepth(plint depth) { andersonDepth = std::max(plint(1), depth); }
    void setConditionTolerance(T condTol) { conditionTol = condTol; }
    void setBeta(T betaVal) { beta = betaVal; }

    /*
     * ┌──────────────────────────────────────────────────────────────────┐
     * │  GETTERS: Functions to read solver state and results             │
     * │  These don't change anything — they just return information.     │
     * └──────────────────────────────────────────────────────────────────┘
     */
    plint getMaxIterations() const { return maxIterations; }
    T getTolerance() const { return tolerance; }
    std::vector<T> getLogK() const { return logK_values; }
    bool didConverge() const { return lastConverged; }
    std::vector<std::string> getSpeciesNames() const { return species_names; }
    std::vector<std::string> getComponentNames() const { return component_names; }
    plint getNumSpecies() const { return species_names.size(); }
    plint getNumComponents() const { return component_names.size(); }
    plint getLastIterations() const { return lastIterations; }
    T getLastResidual() const { return lastResidual; }
    plint getAndersonDepth() const { return andersonDepth; }
    plint getTotalSolves() const { return totalSolves; }
    plint getTotalConverged() const { return totalConverged; }
    plint getTotalDiverged() const { return totalDiverged; }

    /* resetStatistics: Zero out all the counters. */
    void resetStatistics() {
        totalSolves = totalConverged = totalDiverged = 0;
    }

    /*
     * ┌──────────────────────────────────────────────────────────────────┐
     * │  isEquilibriumSpecies                                            │
     * │                                                                  │
     * │  Checks whether species #i participates in equilibrium at all.   │
     * │  A species participates if it has at least one non-zero entry    │
     * │  in its row of the stoichiometry matrix.                         │
     * │                                                                  │
     * │  Species with all-zero rows (like biomass or inert tracers)      │
     * │  are "non-equilibrium" — this solver won't touch them.           │
     * └──────────────────────────────────────────────────────────────────┘
     */
    bool isEquilibriumSpecies(size_t species_idx) const {
        if (species_idx >= stoich_matrix.size()) return false;

        /* Loop through this species' row in the stoichiometry matrix.
         * If ANY entry has absolute value > 1e-10, this species participates. */
        for (size_t j = 0; j < stoich_matrix[species_idx].size(); ++j) {
            if (std::abs(stoich_matrix[species_idx][j]) > 1e-10) return true;
        }
        return false;   /* All entries are essentially zero — not an equilibrium species */
    }

    /* setStoichiometryRow: Set one row of the stoichiometry matrix.
     * Useful for setting up individual species after the matrix is created. */
    void setStoichiometryRow(size_t species_idx, const std::vector<T>& row) {
        if (species_idx < stoich_matrix.size() && row.size() == component_names.size()) {
            stoich_matrix[species_idx] = row;
        }
    }


    /* ╔════════════════════════════════════════════════════════════════════╗
     * ║                                                                    ║
     * ║  SECTION 1: MASS ACTION LAW                                        ║
     * ║  ══════════════════════════                                        ║
     * ║                                                                    ║
     * ║  The mass action law is the fundamental equation of equilibrium    ║
     * ║  chemistry. It says:                                               ║
     * ║                                                                    ║
     * ║    [species_i] = K_i × [comp_1]^S[i][1] × [comp_2]^S[i][2] × ... ║
     * ║                                                                    ║
     * ║  In log10 form (which we use for numerical reasons):               ║
     * ║                                                                    ║
     * ║    log10([species_i]) = logK[i] + S[i][1]*logC[1] + S[i][2]*logC[2] + ...
     * ║                                                                    ║
     * ║  This is just the same equation but with log10 on both sides.      ║
     * ║  It turns multiplication into addition, which is simpler and       ║
     * ║  more numerically stable.                                          ║
     * ║                                                                    ║
     * ╚════════════════════════════════════════════════════════════════════╝
     */

    /*
     * ┌──────────────────────────────────────────────────────────────────┐
     * │  calc_species                                                    │
     * │                                                                  │
     * │  INPUT:  logC = log10 of each component's concentration          │
     * │          (these are the unknowns the solver is trying to find)   │
     * │                                                                  │
     * │          initial_conc = original concentrations of ALL species    │
     * │          (used as fallback for non-equilibrium species like       │
     * │          biomass that this solver should not modify)              │
     * │                                                                  │
     * │  OUTPUT: concentrations of ALL species (in mol/L, not log)       │
     * │                                                                  │
     * │  WHAT IT DOES:                                                   │
     * │  For each species, uses the mass action law to compute its       │
     * │  concentration from the component concentrations.                │
     * │                                                                  │
     * │  EXAMPLE (carbonate system, components = {HCO3-, H+}):           │
     * │    If logC = [-3, -7] meaning [HCO3-] = 10^-3, [H+] = 10^-7:   │
     * │      log10([H2CO3]) = 6.35 + 1*(-3) + 1*(-7) = -3.65           │
     * │      [H2CO3] = 10^-3.65 = 2.24e-4 mol/L                        │
     * └──────────────────────────────────────────────────────────────────┘
     */
    std::vector<T> calc_species(const std::vector<T>& logC,
                                const std::vector<T>& initial_conc = std::vector<T>()) const {
        size_t ns = species_names.size();   /* Total number of species */
        size_t nc = component_names.size(); /* Number of components    */
        std::vector<T> conc(ns);            /* Output: one concentration per species */

        for (size_t i = 0; i < ns; ++i) {  /* Loop over each species */

            /* ── CASE 1: Non-equilibrium species (e.g., biomass, tracers) ──
             * These have all-zero rows in the stoichiometry matrix.
             * We don't change them — just keep their original concentration. */
            if (!isEquilibriumSpecies(i)) {
                conc[i] = (i < initial_conc.size() && initial_conc[i] > MIN_CONC)
                          ? std::min(initial_conc[i], MAX_CONC) : MIN_CONC;
                continue;   /* Skip to the next species */
            }

            /* ── CASE 2: Equilibrium species — apply mass action law ──
             *
             * We compute:  log10([species_i]) = logK[i] + Σⱼ S[i][j] * logC[j]
             *
             * Step by step:
             *   1. Start with logK (the equilibrium constant in log form)
             *   2. Add each component's contribution: S[i][j] * logC[j]
             *      This is like multiplying by [component_j]^S[i][j] in real space
             *   3. Convert back from log: concentration = 10^(result)
             */
            T logConc = (i < logK_values.size()) ? logK_values[i] : 0.0;  /* Start with logK */

            for (size_t j = 0; j < nc; ++j) {   /* Add each component's contribution */
                if (j < stoich_matrix[i].size()) {
                    logConc += stoich_matrix[i][j] * logC[j];
                }
            }

            /* Clamp the log-concentration to prevent:
             *   - overflow: 10^(300) = infinity (bad!)
             *   - underflow: 10^(-300) ≈ 0 (lost precision) */
            logConc = std::max(std::min(logConc, MAX_LOG_C), MIN_LOG_C);

            /* Convert from log-space back to real concentration */
            conc[i] = std::pow(10.0, logConc);  /* 10 raised to the power logConc */

            /* Safety net: if floating-point math produced Not-a-Number or infinity,
             * fall back to the initial concentration */
            if (std::isnan(conc[i]) || std::isinf(conc[i])) {
                conc[i] = (i < initial_conc.size() && initial_conc[i] > MIN_CONC)
                          ? initial_conc[i] : MIN_CONC;
            }

            /* Final clamp to physical bounds */
            conc[i] = std::max(std::min(conc[i], MAX_CONC), MIN_CONC);
        }
        return conc;
    }


    /* ╔════════════════════════════════════════════════════════════════════╗
     * ║                                                                    ║
     * ║  SECTION 2: COMPONENT TOTALS (conserved quantities)                ║
     * ║  ══════════════════════════════════════════════                    ║
     * ║                                                                    ║
     * ║  The "total" of a component is the sum of that component's         ║
     * ║  contribution across ALL species.                                  ║
     * ║                                                                    ║
     * ║  Example (carbonate system, component = carbon ≈ HCO3-):           ║
     * ║    T_carbon = 1*[H2CO3] + 1*[HCO3-] + 1*[CO3--]                  ║
     * ║             = all the carbon, regardless of chemical form           ║
     * ║                                                                    ║
     * ║  These totals are CONSERVED during equilibrium. Transport and      ║
     * ║  kinetics can change them, but equilibrium only shuffles mass      ║
     * ║  between species. The solver must find species concentrations      ║
     * ║  that satisfy mass action AND preserve these totals.               ║
     * ║                                                                    ║
     * ║  Formula: T_total[j] = Σᵢ S[i][j] × [species_i]                  ║
     * ║                                                                    ║
     * ╚════════════════════════════════════════════════════════════════════╝
     */
    std::vector<T> calc_component_totals(const std::vector<T>& species_conc) const {
        size_t nc = component_names.size();   /* Number of components */
        size_t ns = species_names.size();     /* Number of species    */
        std::vector<T> T_total(nc, 0.0);      /* Output: one total per component, start at 0 */

        /* For each component j, sum up contributions from all equilibrium species */
        for (size_t j = 0; j < nc; ++j) {
            for (size_t i = 0; i < ns && i < species_conc.size(); ++i) {
                /* Only include species that participate in equilibrium */
                if (isEquilibriumSpecies(i) && j < stoich_matrix[i].size()) {
                    T conc_safe = std::max(std::min(species_conc[i], MAX_CONC), MIN_CONC);
                    T_total[j] += stoich_matrix[i][j] * conc_safe;
                    /* ↑ Add S[i][j] × concentration of species i */
                }
            }
            T_total[j] = std::max(T_total[j], MIN_CONC);   /* Ensure positive */
        }
        return T_total;
    }


    /* ╔════════════════════════════════════════════════════════════════════╗
     * ║                                                                    ║
     * ║  SECTION 3: THE PCF METHOD                                         ║
     * ║  (Positive Continuous Fraction, Carrayrou et al. 2002)             ║
     * ║  ═════════════════════════════════════════════════                 ║
     * ║                                                                    ║
     * ║  HOW PCF WORKS (simple explanation):                               ║
     * ║                                                                    ║
     * ║  The mass balance equation for component j is:                     ║
     * ║    Σ S[i][j] * c[i] = T_total[j]                                  ║
     * ║                                                                    ║
     * ║  Instead of solving this directly (like Newton-Raphson would),     ║
     * ║  PCF rewrites it as a RATIO of two positive quantities:            ║
     * ║                                                                    ║
     * ║    S_reactive = sum of species that "use" this component           ║
     * ║    S_product  = sum of species that "provide" this component       ║
     * ║                                                                    ║
     * ║  At equilibrium: S_reactive = S_product (the ratio is 1).          ║
     * ║  If they're not equal, the LOG of their ratio tells us             ║
     * ║  which direction to adjust the component concentration.            ║
     * ║                                                                    ║
     * ║                                                 ║
     * ║  - Both sums are guaranteed POSITIVE (no sign issues)              ║
     * ║  - The log-ratio gives a smooth, bounded correction                ║
     * ║  - No derivatives or Jacobian needed (unlike Newton-Raphson)       ║
     * ║  - Always well-defined (no division by zero)                       ║
     * ║                                                                    ║
     * ╚════════════════════════════════════════════════════════════════════╝
     */

    /*
     * ┌──────────────────────────────────────────────────────────────────┐
     * │  calc_reactive_product_sums                                      │
     * │                                                                  │
     * │  For each component j, splits the stoichiometric contributions   │
     * │  into two groups:                                                │
     * │    - "positive contributions" (species that contain/produce j)   │
     * │    - "negative contributions" (species that consume j)           │
     * │                                                                  │
     * │  Then arranges them so both S_reactive and S_product are > 0.    │
     * │                                                                  │
     * │  Returns a pair: (S_reactive, S_product)                         │
     * └──────────────────────────────────────────────────────────────────┘
     */
    std::pair<std::vector<T>, std::vector<T>> calc_reactive_product_sums(
            const std::vector<T>& species_conc, const std::vector<T>& T_total) const {

        size_t nc = component_names.size();
        size_t ns = species_names.size();
        std::vector<T> S_reactive(nc, 0.0), S_product(nc, 0.0);

        for (size_t j = 0; j < nc; ++j) {   /* For each component */

            T sum_positive = 0.0;   /* Species with positive stoichiometry (producers)  */
            T sum_negative = 0.0;   /* Species with negative stoichiometry (consumers)  */

            /* Loop over all species and sort their contributions */
            for (size_t i = 0; i < ns && i < species_conc.size(); ++i) {
                /* Skip species that don't participate in equilibrium */
                if (!isEquilibriumSpecies(i) || j >= stoich_matrix[i].size()) continue;

                T mu_ij = stoich_matrix[i][j];   /* Stoichiometric coefficient */
                T c_i = std::max(species_conc[i], MIN_CONC);   /* Species concentration */

                /* Positive stoichiometry: this species CONTAINS component j
                 * Example: HCO3- has S[HCO3-][carbon] = +1 (it contains carbon) */
                if (mu_ij > 1e-15) {
                    sum_positive += mu_ij * c_i;
                }
                /* Negative stoichiometry: this species CONSUMES component j
                 * (rare, but can happen with certain reaction formulations) */
                else if (mu_ij < -1e-15) {
                    sum_negative += std::abs(mu_ij) * c_i;
                }
            }

            /* Arrange into reactive and product sums.
             * The arrangement depends on whether T_total is positive or negative.
             * The goal: both S_reactive and S_product must be > 0. */
            if (T_total[j] >= 0.0) {
                S_reactive[j] = sum_positive;                /* What we currently have */
                S_product[j] = T_total[j] + sum_negative;    /* What we need to match  */
            } else {
                S_reactive[j] = std::abs(T_total[j]) + sum_positive;
                S_product[j] = sum_negative;
            }

            /* Floor at MIN_CONC to prevent log10(0) = -infinity later */
            S_reactive[j] = std::max(S_reactive[j], MIN_CONC);
            S_product[j] = std::max(S_product[j], MIN_CONC);
        }

        /* Return both vectors as a pair */
        return std::make_pair(S_reactive, S_product);
    }

    /*
     * ┌──────────────────────────────────────────────────────────────────┐
     * │  pcf_residual — THE HEART OF THE PCF METHOD                      │
     * │                                                                  │
     * │  INPUT:  logC = current guess for log10(component concentrations) │
     * │          T_total = conserved component totals                     │
     * │          initial_conc = original species concentrations           │
     * │                                                                  │
     * │  OUTPUT: f = correction vector. If we add f to logC, we get a    │
     * │          better guess. If f ≈ 0, we've found equilibrium.        │
     * │                                                                  │
     * │  THE FORMULA for each component j:                               │
     * │                                                                  │
     * │    f[j] = (1/μ_min) × log10( S_product[j] / S_reactive[j] )     │
     * │                                                                  │
     * │  INTUITION:                                                      │
     * │  - If S_product > S_reactive → f > 0 → increase concentration   │
     * │    "We need more of this component, push it up"                  │
     * │  - If S_product < S_reactive → f < 0 → decrease concentration   │
     * │    "We have too much, push it down"                              │
     * │  - If S_product = S_reactive → f = 0 → EQUILIBRIUM!             │
     * │    "Everything balances perfectly"                                │
     * │                                                                  │
     * │  μ_min is the smallest positive stoichiometry coefficient for    │
     * │  component j. It's a scaling factor to keep units consistent.    │
     * └──────────────────────────────────────────────────────────────────┘
     */
    std::vector<T> pcf_residual(const std::vector<T>& logC, const std::vector<T>& T_total,
                                const std::vector<T>& initial_conc) const {
        size_t nc = component_names.size();

        /* Step 1: From the current guess logC, compute all species concentrations */
        std::vector<T> species_conc = calc_species(logC, initial_conc);

        /* Step 2: Split into reactive and product sums */
        auto sums = calc_reactive_product_sums(species_conc, T_total);
        const std::vector<T>& S_R = sums.first;    /* "What we have" (reactive)  */
        const std::vector<T>& S_P = sums.second;   /* "What we need" (product)   */

        /* Step 3: Compute the correction for each component */
        std::vector<T> f(nc);
        for (size_t j = 0; j < nc; ++j) {

            /* Find the smallest positive stoichiometry coefficient for this component.
             * This is used as a scaling factor (μ_min) in the PCF formula. */
            T mu_i0_j = 1.0;
            for (size_t i = 0; i < stoich_matrix.size(); ++i) {
                if (!isEquilibriumSpecies(i)) continue;
                if (j < stoich_matrix[i].size()) {
                    T mu = stoich_matrix[i][j];
                    if (mu > 1e-10 && mu < mu_i0_j) mu_i0_j = mu;
                }
            }

            /* THE PCF CORRECTION FORMULA:
             *
             *   f[j] = (1 / μ_min) × ( log10(S_product) - log10(S_reactive) )
             *        = (1 / μ_min) × log10( S_product / S_reactive )
             *
             * When S_product == S_reactive → log10(1) = 0 → no correction → equilibrium!
             * When S_product > S_reactive  → positive → increase this component
             * When S_product < S_reactive  → negative → decrease this component
             */
            f[j] = (1.0 / mu_i0_j) * (std::log10(S_P[j]) - std::log10(S_R[j]));

            /* Clamp correction to ±10 orders of magnitude per step.
             * This prevents the solver from making wild jumps. */
            f[j] = std::max(std::min(f[j], 10.0), -10.0);

            /* Replace NaN or infinity with zero (no correction = safe fallback) */
            if (std::isnan(f[j]) || std::isinf(f[j])) f[j] = 0.0;
            /* [EQ-FIX 2026-07, threshold corrected 2026-07-24] Zero-total component guard.
             * A component absent from this water (Fe2, U6, HS, Ac in the initial water) still has
             * trace species, so S_R > 0 while S_P floors to MIN_CONC -> log10(S_P/S_R) is a large
             * spurious residual that never reaches 0 (its free ion is pinned at the MIN_LOG_C floor).
             * That kept the residual norm above tolerance forever even when every REAL component was
             * solved (pH exactly 7), so the solver reported non-convergence and produced garbage.
             * IMPORTANT: calc_component_totals SUMS a MIN_CONC floor over EVERY species that contains
             * the component, so an absent component total is (n_species_with_j)*MIN_CONC, e.g. Fe2 is
             * in 13 species -> 1.3e-29, U6 in ~20 -> ~2e-29. Those are NOT <= MIN_CONC (1e-30), so the
             * old "<= MIN_CONC" test never fired and the residual stalled near 3.0. Use an absolute
             * ABSENT threshold well above (n_species)*MIN_CONC (~1e-28) yet far below any real total
             * (smallest real component here is HPO4 = 1e-5): 1e-18 clears the artifact by 10 orders
             * and sits 13 orders below any physically meaningful total. Once transport/kinetics push a
             * total above 1e-18 the component is solved normally again. Verified: converges in 30
             * iters to pH 7.03 with all real components balanced to ~1e-13 (was: stall at residual 3.09). */
            if (T_total[j] <= 1e-18) f[j] = 0.0;
        }
        return f;
    }


    /* ╔════════════════════════════════════════════════════════════════════╗
     * ║                                                                    ║
     * ║  SECTION 4: HELPER MATH FUNCTIONS                                  ║
     * ║  ════════════════════════════════                                  ║
     * ║                                                                    ║
     * ║  Basic vector operations used throughout the solver.               ║
     * ║  "Vector" here means a list of numbers, not a physical direction.  ║
     * ║                                                                    ║
     * ╚════════════════════════════════════════════════════════════════════╝
     */

    /* norm: Compute the length (Euclidean norm) of a vector.
     * Formula: ||v|| = sqrt(v[0]² + v[1]² + ... + v[n-1]²)
     * Used to check convergence: if norm(residual) < tolerance → done!
     */
    T norm(const std::vector<T>& v) const {
        T s = 0;
        for (T x : v) s += x * x;   /* Sum of squares */
        return std::sqrt(s);          /* Square root of sum */
    }

    /* vec_subtract: Compute a - b element-by-element.
     * Example: {3, 5, 1} - {1, 2, 0} = {2, 3, 1} */
    std::vector<T> vec_subtract(const std::vector<T>& a, const std::vector<T>& b) const {
        std::vector<T> r(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
            r[i] = a[i] - (i < b.size() ? b[i] : 0.0);
        }
        return r;
    }

    /* vec_add: Compute a + b element-by-element.
     * Example: {3, 5, 1} + {1, 2, 0} = {4, 7, 1} */
    std::vector<T> vec_add(const std::vector<T>& a, const std::vector<T>& b) const {
        std::vector<T> r(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
            r[i] = a[i] + (i < b.size() ? b[i] : 0.0);
        }
        return r;
    }

    /* dot_product: Multiply corresponding elements and sum them.
     * Example: {3, 5} · {2, 4} = 3×2 + 5×4 = 26
     * Used in QR decomposition to measure projections between vectors. */
    T dot_product(const std::vector<T>& a, const std::vector<T>& b) const {
        T sum = 0;
        for (size_t i = 0; i < a.size() && i < b.size(); ++i) {
            sum += a[i] * b[i];
        }
        return sum;
    }


    /* ╔════════════════════════════════════════════════════════════════════╗
     * ║                                                                    ║
     * ║  SECTION 5: QR DECOMPOSITION                                       ║
     * ║  ═══════════════════════════                                       ║
     * ║                                                                    ║
     * ║  WHAT IS QR DECOMPOSITION?                                         ║
     * ║  It's a way to break a matrix into two simpler pieces:             ║
     * ║    Matrix = Q × R                                                  ║
     * ║  where:                                                            ║
     * ║    Q = a matrix with perpendicular unit-length columns             ║
     * ║    R = an upper triangular matrix (zeros below the diagonal)       ║
     * ║                                                                    ║
     * ║  WHY DO WE NEED IT?                                                ║
     * ║  Anderson Acceleration needs to solve a "least squares" problem:   ║
     * ║  "Find the best combination of past corrections."                  ║
     * ║  QR decomposition turns this into solving a simple triangular      ║
     * ║  equation, which is fast and numerically stable.                   ║
     * ║                                                                    ║
     * ║  HOW IT WORKS (Gram-Schmidt process):                              ║
     * ║  Given vectors v1, v2, v3, ... we create perpendicular vectors     ║
     * ║  q1, q2, q3, ... by:                                               ║
     * ║    1. q1 = v1, then normalize it (make length = 1)                 ║
     * ║    2. q2 = v2 minus its projection onto q1, then normalize         ║
     * ║    3. q3 = v3 minus its projections onto q1 and q2, then normalize ║
     * ║  The R matrix records how much of each original vector projected   ║
     * ║  onto each q.                                                      ║
     * ║                                                                    ║
     * ║  CONDITION NUMBER:                                                  ║
     * ║  The ratio of the largest to smallest diagonal element of R.       ║
     * ║  If this ratio is huge (>10^10), the problem is "ill-conditioned"  ║
     * ║  and Anderson should use fewer past iterations to stay safe.       ║
     * ║  Reference: Walker & Ni (2011).                                    ║
     * ║                                                                    ║
     * ╚════════════════════════════════════════════════════════════════════╝
     */
    void qr_decomposition(const std::vector<std::vector<T>>& columns,
                          std::vector<std::vector<T>>& Q,
                          std::vector<std::vector<T>>& R, T& condNum) const {
        size_t m = columns.size();          /* Number of columns (= Anderson depth) */
        if (m == 0) { condNum = 1.0; return; }
        size_t n = columns[0].size();       /* Length of each column (= number of components) */

        Q.resize(m);                                    /* Q will have m columns */
        R.resize(m, std::vector<T>(m, 0.0));            /* R is m × m, starts as zeros */
        T R_max = 0.0, R_min = std::numeric_limits<T>::max();   /* For condition number */

        /* Process each column one at a time */
        for (size_t j = 0; j < m; ++j) {
            Q[j] = columns[j];   /* Start with the original column vector */

            /* GRAM-SCHMIDT: Remove projections onto all previous Q vectors.
             * This makes Q[j] perpendicular to Q[0], Q[1], ..., Q[j-1].
             *
             * Think of it like this: if you have two vectors that are almost
             * parallel, subtract the parallel part to keep only the perpendicular part.
             */
            for (size_t i = 0; i < j; ++i) {
                /* R[i][j] = how much of the original column points in the Q[i] direction */
                R[i][j] = dot_product(Q[i], columns[j]);

                /* Subtract that projection from Q[j] */
                for (size_t k = 0; k < n; ++k) {
                    Q[j][k] -= R[i][j] * Q[i][k];
                }
            }

            /* R[j][j] = length of what's left after removing all projections */
            R[j][j] = norm(Q[j]);

            /* Normalize Q[j] to unit length (divide by its length) */
            if (R[j][j] > 1e-15) {   /* Only if the length is nonzero */
                for (size_t k = 0; k < n; ++k) {
                    Q[j][k] /= R[j][j];
                }
            }

            /* Track the largest and smallest diagonal elements of R.
             * Their ratio is the condition number. */
            T absR = std::abs(R[j][j]);
            if (absR > R_max) R_max = absR;
            if (absR > 1e-30 && absR < R_min) R_min = absR;
        }

        /* Condition number = max diagonal / min diagonal.
         * Large value → ill-conditioned → Anderson should reduce depth. */
        condNum = (R_min > 1e-30) ? R_max / R_min : std::numeric_limits<T>::max();
    }

    /*
     * ┌──────────────────────────────────────────────────────────────────┐
     * │  solve_upper_triangular                                          │
     * │                                                                  │
     * │  Solves the equation R × x = b for x, where R is an upper       │
     * │  triangular matrix (all zeros below the diagonal).               │
     * │                                                                  │
     * │  HOW IT WORKS (back substitution):                               │
     * │  Start from the last row (which has only one unknown) and work   │
     * │  upward, plugging in known values as you go.                     │
     * │                                                                  │
     * │  Example for a 2×2 system:                                       │
     * │    R = [[3, 2],   b = [11, 4]                                    │
     * │         [0, 4]]                                                  │
     * │    From row 2: 4*x[1] = 4  →  x[1] = 1                         │
     * │    From row 1: 3*x[0] + 2*1 = 11  →  x[0] = 3                  │
     * │                                                                  │
     * │  This is much faster than general equation solving because       │
     * │  we already know R is triangular (thanks to QR decomposition).   │
     * └──────────────────────────────────────────────────────────────────┘
     */
    std::vector<T> solve_upper_triangular(const std::vector<std::vector<T>>& R,
                                          const std::vector<T>& b) const {
        size_t m = R.size();
        std::vector<T> x(m, 0.0);   /* Solution vector, initialized to zeros */

        /* Work backwards: start from the last row, go to the first */
        for (int i = m - 1; i >= 0; --i) {
            x[i] = b[i];   /* Start with the right-hand side */

            /* Subtract all already-known terms */
            for (size_t j = i + 1; j < m; ++j) {
                x[i] -= R[i][j] * x[j];
            }

            /* Divide by the diagonal element to get x[i] */
            if (std::abs(R[i][i]) > 1e-30) {
                x[i] /= R[i][i];
            }
        }
        return x;
    }


    /* ╔════════════════════════════════════════════════════════════════════╗
     * ║                                                                    ║
     * ║  SECTION 6: THE MAIN SOLVER — Anderson Acceleration + PCF          ║
     * ║  ════════════════════════════════════════════════════              ║
     * ║                                                                    ║
     * ║  This is the core function that combines everything above.         ║
     * ║                                                                    ║
     * ║  PLAIN-ENGLISH ALGORITHM:                                          ║
     * ║                                                                    ║
     * ║  1. Start with a guess for component concentrations                ║
     * ║  2. Use the PCF formula to compute a "correction" (how far         ║
     * ║     off we are from equilibrium)                                   ║
     * ║  3. Check: is the correction tiny? If yes → DONE! (converged)     ║
     * ║  4. If not, use Anderson Acceleration to compute a SMARTER         ║
     * ║     correction by looking at the pattern in recent corrections     ║
     * ║  5. Apply the correction, go back to step 2                        ║
     * ║                                                                    ║
     * ║  WHAT ANDERSON ACCELERATION DOES (analogy):                        ║
     * ║  Imagine you're trying to find the right temperature for your      ║
     * ║  shower. Without AA, you'd adjust the knob a little each time      ║
     * ║  and wait (simple iteration — slow). With AA, you REMEMBER the     ║
     * ║  last few adjustments and their results, notice the pattern, and   ║
     * ║  make a BIGGER, SMARTER adjustment that gets you to the right      ║
     * ║  temperature faster. That's exactly what AA does mathematically.   ║
     * ║                                                                    ║
     * ║  SPECIFICALLY, AA solves a "least-squares" problem each step:      ║
     * ║  "Given the last m residuals and iterates, find the weighted       ║
     * ║  combination that minimizes the next residual."                    ║
     * ║  It does this via QR decomposition (Section 5 above).              ║
     * ║                                                                    ║
     * ║  SAFETY FEATURES:                                                  ║
     * ║  - If the history becomes ill-conditioned (condition number too    ║
     * ║    high), AA automatically reduces its depth (forgets old data)    ║
     * ║  - All values are clamped to physical bounds                       ║
     * ║  - NaN/Inf values are replaced with the previous iterate           ║
     * ║                                                                    ║
     * ╚════════════════════════════════════════════════════════════════════╝
     */
    std::vector<T> solve_equilibrium_anderson(const std::vector<T>& initial_species_conc,
                                              const std::vector<T>& T_total) {
        /* Reset tracking variables for this solve */
        lastConverged = false;
        lastIterations = 0;
        lastResidual = 0.0;
        totalSolves++;   /* Increment the "how many times called" counter */

        size_t nc = component_names.size();   /* Number of components (= number of unknowns) */

        /* If there are no components, there's nothing to solve */
        if (nc == 0) {
            lastConverged = true;
            totalConverged++;
            return initial_species_conc;
        }

        /*
         * ──── STEP 1: Build the initial guess ────
         *
         * The unknowns are ω[j] = log10([component_j]).
         * We initialize them from the current species concentrations.
         * For each component, find its concentration in the species list
         * and take the log10. If not found, default to 10^-7 (pH 7 guess).
         */
        std::vector<T> omega(nc);   /* ω = the current guess (log-concentrations) */

        for (size_t j = 0; j < nc; ++j) {
            /* Search for this component name in the species list */
            auto it = std::find(species_names.begin(), species_names.end(), component_names[j]);
            T comp_conc = 1e-7;   /* Default guess if not found */

            if (it != species_names.end()) {   /* Found it! */
                size_t idx = std::distance(species_names.begin(), it);   /* Its index */
                if (idx < initial_species_conc.size()) {
                    comp_conc = initial_species_conc[idx];   /* Use actual concentration */
                }
            }

            /* Clamp to physical bounds, then take log10 */
            comp_conc = std::max(std::min(comp_conc, MAX_CONC), MIN_CONC);
            omega[j] = std::log10(comp_conc);
        }

        /*
         * ──── STEP 2: First PCF step (plain fixed-point, no AA yet) ────
         *
         * Before Anderson Acceleration can work, it needs at least one
         * previous iterate in its history. So we do one simple PCF step:
         *   ω₁ = ω₀ + f(ω₀)
         * where f is the PCF correction computed above.
         */
        std::vector<T> f0 = pcf_residual(omega, T_total, initial_species_conc);
        std::vector<T> omega_new = vec_add(omega, f0);   /* New guess = old + correction */

        /* Clamp to bounds */
        for (size_t j = 0; j < nc; ++j) {
            omega_new[j] = std::max(std::min(omega_new[j], MAX_LOG_C), MIN_LOG_C);
        }

        /*
         * ──── STEP 3: Set up history storage for Anderson Acceleration ────
         *
         * AA remembers the last few iterates (ω values) and their residuals (f values).
         * It uses this history to extrapolate a better next guess.
         */
        std::vector<std::vector<T>> omega_history, f_history;
        omega_history.push_back(omega);   /* Remember ω₀ */
        f_history.push_back(f0);          /* Remember f(ω₀) */
        omega = omega_new;                /* Current iterate is now ω₁ */

        /*
         * ──── STEP 4: Main iteration loop ────
         *
         * Repeat until converged or max iterations reached.
         * Each iteration:
         *   a) Compute PCF residual at current guess
         *   b) Check convergence
         *   c) Use Anderson Acceleration (or simple PCF if not enough history)
         *   d) Update the guess
         */
        for (plint iter = 1; iter < maxIterations; ++iter) {
            lastIterations = iter;

            /* (a) Compute the PCF residual: "how far are we from equilibrium?" */
            std::vector<T> f_k = pcf_residual(omega, T_total, initial_species_conc);
            T f_norm = norm(f_k);   /* The length of the residual vector */
            lastResidual = f_norm;

            /* (b) CONVERGENCE CHECK:
             * If the residual is smaller than our tolerance (1e-8 by default),
             * the mass balance is satisfied to sufficient precision.
             * Convert log-concentrations back to real concentrations and return. */
            if (f_norm < tolerance) {
                lastConverged = true;
                totalConverged++;
                return calc_species(omega, initial_species_conc);
            }

            /* Save current iterate and residual to history */
            omega_history.push_back(omega);
            f_history.push_back(f_k);

            /* Determine Anderson depth for this step: use min(requested depth, available history) */
            plint m_k = std::min(andersonDepth, (plint)(f_history.size() - 1));

            if (m_k < 1) {
                /*
                 * NOT ENOUGH HISTORY: Fall back to simple PCF fixed-point.
                 * ω_new = ω + f(ω)
                 * This is the plain iteration without acceleration.
                 */
                omega_new = vec_add(omega, f_k);
            } else {
                /*
                 * ══════════════════════════════════════════════════════════
                 *  ANDERSON ACCELERATION
                 * ══════════════════════════════════════════════════════════
                 *
                 *  AA uses the DIFFERENCES between consecutive iterates and
                 *  residuals to build a better update.
                 */

                /* Build difference vectors from history:
                 *   ΔF[i] = f_{k-m+i+1} - f_{k-m+i}   (how residuals changed)
                 *   ΔX[i] = ω_{k-m+i+1} - ω_{k-m+i}   (how iterates changed)
                 *
                 * These tell AA "in what direction did the solver move,
                 * and how did that affect the residual?"
                 */
                std::vector<std::vector<T>> Delta_F, Delta_X;
                size_t hist_size = f_history.size();

                for (plint i = 0; i < m_k; ++i) {
                    size_t idx = hist_size - m_k - 1 + i;
                    Delta_F.push_back(vec_subtract(f_history[idx + 1], f_history[idx]));
                    Delta_X.push_back(vec_subtract(omega_history[idx + 1], omega_history[idx]));
                }

                /* QR-decompose the ΔF matrix (see Section 5 for what QR does).
                 * Also compute the condition number to check numerical health. */
                std::vector<std::vector<T>> Q, R;
                T condNum;
                qr_decomposition(Delta_F, Q, R, condNum);

                /* SELF-HEALING: If condition number is too large, the history
                 * vectors are nearly parallel (numerically dangerous).
                 * Solution: drop the oldest history entries until it's safe.
                 * This is the "safety valve" that makes AA robust. */
                while (condNum > conditionTol && Delta_F.size() > 1) {
                    Delta_F.erase(Delta_F.begin());   /* Drop oldest */
                    Delta_X.erase(Delta_X.begin());
                    qr_decomposition(Delta_F, Q, R, condNum);   /* Recompute */
                }

                /* SOLVE THE LEAST-SQUARES PROBLEM:
                 * "Find weights γ that minimize ||f_k - Σ γᵢ × ΔFᵢ||²"
                 *
                 * With QR, this becomes: R × γ = Q^T × f_k
                 * (a triangular system, solved by back substitution) */

                /* First, compute Q^T × f_k (project f_k onto each Q column) */
                std::vector<T> Qt_fk(Delta_F.size());
                for (size_t i = 0; i < Delta_F.size(); ++i) {
                    Qt_fk[i] = dot_product(Q[i], f_k);
                }

                /* Then solve the triangular system for γ */
                std::vector<T> gamma = solve_upper_triangular(R, Qt_fk);

                /* COMPUTE THE ANDERSON UPDATE:
                 *
                 *   ω_new = ω_k - Σ γᵢ×ΔXᵢ + β × (f_k - Σ γᵢ×ΔFᵢ)
                 *           ─────────────────   ──────────────────────────
                 *           "correct position    "apply optimized residual
                 *            using past steps"     correction"
                 *
                 * In plain English: "Based on the pattern in past iterations,
                 * extrapolate where the residual would be zero."
                 *
                 * β (beta) is the mixing parameter. At β=1.0, we use the full
                 * Anderson step. Values < 1 would add damping.
                 */
                std::vector<T> DX_gamma(nc, 0.0);   /* Σ γᵢ × ΔXᵢ */
                std::vector<T> DF_gamma(nc, 0.0);   /* Σ γᵢ × ΔFᵢ */

                for (size_t i = 0; i < Delta_X.size() && i < gamma.size(); ++i) {
                    for (size_t j = 0; j < nc; ++j) {
                        DX_gamma[j] += Delta_X[i][j] * gamma[i];
                        DF_gamma[j] += Delta_F[i][j] * gamma[i];
                    }
                }

                omega_new.resize(nc);
                for (size_t j = 0; j < nc; ++j) {
                    omega_new[j] = omega[j] - DX_gamma[j] + beta * (f_k[j] - DF_gamma[j]);
                }
            }

            /* Clamp the new guess to physical bounds and replace any NaN/Inf */
            for (size_t j = 0; j < nc; ++j) {
                omega_new[j] = std::max(std::min(omega_new[j], MAX_LOG_C), MIN_LOG_C);
                if (std::isnan(omega_new[j]) || std::isinf(omega_new[j])) {
                    omega_new[j] = omega[j];   /* Fall back to previous value */
                }
            }
            omega = omega_new;   /* Accept the new guess */

            /* Trim history: only keep the most recent (depth + 1) entries.
             * Older history is no longer useful and wastes memory. */
            while (omega_history.size() > (size_t)(andersonDepth + 1)) {
                omega_history.erase(omega_history.begin());   /* Remove oldest */
                f_history.erase(f_history.begin());
            }
        }

        /* If we get here, we used up all iterations without converging.
         * Return the best guess we have. The caller can check didConverge(). */
        totalDiverged++;
        return calc_species(omega, initial_species_conc);
    }


    /* ╔════════════════════════════════════════════════════════════════════╗
     * ║                                                                    ║
     * ║  SECTION 7: MAIN ENTRY POINT                                       ║
     * ║  ═══════════════════════════                                       ║
     * ║                                                                    ║
     * ║  This is THE function that gets called from the simulation.         ║
     * ║  It takes current species concentrations, computes the conserved   ║
     * ║  totals, runs the Anderson+PCF solver, and returns equilibrated    ║
     * ║  concentrations.                                                   ║
     * ║                                                                    ║
     * ╚════════════════════════════════════════════════════════════════════╝
     */
    std::vector<T> calculate_species_concentrations(const std::vector<T>& initial_conc) {
        size_t nc = component_names.size();
        size_t ns = species_names.size();

        /* If there are no species or components, nothing to solve */
        if (nc == 0 || ns == 0) return initial_conc;

        /* Step 1: Compute the conserved component totals.
         * These are the "target" values that the solver must preserve. */
        std::vector<T> T_total = calc_component_totals(initial_conc);

        /* Step 2: Run the Anderson+PCF solver.
         * This finds new species concentrations that satisfy both
         * the mass action law AND the conserved totals. */
        std::vector<T> result = solve_equilibrium_anderson(initial_conc, T_total);

        /* Step 3: Validate the result — catch any numerical problems */
        for (size_t i = 0; i < result.size(); ++i) {
            if (std::isnan(result[i]) || std::isinf(result[i])) {
                /* If any species has a bad value, use its original concentration */
                result[i] = (i < initial_conc.size() && initial_conc[i] > MIN_CONC)
                            ? std::min(initial_conc[i], MAX_CONC) : MIN_CONC;
            }
            /* Final clamp to physical bounds */
            result[i] = std::max(std::min(result[i], MAX_CONC), MIN_CONC);
        }
        return result;
    }

    /*
     * ┌──────────────────────────────────────────────────────────────────┐
     * │  printStatistics: Print a summary of how the solver performed.   │
     * │  Shows total solves, convergence count, and failure count.       │
     * │  Called at the end of the simulation for diagnostic purposes.    │
     * └──────────────────────────────────────────────────────────────────┘
     */
    void printStatistics() const {
        pcout << "╔══════════════════════════════════════════════════════════════════╗\n";
        pcout << "║  ANDERSON ACCELERATION + PCF SOLVER STATISTICS                   ║\n";
        pcout << "╠══════════════════════════════════════════════════════════════════╣\n";
        pcout << "║  Total solves:    " << std::setw(10) << totalSolves
              << "                                ║\n";
        pcout << "║  Converged:       " << std::setw(10) << totalConverged << " ("
              << std::fixed << std::setprecision(1)
              << (totalSolves > 0 ? 100.0 * totalConverged / totalSolves : 0.0)
              << "%)                       ║\n";
        pcout << "║  Did not converge:" << std::setw(10) << totalDiverged
              << "                                ║\n";
        pcout << "╚══════════════════════════════════════════════════════════════════╝\n";
    }

/*
 * ┌──────────────────────────────────────────────────────────────────┐
 * │  PRIVATE DATA MEMBERS                                            │
 * │                                                                  │
 * │  "private:" means these variables can only be accessed from      │
 * │  inside this class, not from outside code.                       │
 * └──────────────────────────────────────────────────────────────────┘
 */
private:
    std::vector<std::string> species_names;    /* List of all species names           */
    std::vector<std::string> component_names;  /* List of component (master) names    */
    std::vector<T> logK_values;                /* log10(K) for each species           */
    std::vector<std::vector<T>> stoich_matrix; /* Stoichiometry matrix S[species][component] */

    plint maxIterations;   /* Maximum iterations before giving up (default 200) */
    T tolerance;           /* Convergence criterion (default 1e-8)              */

    /* "mutable" means these CAN be changed even in "const" functions.
     * This is because the solver tracks statistics about its performance
     * even when called from a const context. */
    mutable bool lastConverged;     /* Did the most recent solve converge?       */
    mutable plint lastIterations;   /* How many iterations did the last solve take? */
    mutable T lastResidual;         /* What was the final residual of the last solve? */

    bool verbose;   /* If true, print detailed output during iteration */

    mutable plint totalSolves;     /* Total number of solves attempted     */
    mutable plint totalConverged;  /* Total number of successful solves    */
    mutable plint totalDiverged;   /* Total number of failed solves        */

    plint andersonDepth;   /* How many past iterations AA remembers (default 4) */
    T conditionTol;        /* Max allowed condition number (default 1e10)       */
    T beta;                /* Mixing/damping parameter (default 1.0 = no damping) */
};


/* ╔════════════════════════════════════════════════════════════════════════════╗
 * ║                                                                          ║
 * ║  SECTION 8: LATTICE BOLTZMANN DATA PROCESSORS                           ║
 * ║  ═════════════════════════════════════════════                           ║
 * ║                                                                          ║
 * ║  These classes are "data processors" in the Palabos framework.           ║
 * ║  They are the bridge between the equilibrium solver (above) and the      ║
 * ║  LBM simulation grid. A data processor is a piece of code that Palabos   ║
 * ║  applies to every voxel (3D grid cell) in the domain.                    ║
 * ║                                                                          ║
 * ║  HOW PALABOS LATTICES WORK (simplified):                                 ║
 * ║  - Each chemical species has its own 3D grid (a "lattice").              ║
 * ║  - At each grid cell, the lattice stores 7 numbers called "populations"  ║
 * ║    (for D3Q7: 1 rest + 6 directions: +x, -x, +y, -y, +z, -z).         ║
 * ║  - The SUM of these 7 populations = the concentration at that cell.      ║
 * ║  - To change the concentration by dC, we add dC*weight to each           ║
 * ║    population, where the weights are: 1/4 (rest), 1/8 (each direction). ║
 * ║  - These weights are the D3Q7 equilibrium weights, ensuring the          ║
 * ║    concentration change is isotropic (same in all directions).           ║
 * ║                                                                          ║
 * ╚════════════════════════════════════════════════════════════════════════════╝
 */


/*
 * ┌──────────────────────────────────────────────────────────────────────┐
 * │  run_equilibrium_biotic — RECOMMENDED equilibrium data processor     │
 * │                                                                      │
 * │  This processor loops over every pore voxel and:                     │
 * │    1. Reads substrate concentrations from the lattice populations    │
 * │    2. Calls the Anderson+PCF solver to equilibrate them              │
 * │    3. Writes the concentration CHANGES back into the lattice         │
 * │                                                                      │
 * │  LATTICE LAYOUT:                                                     │
 * │    lattices[0]            = first substrate (e.g., H2CO3)            │
 * │    lattices[1]            = second substrate (e.g., HCO3-)           │
 * │    ...                                                               │
 * │    lattices[subsNum-1]    = last substrate                           │
 * │    lattices[subsNum]      = mask lattice (tells us solid vs pore)    │
 * │                                                                      │
 * │  WHEN IS THIS CALLED?                                                │
 * │  After transport and kinetics, at each timestep:                     │
 * │    1. Transport (LBM advection-diffusion)                            │
 * │    2. Kinetics (run_kinetics — biological reactions)                 │
 * │    3. *** Equilibrium (THIS processor) ***                           │
 * └──────────────────────────────────────────────────────────────────────┘
 */
template<typename T, template<typename U> class Descriptor>
class run_equilibrium_biotic : public LatticeBoxProcessingFunctional3D<T, Descriptor>
{
public:
    /*
     * CONSTRUCTOR: Called once when the processor is created.
     *
     * Parameters:
     *   nx_       = total domain size in the x-direction (flow direction)
     *   subsNum_  = number of substrates (chemical species with lattices)
     *   eqChem_   = reference to the equilibrium solver object (configured earlier)
     *   solid_    = mask value that means "solid wall" (don't compute here)
     *   bb_       = mask value that means "bounce-back" (don't compute here)
     *
     * maskLloc = subsNum_ means the mask lattice is at index subsNum
     * (right after all the substrate lattices).
     */
    run_equilibrium_biotic(plint nx_, plint subsNum_, EquilibriumChemistry<T>& eqChem_,
                           plint solid_, plint bb_)
        : nx(nx_), subsNum(subsNum_), eqChem(eqChem_), solid(solid_), bb(bb_),
          maskLloc(subsNum_)
    {}

    /*
     * PROCESS: The main function that Palabos calls for each sub-block.
     *
     * "domain" = the rectangular region of grid cells this call handles.
     * "lattices" = the list of all lattices (substrates + mask).
     *
     * Palabos may split the domain across multiple CPUs. Each CPU calls
     * process() on its own sub-block. The "absoluteOffset" corrects for
     * the fact that each sub-block's local coordinates start at (0,0,0)
     * but the global coordinates may be different.
     */
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {

        /* Get the global offset so we know where this sub-block sits
         * in the overall simulation domain. */
        Dot3D absoluteOffset = lattices[0]->getLocation();

        /* ═══════════════════════════════════════════════════════════════
         * MAIN TRIPLE LOOP: Visit every voxel in this sub-block.
         * For each PORE voxel:
         *   1. Read concentrations
         *   2. Solve equilibrium
         *   3. Write changes back
         * ═══════════════════════════════════════════════════════════════ */
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;   /* Global x coordinate */

            /* Skip the inlet (x=0) and outlet (x=nx-1) boundary planes.
             * These are handled separately by boundary conditions. */
            if (absX > 0 && absX < nx - 1) {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                    for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {

                        /* READ THE MASK VALUE for this voxel.
                         * The mask tells us what type of node this is:
                         *   solid → skip (inside grain)
                         *   bb    → skip (bounce-back, at grain surface)
                         *   other → pore space, process it!
                         *
                         * "computeRelativeDisplacement" aligns the mask lattice
                         * with the substrate lattice (they may have different
                         * memory layouts in Palabos multi-block). */
                        Dot3D maskOffset = computeRelativeDisplacement(*lattices[0], *lattices[maskLloc]);
                        plint mask = util::roundToInt(lattices[maskLloc]->get(
                            iX + maskOffset.x, iY + maskOffset.y, iZ + maskOffset.z).computeDensity());

                        /* Only process PORE voxels */
                        if (mask != solid && mask != bb) {

                            /* Compute memory offsets for all lattices.
                             * Each lattice may be stored at a different memory location,
                             * so we need to know the offset for each one. */
                            std::vector<Dot3D> vec_offset;
                            for (plint iT = 0; iT <= maskLloc; ++iT) {
                                vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT]));
                            }

                            /* ── STEP 1: READ current substrate concentrations ──
                             * For each substrate, read its concentration at this voxel.
                             * "computeDensity()" sums the 7 D3Q7 populations to get
                             * the macroscopic concentration value. */
                            std::vector<T> conc(subsNum);
                            for (plint iS = 0; iS < subsNum; ++iS) {
                                plint iXs = iX + vec_offset[iS].x;
                                plint iYs = iY + vec_offset[iS].y;
                                plint iZs = iZ + vec_offset[iS].z;
                                T c0 = lattices[iS]->get(iXs, iYs, iZs).computeDensity();
                                conc[iS] = std::max(c0, EquilibriumChemistry<T>::MIN_CONC);
                            }

                            /* ── STEP 2: SOLVE equilibrium at this voxel ──
                             * Pass current concentrations to the solver.
                             * It returns new equilibrium concentrations. */
                            std::vector<T> eq_conc = eqChem.calculate_species_concentrations(conc);

                            /* ── STEP 3: WRITE concentration changes back ── */
                            for (plint iS = 0; iS < subsNum; ++iS) {
                                T dC = eq_conc[iS] - conc[iS];   /* Change = new - old */

                                /* ─── SAFEGUARD: Limit maximum change per timestep ───
                                 *
                                 * Why? If the equilibrium solver computes a huge change
                                 * (e.g., concentration should double), applying it all
                                 * at once could destabilize the LBM transport solver.
                                 * So we limit each step to at most:
                                 *   - 10% of the current concentration, OR
                                 *   - 1e-4 mol/L absolute
                                 * whichever is larger. The remaining change will happen
                                 * over subsequent timesteps.
                                 */
                                const T MAX_RELATIVE_CHANGE = 0.1;   /* 10% max relative change */
                                const T MAX_ABSOLUTE_CHANGE = 1e-9;  /* [BUG4] tiny floor (mol/L): micromolar pools (e.g. U ~2e-5) governed by the 10% relative cap, not overshot by a 1e-4 absolute step */

                                T max_allowed = std::max(MAX_ABSOLUTE_CHANGE,
                                                         MAX_RELATIVE_CHANGE * std::abs(conc[iS]));
                                if (std::abs(dC) > max_allowed) {
                                    dC = (dC > 0) ? max_allowed : -max_allowed;
                                }

                                /* Ensure the new concentration doesn't go negative */
                                T new_conc = conc[iS] + dC;
                                if (new_conc < 1e-20) {
                                    dC = 1e-20 - conc[iS];   /* Adjust dC to floor at 1e-20 */
                                }

                                /* Only update if the change is significant.
                                 * "thrd" = 1e-12, a tiny threshold to avoid
                                 * unnecessary floating-point operations. */
                                if (std::abs(dC) > thrd) {
                                    Array<T, 7> g;   /* Array to hold 7 D3Q7 populations */

                                    /* Get the lattice coordinates for this substrate */
                                    plint iXt = iX + vec_offset[iS].x;
                                    plint iYt = iY + vec_offset[iS].y;
                                    plint iZt = iZ + vec_offset[iS].z;

                                    /* Read the current 7 populations */
                                    lattices[iS]->get(iXt, iYt, iZt).getPopulations(g);

                                    /* ADD the change to each population, weighted by D3Q7 weights.
                                     *
                                     * D3Q7 weights:
                                     *   g[0] = rest population (stays at this node)   → weight 1/4
                                     *   g[1] = +x direction                           → weight 1/8
                                     *   g[2] = -x direction                           → weight 1/8
                                     *   g[3] = +y direction                           → weight 1/8
                                     *   g[4] = -y direction                           → weight 1/8
                                     *   g[5] = +z direction                           → weight 1/8
                                     *   g[6] = -z direction                           → weight 1/8
                                     *                           Total: 1/4 + 6×(1/8) = 1.0 ✓
                                     *
                                     * By distributing dC according to these weights, the
                                     * total (= concentration) changes by exactly dC. */
                                    g[0] += dC / 4.0;   /* Rest: stays here             */
                                    g[1] += dC / 8.0;   /* Streams to +x neighbor       */
                                    g[2] += dC / 8.0;   /* Streams to -x neighbor       */
                                    g[3] += dC / 8.0;   /* Streams to +y neighbor       */
                                    g[4] += dC / 8.0;   /* Streams to -y neighbor       */
                                    g[5] += dC / 8.0;   /* Streams to +z neighbor       */
                                    g[6] += dC / 8.0;   /* Streams to -z neighbor       */

                                    /* Write the updated populations back to the lattice */
                                    lattices[iS]->get(iXt, iYt, iZt).setPopulations(g);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /* appliesTo: Tells Palabos this processor works on the bulk AND envelope
     * (the envelope is a thin layer of ghost cells shared between CPUs). */
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }

    /* clone: Creates a copy of this processor (required by Palabos). */
    virtual run_equilibrium_biotic<T, Descriptor>* clone() const {
        return new run_equilibrium_biotic<T, Descriptor>(*this);
    }

    /* getTypeOfModification: Tells Palabos which lattices this processor modifies.
     * Substrates are modified (staticVariables); the mask is read-only (nothing). */
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iS = 0; iS < subsNum; ++iS) {
            modified[iS] = modif::staticVariables;   /* We modify substrate populations */
        }
        modified[maskLloc] = modif::nothing;         /* We only read the mask */
    }

private:
    plint nx;                          /* Domain size in x (flow direction)       */
    plint subsNum;                     /* Number of substrate species             */
    EquilibriumChemistry<T>& eqChem;   /* Reference to the equilibrium solver     */
    plint solid, bb;                   /* Mask values for solid and bounce-back   */
    plint maskLloc;                    /* Index of the mask lattice in the array  */
};


/*
 * ┌──────────────────────────────────────────────────────────────────────┐
 * │  run_equilibrium_full — Alternative version with DELTA LATTICES      │
 * │                                                                      │
 * │  Instead of updating substrates directly, this processor stores      │
 * │  the concentration CHANGES in separate "delta" lattices.             │
 * │  A second processor (update_equilibriumLattices) then applies them.  │
 * │                                                                      │
 * │  When to use: If you need to accumulate equilibrium changes before   │
 * │  applying them (e.g., for operator-splitting accuracy checks).       │
 * │  For most cases, use run_equilibrium_biotic instead.                 │
 * │                                                                      │
 * │  LATTICE LAYOUT:                                                     │
 * │    lattices[0 ... subsNum-1]           = substrate lattices           │
 * │    lattices[subsNum ... 2*subsNum-1]   = delta lattices (store dC)   │
 * │    lattices[2*subsNum]                 = mask lattice                 │
 * └──────────────────────────────────────────────────────────────────────┘
 */
template<typename T, template<typename U> class Descriptor>
class run_equilibrium_full : public LatticeBoxProcessingFunctional3D<T, Descriptor>
{
public:
    run_equilibrium_full(plint nx_, plint subsNum_, EquilibriumChemistry<T>& eqChem_,
                         plint solid_, plint bb_)
        : nx(nx_), subsNum(subsNum_), eqChem(eqChem_), solid(solid_), bb(bb_),
          dCloc(subsNum_),               /* Delta lattices start at index subsNum   */
          maskLloc(2 * subsNum_)         /* Mask lattice is at index 2*subsNum      */
    {}

    virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
        Dot3D absoluteOffset = lattices[0]->getLocation();

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;

            if (absX > 0 && absX < nx - 1) {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                    for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                        Dot3D maskOffset = computeRelativeDisplacement(*lattices[0], *lattices[maskLloc]);
                        plint mask = util::roundToInt(lattices[maskLloc]->get(
                            iX + maskOffset.x, iY + maskOffset.y, iZ + maskOffset.z).computeDensity());

                        if (mask != solid && mask != bb) {
                            std::vector<Dot3D> vec_offset;
                            for (plint iT = 0; iT <= maskLloc; ++iT) {
                                vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT]));
                            }

                            /* Read current concentrations from substrate lattices */
                            std::vector<T> conc(subsNum);
                            for (plint iS = 0; iS < subsNum; ++iS) {
                                plint iXs = iX + vec_offset[iS].x;
                                plint iYs = iY + vec_offset[iS].y;
                                plint iZs = iZ + vec_offset[iS].z;
                                T c0 = lattices[iS]->get(iXs, iYs, iZs).computeDensity();
                                conc[iS] = std::max(c0, EquilibriumChemistry<T>::MIN_CONC);
                            }

                            /* Solve equilibrium */
                            std::vector<T> eq_conc = eqChem.calculate_species_concentrations(conc);

                            /* Store changes in DELTA lattices (not substrate lattices!) */
                            for (plint iS = 0; iS < subsNum; ++iS) {
                                T dC = eq_conc[iS] - conc[iS];

                                if (std::abs(dC) > thrd) {
                                    Array<T, 7> g;
                                    /* Note: writing to lattice[iS + dCloc] = the delta lattice */
                                    plint iXd = iX + vec_offset[iS + dCloc].x;
                                    plint iYd = iY + vec_offset[iS + dCloc].y;
                                    plint iZd = iZ + vec_offset[iS + dCloc].z;
                                    lattices[iS + dCloc]->get(iXd, iYd, iZd).getPopulations(g);

                                    g[0] += dC / 4.0;
                                    g[1] += dC / 8.0;
                                    g[2] += dC / 8.0;
                                    g[3] += dC / 8.0;
                                    g[4] += dC / 8.0;
                                    g[5] += dC / 8.0;
                                    g[6] += dC / 8.0;

                                    lattices[iS + dCloc]->get(iXd, iYd, iZd).setPopulations(g);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }

    virtual run_equilibrium_full<T, Descriptor>* clone() const {
        return new run_equilibrium_full<T, Descriptor>(*this);
    }

    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iS = 0; iS < subsNum; ++iS) {
            modified[iS] = modif::nothing;              /* Substrates: read-only here */
        }
        for (plint iS = 0; iS < subsNum; ++iS) {
            modified[iS + dCloc] = modif::staticVariables;  /* Deltas: modified */
        }
        modified[maskLloc] = modif::nothing;            /* Mask: read-only */
    }

private:
    plint nx;
    plint subsNum;
    EquilibriumChemistry<T>& eqChem;
    plint solid, bb;
    plint dCloc, maskLloc;
};


/*
 * ┌──────────────────────────────────────────────────────────────────────┐
 * │  update_equilibriumLattices — Apply stored delta changes             │
 * │                                                                      │
 * │  Companion to run_equilibrium_full. Reads the delta lattices and     │
 * │  adds their values to the substrate lattices.                        │
 * │                                                                      │
 * │  LATTICE LAYOUT: same as run_equilibrium_full                        │
 * │    lattices[0 ... subsNum-1]           = substrate lattices           │
 * │    lattices[subsNum ... 2*subsNum-1]   = delta lattices              │
 * │    lattices[2*subsNum]                 = mask lattice                 │
 * └──────────────────────────────────────────────────────────────────────┘
 */
template<typename T, template<typename U> class Descriptor>
class update_equilibriumLattices : public LatticeBoxProcessingFunctional3D<T, Descriptor>
{
public:
    update_equilibriumLattices(plint nx_, plint subsNum_, plint solid_, plint bb_)
        : nx(nx_), subsNum(subsNum_), solid(solid_), bb(bb_),
          dCloc(subsNum_), maskLloc(2 * subsNum_)
    {}

    virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
        Dot3D absoluteOffset = lattices[0]->getLocation();

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;

            if (absX > 0 && absX < nx - 1) {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                    for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                        Dot3D maskOffset = computeRelativeDisplacement(*lattices[0], *lattices[maskLloc]);
                        plint mask = util::roundToInt(lattices[maskLloc]->get(
                            iX + maskOffset.x, iY + maskOffset.y, iZ + maskOffset.z).computeDensity());

                        if (mask != solid && mask != bb) {
                            std::vector<Dot3D> vec_offset;
                            for (plint iT = 0; iT <= maskLloc; ++iT) {
                                vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT]));
                            }

                            for (plint iS = 0; iS < subsNum; ++iS) {
                                /* Read dC from the delta lattice */
                                plint iXd = iX + vec_offset[iS + dCloc].x;
                                plint iYd = iY + vec_offset[iS + dCloc].y;
                                plint iZd = iZ + vec_offset[iS + dCloc].z;
                                T dC = lattices[iS + dCloc]->get(iXd, iYd, iZd).computeDensity();

                                if (std::abs(dC) > thrd) {
                                    /* Apply dC to the substrate lattice */
                                    Array<T, 7> g;
                                    plint iXs = iX + vec_offset[iS].x;
                                    plint iYs = iY + vec_offset[iS].y;
                                    plint iZs = iZ + vec_offset[iS].z;
                                    lattices[iS]->get(iXs, iYs, iZs).getPopulations(g);

                                    g[0] += dC / 4.0;
                                    g[1] += dC / 8.0;
                                    g[2] += dC / 8.0;
                                    g[3] += dC / 8.0;
                                    g[4] += dC / 8.0;
                                    g[5] += dC / 8.0;
                                    g[6] += dC / 8.0;

                                    lattices[iS]->get(iXs, iYs, iZs).setPopulations(g);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }

    virtual update_equilibriumLattices<T, Descriptor>* clone() const {
        return new update_equilibriumLattices<T, Descriptor>(*this);
    }

    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iS = 0; iS < subsNum; ++iS) {
            modified[iS] = modif::staticVariables;     /* Substrates: modified */
        }
        for (plint iS = 0; iS < subsNum; ++iS) {
            modified[iS + dCloc] = modif::nothing;     /* Deltas: read-only here */
        }
        modified[maskLloc] = modif::nothing;           /* Mask: read-only */
    }

private:
    plint nx;
    plint subsNum;
    plint solid, bb;
    plint dCloc, maskLloc;
};


/*
 * ┌──────────────────────────────────────────────────────────────────────┐
 * │  reset_deltaLattices — Zero out delta lattices                       │
 * │                                                                      │
 * │  After applying deltas, they need to be reset to zero for the next   │
 * │  timestep. This processor sets all populations to the "zero density" │
 * │  state of D3Q7: g[0]=-0.25, g[1..6]=-0.125                          │
 * │  (these sum to -1, and computeDensity adds 1, giving density = 0).  │
 * │                                                                      │
 * │  LATTICE LAYOUT:                                                     │
 * │    lattices[0 ... numDelta-1] = delta lattices to reset              │
 * │    lattices[numDelta]         = mask lattice                         │
 * └──────────────────────────────────────────────────────────────────────┘
 */
template<typename T, template<typename U> class Descriptor>
class reset_deltaLattices : public LatticeBoxProcessingFunctional3D<T, Descriptor>
{
public:
    reset_deltaLattices(plint nx_, plint numDelta_, plint solid_, plint bb_)
        : nx(nx_), numDelta(numDelta_), solid(solid_), bb(bb_), maskLloc(numDelta_)
    {}

    virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
        Dot3D absoluteOffset = lattices[0]->getLocation();

        /* D3Q7 populations that give density = 0.
         * In Palabos, density = sum(g) + 1, so for density = 0:
         *   g[0] + g[1] + ... + g[6] = -1
         *   With D3Q7 weights (1/4 and 1/8): g[0] = -1/4, g[1..6] = -1/8 */
        Array<T, 7> g_zero;
        g_zero[0] = -0.25;
        g_zero[1] = g_zero[2] = g_zero[3] = g_zero[4] = g_zero[5] = g_zero[6] = -0.125;

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;

            if (absX > 0 && absX < nx - 1) {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                    for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                        Dot3D maskOffset = computeRelativeDisplacement(*lattices[0], *lattices[maskLloc]);
                        plint mask = util::roundToInt(lattices[maskLloc]->get(
                            iX + maskOffset.x, iY + maskOffset.y, iZ + maskOffset.z).computeDensity());

                        if (mask != solid && mask != bb) {
                            /* Set each delta lattice to zero at this voxel */
                            for (plint iD = 0; iD < numDelta; ++iD) {
                                Dot3D offset = computeRelativeDisplacement(*lattices[0], *lattices[iD]);
                                plint iXd = iX + offset.x;
                                plint iYd = iY + offset.y;
                                plint iZd = iZ + offset.z;
                                lattices[iD]->get(iXd, iYd, iZd).setPopulations(g_zero);
                            }
                        }
                    }
                }
            }
        }
    }

    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }

    virtual reset_deltaLattices<T, Descriptor>* clone() const {
        return new reset_deltaLattices<T, Descriptor>(*this);
    }

    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iD = 0; iD < numDelta; ++iD) {
            modified[iD] = modif::staticVariables;   /* Deltas: modified (reset to zero) */
        }
        modified[maskLloc] = modif::nothing;         /* Mask: read-only */
    }

private:
    plint nx;
    plint numDelta;
    plint solid, bb;
    plint maskLloc;
};


/* ╔════════════════════════════════════════════════════════════════════════════╗
 * ║  SUMMARY OF THIS FILE                                                    ║
 * ╠══════════════════════════════════════════════════════════════════════════╣
 * ║                                                                          ║
 * ║  EquilibriumChemistry<T>           The solver itself                     ║
 * ║    .setSpeciesNames()              Tell it what species exist            ║
 * ║    .setComponentNames()            Tell it which are components          ║
 * ║    .setLogK()                      Tell it the equilibrium constants     ║
 * ║    .setStoichiometryMatrix()       Tell it the stoichiometry             ║
 * ║    .calculate_species_concentrations()  <-- CALL THIS to solve          ║
 * ║    .printStatistics()              Print convergence summary             ║
 * ║                                                                          ║
 * ║  run_equilibrium_biotic            RECOMMENDED LBM processor             ║
 * ║    Lattice layout: [substrates..., mask]                                 ║
 * ║    Directly updates substrate populations                                ║
 * ║                                                                          ║
 * ║  run_equilibrium_full              Alternative with delta lattices        ║
 * ║    Lattice layout: [substrates..., deltas..., mask]                      ║
 * ║    Stores changes in deltas (applied later)                              ║
 * ║                                                                          ║
 * ║  update_equilibriumLattices        Apply stored deltas to substrates     ║
 * ║  reset_deltaLattices               Zero out delta lattices               ║
 * ║                                                                          ║
 * ╚════════════════════════════════════════════════════════════════════════════╝
 */

#endif /* COMPLAB3D_PROCESSORS_PART4_HH */
