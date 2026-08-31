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


/* ╔════════════════════════════════════════════════════════════════════════════╗
 * ║                                                                          ║
 * ║  FILE: complab3d_processors_part4_eqsolver_NR_reserve.hh                ║
 * ║  STATUS: *** DEPRECATED / RESERVE COPY — NOT COMPILED ***               ║
 * ║                                                                          ║
 * ║  Author:      Shahram Asgari                                             ║
 * ║  Advisor:     Dr. Christof Meile                                         ║
 * ║  Laboratory:  Meile Lab, University of Georgia (UGA)                     ║
 * ║                                                                          ║
 * ║  ─────────────────────────────────────────────────────────────────────── ║
 * ║                                                                          ║
 * ║  WHAT THIS FILE IS:                                                      ║
 * ║  This file contains the ORIGINAL Newton-Raphson equilibrium solver       ║
 * ║  that was developed first. It has been REPLACED by the Anderson          ║
 * ║  Acceleration + PCF method (in complab3d_processors_part4_eqsolver.hh)  ║
 * ║  and is kept here as a reserve / reference copy.                         ║
 * ║                                                                          ║
 * ║  ALL CODE IN THIS FILE IS COMMENTED OUT and will not be compiled.        ║
 * ║  To reactivate, remove the comment markers.                              ║
 * ║                                                                          ║
 * ║  ─────────────────────────────────────────────────────────────────────── ║
 * ║                                                                          ║
 * ║  WHY WAS NEWTON-RAPHSON REPLACED?                                        ║
 * ║                                                                          ║
 * ║  Newton-Raphson (NR) is the "textbook" method for solving nonlinear      ║
 * ║  equations. It uses derivatives (the Jacobian matrix) to find the        ║
 * ║  solution, and converges very fast (quadratically) when it works.        ║
 * ║  However, for geochemical equilibrium problems, NR has weaknesses:       ║
 * ║                                                                          ║
 * ║    1. CONVERGENCE RATE: NR converges ~85-90% of the time for stiff      ║
 * ║       geochemical systems. PCF+AA achieves ~99.5%.                       ║
 * ║       [Awada et al. 2025, Carrayrou et al. 2002]                        ║
 * ║                                                                          ║
 * ║    2. DERIVATIVES REQUIRED: NR needs to compute the Jacobian matrix     ║
 * ║       (partial derivatives of residuals with respect to unknowns).      ║
 * ║       PCF is derivative-free — simpler and cheaper per iteration.       ║
 * ║                                                                          ║
 * ║    3. ROBUSTNESS: NR can diverge for stiff systems (large logK          ║
 * ║       spreads) or poor initial guesses. PCF+AA handles these better.    ║
 * ║                                                                          ║
 * ║    4. JACOBIAN SINGULARITY: If the Jacobian is nearly singular          ║
 * ║       (ill-conditioned), NR's linear solve can blow up. PCF+AA's        ║
 * ║       condition monitoring auto-handles this.                            ║
 * ║                                                                          ║
 * ║  REFERENCES:                                                             ║
 * ║  [1] Bethke (2007) - Geochemical and Biogeochemical Reaction Modeling   ║
 * ║  [2] Steefel & MacQuarrie (1996) - Reactive transport formulation       ║
 * ║  [3] Carrayrou et al. (2002) - Comparison of NR vs PCF                  ║
 * ║  [4] Awada et al. (2025) - Anderson + PCF, 99.5% convergence           ║
 * ║                                                                          ║
 * ╚════════════════════════════════════════════════════════════════════════════╝
 *
 *
 * ┌──────────────────────────────────────────────────────────────────────────┐
 * │                                                                          │
 * │  HOW NEWTON-RAPHSON WORKS (for non-specialists):                         │
 * │  ───────────────────────────────────────────────                         │
 * │                                                                          │
 * │  Imagine you're trying to find where a curve crosses zero (the "root").  │
 * │  Newton-Raphson says:                                                    │
 * │                                                                          │
 * │    1. Start with a guess x₀                                              │
 * │    2. Draw the tangent line at x₀ (the tangent = the derivative)         │
 * │    3. Follow the tangent line to where it crosses zero → that's x₁       │
 * │    4. Repeat from x₁                                                     │
 * │                                                                          │
 * │  For our chemistry problem:                                              │
 * │    - The "curve" is the mass balance residual F (should be zero)          │
 * │    - The "tangent" is the Jacobian matrix J (derivatives of F)           │
 * │    - The "step" is: x_new = x - J⁻¹ × F (solve J × step = F, then     │
 * │      subtract step from x)                                               │
 * │                                                                          │
 * │  PROS: When it works, it converges very fast (error squares each step)   │
 * │  CONS: Needs derivatives; can overshoot and diverge; needs line search   │
 * │                                                                          │
 * │  The PCF+AA method (in the active solver file) achieves similar speed    │
 * │  WITHOUT needing derivatives, which is why we switched to it.            │
 * │                                                                          │
 * └──────────────────────────────────────────────────────────────────────────┘
 */


/*
 * Guard against double-inclusion (same concept as the active solver file).
 */
#ifndef COMPLAB3D_PROCESSORS_PART4_NR_RESERVE_HH
#define COMPLAB3D_PROCESSORS_PART4_NR_RESERVE_HH

#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>

using namespace plb;
typedef double T;


/* ══════════════════════════════════════════════════════════════════════════
 *
 *  ██████  EVERYTHING BELOW IS COMMENTED OUT  ██████
 *
 *  The entire Newton-Raphson solver and its LBM processor are wrapped
 *  in block comments.  To reactivate:
 *    1. Remove the "/*" on the line marked "BEGIN COMMENTED-OUT"
 *    2. Remove the "/*" on the line marked "END COMMENTED-OUT"
 *    3. Also uncomment the run_equilibrium_NR processor below
 *
 * ══════════════════════════════════════════════════════════════════════════
 */


/*  <<<< BEGIN COMMENTED-OUT NEWTON-RAPHSON SOLVER >>>>


// ╔════════════════════════════════════════════════════════════════════════╗
// ║  CLASS: EquilibriumChemistry_NR                                       ║
// ║                                                                        ║
// ║  Same purpose as EquilibriumChemistry in the active file, but uses    ║
// ║  Newton-Raphson instead of PCF+AA. See the active file for general   ║
// ║  vocabulary and concepts (species, components, logK, stoichiometry).  ║
// ╚════════════════════════════════════════════════════════════════════════╝

template<typename T>
class EquilibriumChemistry_NR
{
public:
    // Same physical bounds as the PCF+AA version
    static constexpr T MIN_CONC = 1e-30;
    static constexpr T MAX_CONC = 10.0;
    static constexpr T MIN_LOG_C = -30.0;
    static constexpr T MAX_LOG_C = 1.0;

    EquilibriumChemistry_NR()
        : maxIterations(100),       // NR usually converges in fewer iterations (but fails more often)
          tolerance(1e-12),         // Tighter tolerance since NR converges quadratically
          lastConverged(true),
          lastIterations(0), lastResidual(0.0),
          verbose(false),
          totalSolves(0), totalConverged(0), totalDiverged(0)
    {}

    // ════════════════════════════════════════════════════════════════════
    // SETTERS & GETTERS: Same interface as the PCF+AA version.
    // This is intentional — it means you can swap one solver for the
    // other without changing any calling code.
    // ════════════════════════════════════════════════════════════════════

    void setSpeciesNames(const std::vector<std::string>& names) {
        species_names.clear();
        for (const auto& name : names) {
            if (name != "H2O" && name != "h2o") species_names.push_back(name);
        }
    }

    void setComponentNames(const std::vector<std::string>& names) { component_names = names; }
    void setLogK(const std::vector<T>& logK) { logK_values = logK; }
    void setStoichiometryMatrix(const std::vector<std::vector<T>>& S) { stoich_matrix = S; }
    void setMaxIterations(plint maxIter) { maxIterations = maxIter; }
    void setTolerance(T tol) { tolerance = tol; }
    void setVerbose(bool v) { verbose = v; }

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
    plint getTotalSolves() const { return totalSolves; }
    plint getTotalConverged() const { return totalConverged; }
    plint getTotalDiverged() const { return totalDiverged; }
    void resetStatistics() { totalSolves = totalConverged = totalDiverged = 0; }

    // Same species check as PCF+AA version
    bool isEquilibriumSpecies(size_t species_idx) const {
        if (species_idx >= stoich_matrix.size()) return false;
        for (size_t j = 0; j < stoich_matrix[species_idx].size(); ++j) {
            if (std::abs(stoich_matrix[species_idx][j]) > 1e-10) return true;
        }
        return false;
    }

    void setStoichiometryRow(size_t species_idx, const std::vector<T>& row) {
        if (species_idx < stoich_matrix.size() && row.size() == component_names.size()) {
            stoich_matrix[species_idx] = row;
        }
    }


    // ╔════════════════════════════════════════════════════════════════════╗
    // ║  calc_species — MASS ACTION LAW (same as PCF+AA version)          ║
    // ║                                                                    ║
    // ║  [species_i] = 10^( logK[i] + Σⱼ S[i][j] × logC[j] )            ║
    // ║                                                                    ║
    // ║  See the active solver file for the detailed explanation.          ║
    // ╚════════════════════════════════════════════════════════════════════╝
    std::vector<T> calc_species(const std::vector<T>& logC,
                                const std::vector<T>& initial_conc = std::vector<T>()) const {
        size_t ns = species_names.size();
        size_t nc = component_names.size();
        std::vector<T> conc(ns);

        for (size_t i = 0; i < ns; ++i) {
            if (!isEquilibriumSpecies(i)) {
                conc[i] = (i < initial_conc.size() && initial_conc[i] > MIN_CONC)
                          ? std::min(initial_conc[i], MAX_CONC) : MIN_CONC;
                continue;
            }
            T logConc = (i < logK_values.size()) ? logK_values[i] : 0.0;
            for (size_t j = 0; j < nc; ++j) {
                if (j < stoich_matrix[i].size()) {
                    logConc += stoich_matrix[i][j] * logC[j];
                }
            }
            logConc = std::max(std::min(logConc, MAX_LOG_C), MIN_LOG_C);
            conc[i] = std::pow(10.0, logConc);
            if (std::isnan(conc[i]) || std::isinf(conc[i])) {
                conc[i] = (i < initial_conc.size() && initial_conc[i] > MIN_CONC)
                          ? initial_conc[i] : MIN_CONC;
            }
            conc[i] = std::max(std::min(conc[i], MAX_CONC), MIN_CONC);
        }
        return conc;
    }

    // Same as PCF+AA version
    std::vector<T> calc_component_totals(const std::vector<T>& species_conc) const {
        size_t nc = component_names.size();
        size_t ns = species_names.size();
        std::vector<T> T_total(nc, 0.0);
        for (size_t j = 0; j < nc; ++j) {
            for (size_t i = 0; i < ns && i < species_conc.size(); ++i) {
                if (isEquilibriumSpecies(i) && j < stoich_matrix[i].size()) {
                    T conc_safe = std::max(std::min(species_conc[i], MAX_CONC), MIN_CONC);
                    T_total[j] += stoich_matrix[i][j] * conc_safe;
                }
            }
            T_total[j] = std::max(T_total[j], MIN_CONC);
        }
        return T_total;
    }


    // ╔════════════════════════════════════════════════════════════════════╗
    // ║                                                                    ║
    // ║  NR-SPECIFIC: MASS BALANCE RESIDUAL                                ║
    // ║  ══════════════════════════════════                                ║
    // ║                                                                    ║
    // ║  F[j] = Σᵢ S[i][j] × [species_i] − T_total[j]                   ║
    // ║                                                                    ║
    // ║  At equilibrium, F = 0 (all mass balances satisfied).              ║
    // ║                                                                    ║
    // ║  In plain English: "For each component, add up how much of it     ║
    // ║  exists across all species, then subtract the target total.        ║
    // ║  The difference is how far off we are."                            ║
    // ║                                                                    ║
    // ║  NOTE: The PCF+AA version doesn't use this function at all.       ║
    // ║  PCF computes its own "residual" using the ratio of reactive      ║
    // ║  and product sums (see pcf_residual in the active file).          ║
    // ╚════════════════════════════════════════════════════════════════════╝
    std::vector<T> calc_residual(const std::vector<T>& logC,
                                  const std::vector<T>& T_total,
                                  const std::vector<T>& initial_conc) const {
        // First, compute all species concentrations from the current guess
        std::vector<T> c = calc_species(logC, initial_conc);
        size_t nc = component_names.size();
        std::vector<T> F(nc);

        for (size_t j = 0; j < nc; ++j) {
            F[j] = -T_total[j];   // Start with −T_total (the target, with negative sign)
            for (size_t i = 0; i < species_names.size(); ++i) {
                if (isEquilibriumSpecies(i) && j < stoich_matrix[i].size()) {
                    F[j] += stoich_matrix[i][j] * c[i];   // Add S^T × c
                }
            }
            // F[j] = (S^T × c)[j] − T_total[j]
            // When F = 0, mass balance is satisfied ✓
        }
        return F;
    }


    // ╔════════════════════════════════════════════════════════════════════╗
    // ║                                                                    ║
    // ║  NR-SPECIFIC: THE JACOBIAN MATRIX                                  ║
    // ║  ════════════════════════════════                                  ║
    // ║                                                                    ║
    // ║  The Jacobian is a matrix of partial derivatives:                   ║
    // ║    J[j][k] = ∂F[j] / ∂(logC[k])                                  ║
    // ║                                                                    ║
    // ║  In plain English: "How much does residual j change when we         ║
    // ║  slightly change log-concentration k?"                             ║
    // ║                                                                    ║
    // ║  The formula (derived using the chain rule from calculus):          ║
    // ║    J[j][k] = ln(10) × Σᵢ S[i][j] × S[i][k] × [species_i]       ║
    // ║                                                                    ║
    // ║  The ln(10) ≈ 2.303 comes from differentiating 10^x:              ║
    // ║    d/dx (10^x) = ln(10) × 10^x                                    ║
    // ║                                                                    ║
    // ║  In matrix form: J = ln(10) × S^T × diag(c) × S                  ║
    // ║  This matrix is always symmetric (J[j][k] = J[k][j]) and          ║
    // ║  positive semi-definite (good numerical properties).               ║
    // ║                                                                    ║
    // ║  NOTE: This is the main DISADVANTAGE of NR — computing this       ║
    // ║  Jacobian requires O(nc² × ns) operations, and then we need       ║
    // ║  O(nc³) to solve the linear system. PCF avoids this entirely.     ║
    // ╚════════════════════════════════════════════════════════════════════╝
    std::vector<std::vector<T>> calc_jacobian(const std::vector<T>& logC,
                                               const std::vector<T>& initial_conc) const {
        std::vector<T> c = calc_species(logC, initial_conc);
        size_t nc = component_names.size();
        T ln10 = std::log(10.0);   // ≈ 2.302585

        // Create an nc × nc matrix of zeros
        std::vector<std::vector<T>> J(nc, std::vector<T>(nc, 0.0));

        // Fill in J[j][k] = ln(10) × Σᵢ S[i][j] × S[i][k] × c[i]
        for (size_t j = 0; j < nc; ++j) {
            for (size_t k = 0; k < nc; ++k) {
                for (size_t i = 0; i < species_names.size(); ++i) {
                    if (isEquilibriumSpecies(i) &&
                        j < stoich_matrix[i].size() && k < stoich_matrix[i].size()) {
                        J[j][k] += ln10 * stoich_matrix[i][j] * stoich_matrix[i][k] * c[i];
                    }
                }
            }
        }
        return J;
    }


    // Euclidean norm (same as PCF+AA version)
    T norm(const std::vector<T>& v) const {
        T s = 0;
        for (T x : v) s += x * x;
        return std::sqrt(s);
    }


    // ╔════════════════════════════════════════════════════════════════════╗
    // ║                                                                    ║
    // ║  NR-SPECIFIC: GAUSSIAN ELIMINATION                                 ║
    // ║  ═════════════════════════════════                                 ║
    // ║                                                                    ║
    // ║  Solves the linear equation  A × x = b  for x.                    ║
    // ║                                                                    ║
    // ║  In the NR context: A is the Jacobian J, b is −F (negative        ║
    // ║  residual), and x is the Newton step delta.                        ║
    // ║                                                                    ║
    // ║  HOW IT WORKS (two phases):                                        ║
    // ║                                                                    ║
    // ║  Phase 1 — FORWARD ELIMINATION:                                    ║
    // ║    Turn the matrix into an upper triangle (zeros below diagonal)   ║
    // ║    by subtracting multiples of each row from the rows below it.    ║
    // ║    "Partial pivoting" means we first swap rows to put the          ║
    // ║    largest element on the diagonal, which reduces rounding errors. ║
    // ║                                                                    ║
    // ║  Phase 2 — BACK SUBSTITUTION:                                      ║
    // ║    Start from the last row (1 unknown) and work upward,            ║
    // ║    plugging known values into earlier rows.                        ║
    // ║                                                                    ║
    // ║  NOTE: The PCF+AA solver doesn't need this at all!                ║
    // ║  It uses QR decomposition (a different, more stable approach)     ║
    // ║  for a smaller problem (the Anderson acceleration weights).        ║
    // ╚════════════════════════════════════════════════════════════════════╝
    std::vector<T> gaussian_elimination(std::vector<std::vector<T>> A,
                                         std::vector<T> b, bool& success) const {
        success = true;
        size_t n = A.size();

        if (n == 0 || A[0].size() != n || b.size() != n) {
            success = false;
            return std::vector<T>(n, 0.0);
        }

        const T SINGULAR_THRESHOLD = 1e-30;

        // ── Phase 1: Forward elimination with partial pivoting ──
        for (size_t i = 0; i < n; ++i) {
            // Find the row below (including i) with the largest |A[k][i]|
            // and swap it to position i. This is "partial pivoting."
            size_t maxRow = i;
            T maxVal = std::abs(A[i][i]);
            for (size_t k = i + 1; k < n; ++k) {
                if (std::abs(A[k][i]) > maxVal) {
                    maxVal = std::abs(A[k][i]);
                    maxRow = k;
                }
            }

            // If the largest pivot is essentially zero, the matrix is singular
            // (no unique solution exists). This is one reason NR can fail.
            if (maxVal < SINGULAR_THRESHOLD) {
                success = false;
                return std::vector<T>(n, 0.0);
            }

            // Swap rows i and maxRow
            if (maxRow != i) {
                std::swap(A[i], A[maxRow]);
                std::swap(b[i], b[maxRow]);
            }

            // Eliminate all entries below the pivot
            for (size_t k = i + 1; k < n; ++k) {
                T factor = A[k][i] / A[i][i];
                for (size_t j = i; j < n; ++j) {
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }

        // ── Phase 2: Back substitution ──
        std::vector<T> x(n, 0.0);
        for (int i = n - 1; i >= 0; --i) {
            if (std::abs(A[i][i]) < SINGULAR_THRESHOLD) {
                success = false;
                return std::vector<T>(n, 0.0);
            }
            x[i] = b[i];
            for (size_t j = i + 1; j < n; ++j) {
                x[i] -= A[i][j] * x[j];
            }
            x[i] /= A[i][i];

            if (std::isnan(x[i]) || std::isinf(x[i])) {
                success = false;
                return std::vector<T>(n, 0.0);
            }
        }
        return x;
    }


    // ╔════════════════════════════════════════════════════════════════════╗
    // ║                                                                    ║
    // ║  THE MAIN NR SOLVER                                                ║
    // ║  ══════════════════                                                ║
    // ║                                                                    ║
    // ║  ALGORITHM:                                                        ║
    // ║    1. Start with initial guess logC                                ║
    // ║    2. Compute residual F = S^T × c(logC) − T_total                ║
    // ║    3. Is ||F|| < tolerance? → DONE (converged)                    ║
    // ║    4. Compute Jacobian J                                          ║
    // ║    5. Solve J × delta = −F (using Gaussian elimination)           ║
    // ║    6. LINE SEARCH: try step sizes α = 1, 0.5, 0.25, ...           ║
    // ║       until ||F(logC + α×delta)|| < ||F(logC)||                   ║
    // ║       This prevents overshooting.                                  ║
    // ║    7. Update: logC = logC + α × delta                             ║
    // ║    8. Go to step 2                                                ║
    // ║                                                                    ║
    // ║  THE LINE SEARCH (step 6) is critical:                             ║
    // ║  Without it, NR can take a huge step that actually makes things    ║
    // ║  worse. The line search halves the step size repeatedly until      ║
    // ║  the residual actually decreases. This is called "backtracking     ║
    // ║  with Armijo condition."                                           ║
    // ║                                                                    ║
    // ║  COMPARE WITH PCF+AA:                                              ║
    // ║  - NR needs steps 4-6 (Jacobian + solve + line search) = expensive║
    // ║  - PCF+AA replaces all that with one PCF residual + AA weights    ║
    // ║  - PCF+AA is simpler, cheaper, and more robust                    ║
    // ╚════════════════════════════════════════════════════════════════════╝
    std::vector<T> solve_equilibrium_NR(const std::vector<T>& initial_species_conc,
                                         const std::vector<T>& T_total) {
        lastConverged = false;
        lastIterations = 0;
        lastResidual = 0.0;
        totalSolves++;

        size_t nc = component_names.size();
        if (nc == 0) {
            lastConverged = true;
            totalConverged++;
            return initial_species_conc;
        }

        // ── Step 1: Build initial guess ──
        // (Same approach as PCF+AA: find component concentrations, take log10)
        std::vector<T> logC(nc);
        for (size_t j = 0; j < nc; ++j) {
            auto it = std::find(species_names.begin(), species_names.end(), component_names[j]);
            T comp_conc = 1e-7;
            if (it != species_names.end()) {
                size_t idx = std::distance(species_names.begin(), it);
                if (idx < initial_species_conc.size()) {
                    comp_conc = initial_species_conc[idx];
                }
            }
            comp_conc = std::max(std::min(comp_conc, MAX_CONC), MIN_CONC);
            logC[j] = std::log10(comp_conc);
        }

        // ── Steps 2-8: Iterate ──
        for (plint iter = 0; iter < maxIterations; ++iter) {
            lastIterations = iter + 1;

            // Step 2: Compute residual
            std::vector<T> F = calc_residual(logC, T_total, initial_species_conc);
            T Fnorm = norm(F);
            lastResidual = Fnorm;

            // Step 3: Check convergence
            if (Fnorm < tolerance) {
                lastConverged = true;
                totalConverged++;
                return calc_species(logC, initial_species_conc);
            }

            // Step 4: Compute Jacobian (THIS IS WHAT PCF AVOIDS)
            std::vector<std::vector<T>> J = calc_jacobian(logC, initial_species_conc);

            // Step 5: Solve J × delta = −F
            bool solveSuccess;
            std::vector<T> negF(F.size());
            for (size_t i = 0; i < F.size(); ++i) negF[i] = -F[i];
            std::vector<T> delta = gaussian_elimination(J, negF, solveSuccess);

            if (!solveSuccess) {
                // Jacobian is singular — NR can't proceed.
                // This is one of the failure modes that PCF+AA avoids.
                totalDiverged++;
                return calc_species(logC, initial_species_conc);
            }

            // Step 6: LINE SEARCH (Armijo backtracking)
            // Try full step first (α=1). If residual doesn't improve enough,
            // halve α and try again. Repeat up to 40 times.
            T alpha = 1.0;
            std::vector<T> logC_new(nc);
            T Fnorm_new;

            for (int ls = 0; ls < 40; ++ls) {
                // Trial point: logC_new = logC + α × delta
                for (size_t j = 0; j < nc; ++j) {
                    logC_new[j] = logC[j] + alpha * delta[j];
                    logC_new[j] = std::max(std::min(logC_new[j], MAX_LOG_C), MIN_LOG_C);
                }

                // Check: did the residual decrease enough?
                std::vector<T> F_new = calc_residual(logC_new, T_total, initial_species_conc);
                Fnorm_new = norm(F_new);

                // Armijo condition: ||F_new|| < ||F|| × (1 − 0.0001 × α)
                // Or accept if α is already tiny (step size too small to matter)
                if (Fnorm_new < Fnorm * (1.0 - 0.0001 * alpha) || alpha < 1e-15) {
                    break;
                }
                alpha *= 0.5;   // Halve step size and try again
            }

            // Step 7: Accept the update
            logC = logC_new;

            // Check for stalled convergence
            T step = alpha * norm(delta);
            if (step < tolerance * 1e-3 && Fnorm_new < tolerance * 100) {
                lastConverged = true;
                totalConverged++;
                return calc_species(logC, initial_species_conc);
            }
        }

        // Max iterations without convergence
        totalDiverged++;
        return calc_species(logC, initial_species_conc);
    }


    // Main entry point (same interface as PCF+AA)
    std::vector<T> calculate_species_concentrations(const std::vector<T>& initial_conc) {
        size_t nc = component_names.size();
        size_t ns = species_names.size();
        if (nc == 0 || ns == 0) return initial_conc;

        std::vector<T> T_total = calc_component_totals(initial_conc);
        std::vector<T> result = solve_equilibrium_NR(initial_conc, T_total);

        for (size_t i = 0; i < result.size(); ++i) {
            if (std::isnan(result[i]) || std::isinf(result[i])) {
                result[i] = (i < initial_conc.size() && initial_conc[i] > MIN_CONC)
                            ? std::min(initial_conc[i], MAX_CONC) : MIN_CONC;
            }
            result[i] = std::max(std::min(result[i], MAX_CONC), MIN_CONC);
        }
        return result;
    }


    void printStatistics() const {
        pcout << "╔══════════════════════════════════════════════════════════════════╗\n";
        pcout << "║  NEWTON-RAPHSON EQUILIBRIUM SOLVER STATISTICS (RESERVE)          ║\n";
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

private:
    std::vector<std::string> species_names, component_names;
    std::vector<T> logK_values;
    std::vector<std::vector<T>> stoich_matrix;
    plint maxIterations;
    T tolerance;
    mutable bool lastConverged;
    mutable plint lastIterations;
    mutable T lastResidual;
    bool verbose;
    mutable plint totalSolves, totalConverged, totalDiverged;
};


<<<< END COMMENTED-OUT NEWTON-RAPHSON SOLVER >>>> */



/* ══════════════════════════════════════════════════════════════════════════
 *  NR-BASED LBM PROCESSOR (also commented out)
 * ══════════════════════════════════════════════════════════════════════════
 */


/* <<<< BEGIN COMMENTED-OUT NR PROCESSOR >>>>

// ┌──────────────────────────────────────────────────────────────────────┐
// │  run_equilibrium_NR — LBM data processor using Newton-Raphson        │
// │                                                                      │
// │  Same structure as run_equilibrium_biotic in the active file.        │
// │  Loops over every pore voxel, reads concentrations, calls NR solver, │
// │  writes changes back to lattice populations.                         │
// │                                                                      │
// │  LATTICE LAYOUT:                                                     │
// │    lattices[0 ... subsNum-1] = substrate lattices                    │
// │    lattices[subsNum]         = mask lattice                          │
// └──────────────────────────────────────────────────────────────────────┘

template<typename T, template<typename U> class Descriptor>
class run_equilibrium_NR : public LatticeBoxProcessingFunctional3D<T, Descriptor>
{
public:
    run_equilibrium_NR(plint nx_, plint subsNum_, EquilibriumChemistry_NR<T>& eqChem_,
                       plint solid_, plint bb_)
        : nx(nx_), subsNum(subsNum_), eqChem(eqChem_), solid(solid_), bb(bb_),
          maskLloc(subsNum_)
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

                            // Read concentrations
                            std::vector<T> conc(subsNum);
                            for (plint iS = 0; iS < subsNum; ++iS) {
                                plint iXs = iX + vec_offset[iS].x;
                                plint iYs = iY + vec_offset[iS].y;
                                plint iZs = iZ + vec_offset[iS].z;
                                T c0 = lattices[iS]->get(iXs, iYs, iZs).computeDensity();
                                conc[iS] = std::max(c0, EquilibriumChemistry_NR<T>::MIN_CONC);
                            }

                            // Solve equilibrium using Newton-Raphson
                            std::vector<T> eq_conc = eqChem.calculate_species_concentrations(conc);

                            // Only apply if converged
                            if (eqChem.didConverge()) {
                                for (plint iS = 0; iS < subsNum; ++iS) {
                                    T dC = eq_conc[iS] - conc[iS];
                                    if (std::abs(dC) > thrd) {
                                        Array<T, 7> g;
                                        plint iXt = iX + vec_offset[iS].x;
                                        plint iYt = iY + vec_offset[iS].y;
                                        plint iZt = iZ + vec_offset[iS].z;
                                        lattices[iS]->get(iXt, iYt, iZt).getPopulations(g);

                                        // D3Q7 weights: 1/4 rest, 1/8 each direction
                                        g[0] += dC / 4.0;
                                        g[1] += dC / 8.0;
                                        g[2] += dC / 8.0;
                                        g[3] += dC / 8.0;
                                        g[4] += dC / 8.0;
                                        g[5] += dC / 8.0;
                                        g[6] += dC / 8.0;

                                        lattices[iS]->get(iXt, iYt, iZt).setPopulations(g);
                                    }
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

    virtual run_equilibrium_NR<T, Descriptor>* clone() const {
        return new run_equilibrium_NR<T, Descriptor>(*this);
    }

    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iS = 0; iS < subsNum; ++iS) {
            modified[iS] = modif::staticVariables;
        }
        modified[maskLloc] = modif::nothing;
    }

private:
    plint nx;
    plint subsNum;
    EquilibriumChemistry_NR<T>& eqChem;
    plint solid, bb;
    plint maskLloc;
};

<<<< END COMMENTED-OUT NR PROCESSOR >>>> */



/* ╔════════════════════════════════════════════════════════════════════════════╗
 * ║                                                                          ║
 * ║  COMPARISON TABLE: Newton-Raphson vs Anderson Acceleration + PCF         ║
 * ║                                                                          ║
 * ╠══════════════════════════════════════════════════════════════════════════╣
 * ║                                                                          ║
 * ║  Feature                   │ Newton-Raphson       │ PCF + Anderson Accel. ║
 * ║  ──────────────────────────┼──────────────────────┼──────────────────── ║
 * ║  Convergence rate          │ ~85-90%              │ ~99.5%               ║
 * ║  Requires derivatives?     │ YES (Jacobian)       │ NO (derivative-free) ║
 * ║  Cost per iteration        │ O(nc²×ns) Jacobian   │ O(nc×ns) PCF residual║
 * ║                            │ + O(nc³) linear solve│ + O(m²×nc) QR        ║
 * ║  Convergence speed         │ Quadratic (locally)  │ Superlinear with AA  ║
 * ║  Robustness                │ Needs line search    │ Very robust           ║
 * ║  Ill-conditioning handling │ Can fail silently     │ Auto-reduces depth   ║
 * ║  Implementation complexity │ Moderate             │ Simpler              ║
 * ║  Memory                    │ O(nc²) for Jacobian  │ O(m×nc) for history  ║
 * ║                                                                          ║
 * ║  nc = number of components, ns = number of species, m = Anderson depth   ║
 * ║                                                                          ║
 * ╚════════════════════════════════════════════════════════════════════════════╝
 */


#endif /* COMPLAB3D_PROCESSORS_PART4_NR_RESERVE_HH */
