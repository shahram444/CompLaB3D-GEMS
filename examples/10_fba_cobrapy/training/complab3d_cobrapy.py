"""
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

===================================================================================================
 Python half of the optional CompLaB3D COBRApy bridge (C++ half: src/complab3d_pythonAPI.hh).
 Ported from the 2-D src/complab_cobrapy.py.

 Requires cobrapy >= 0.13 (optlang-based objective interface, Reaction.bounds). See the note on
 the objective assignment in feed_cobraModel() -- the deprecated Reaction.objective_coefficient
 attribute used by the 2-D version is NOT used here, which is what sets that minimum version.

===================================================================================================
 FLUX-VECTOR CONTRACT (C++ <-> complab3d_cobrapy.py)
 ---------------------------------------------------------------------------------------------
 THIS BLOCK IS DUPLICATED VERBATIM IN complab3d_pythonAPI.hh. If you change one, change both.

 Argument vector handed to run_FBA() by optimize_cobrapy() (all values are Python floats, models
 appended last as opaque handles):

     args[0]                                    = M            number of microbes present in this voxel
     args[1 + i]                                = n_i          number of USED exchange reactions of microbe i
     args[1 + M + off_i + k]                    = loc[i][k]    reaction INDEX inside microbe i's own model
     args[1 + M + N + off_i + k]                = lb[i][k]     lower bound for reaction loc[i][k]
     args[1 + M + 2N + off_i + k]               = ub[i][k]     upper bound for reaction loc[i][k]
     args[len(args) - M + i]                    = model handle of microbe i (returned by feed_cobraModel)

     with N = sum_i n_i,  off_i = sum_{j<i} n_j,  k = 0 .. n_i-1.

 Return value of run_FBA(): the list [grate, flux, status] with

     grate[i]     float   objective value of microbe i, or 0.0 when the LP did not solve
     flux[i]      list of EXACTLY n_i floats.  flux[i][k] is the flux of reaction loc[i][k].
                  ==> flux[i] is indexed by the PACKED index k, in the SAME ORDER as loc[i] was
                      supplied. It is NOT indexed by the global substrate id iS of the calling
                      data processor. The C++ side calls expandFluxToSubstrates() to convert.
     status[i]    str     solver status of microbe i ('optimal' when the LP solved, anything else
                          is treated as failure and hard-zeroed by the C++ side)

 Why this is spelled out: in the 2-D code the data processor read vec2_flux[iB][iS] with iS running
 over ALL substrates while this function returned only the microbe's used exchange reactions, so the
 fluxes were silently attributed to the wrong solutes (and read out of bounds) as soon as a microbe
 did not exchange every substrate. Length and order are now pinned here and verified in C++.
===================================================================================================
"""

import cobra
from cobra import Model, Reaction, Metabolite

# Status string that the C++ side accepts as "solution usable".
STATUS_OPTIMAL = "optimal"


def feed_cobraModel(*args):
    """Build one cobra model from the dense arrays shipped by prep_cobrapy().

    Wire format (unchanged from 2-D):
        args[0]                              = nmets
        args[1]                              = nrxns
        args[2 .. 2+nmets*nrxns-1]           = S, metabolite-major: S[i*nrxns + j] = metabolite i in reaction j
        args[2+ nmets   *nrxns .. +nrxns-1]  = lower bounds
        args[2+(nmets+1)*nrxns .. +nrxns-1]  = upper bounds
        args[2+(nmets+2)*nrxns .. +nrxns-1]  = objective coefficients
        args[-1]                             = objective sense: -1 maximize, +1 minimize
                                               ([FIX-3D] new; <objective_direction> used to reach
                                               GLPK but not COBRApy, so `minimize` was silently
                                               maximised here)

    Metabolites are named "0".."nmets-1" and reactions "0".."nrxns-1", so the reaction indices that
    the data processor stores in vec_EX_loc address model.reactions positionally, as in 2-D.

    Returns a one-element list [model]; the C++ side keeps it as an opaque handle and hands it back
    to run_FBA(), which unwraps it with [0]. The wrapper is kept for wire compatibility with 2-D.
    """
    nmets = int(args[0])
    nrxns = int(args[1])

    S = args[2:2 + nmets * nrxns]
    lb = args[2 + nmets * nrxns: 2 + (nmets + 1) * nrxns]
    ub = args[2 + (nmets + 1) * nrxns: 2 + (nmets + 2) * nrxns]
    c = args[2 + (nmets + 2) * nrxns: 2 + (nmets + 3) * nrxns]
    # [FIX-3D] trailing objective sense. Tolerate the old, shorter wire format.
    sense = int(args[-1]) if len(args) > 2 + (nmets + 3) * nrxns else -1

    model = Model('0')

    # ---------------------------------------------------------------------------------------------
    # PERFORMANCE: the 2-D version created a fresh Metabolite object inside the reaction loop and
    # called reaction.add_metabolites({met: coeff}) once per (metabolite, reaction) pair -- including
    # the ~99% of pairs whose stoichiometric coefficient is exactly zero. That is O(nmets*nrxns)
    # Python-level calls plus nmets*nrxns throw-away Metabolite objects, and every zero coefficient
    # still ends up in the model's sparse matrix as an explicit structural zero.
    # Here the Metabolite objects are created ONCE up front, each reaction gets a single
    # add_metabolites() call with a sparse dict, and zero coefficients are skipped entirely. For a
    # genome-scale network (nmets ~ 1000, nrxns ~ 2000, density ~0.2%) this turns ~2e6 calls into
    # ~2e3 calls handling ~4e3 nonzeros: three orders of magnitude fewer operations, and model
    # construction drops from minutes to well under a second.
    # ---------------------------------------------------------------------------------------------
    metabolites = [Metabolite(str(i)) for i in range(nmets)]

    reactions = []
    for j in range(nrxns):
        reaction = Reaction(str(j))
        # Reaction.bounds sets both bounds atomically; assigning lower_bound first can raise when the
        # new lower bound exceeds the still-default upper bound.
        reaction.bounds = (float(lb[j]), float(ub[j]))
        sparse = {}
        for i in range(nmets):
            coeff = S[i * nrxns + j]
            if coeff != 0.0:
                sparse[metabolites[i]] = float(coeff)
        if sparse:
            reaction.add_metabolites(sparse)
        reactions.append(reaction)

    # One batched call instead of nrxns separate model.add_reactions([reaction]) calls.
    model.add_reactions(reactions)

    # Objective. MODERN COBRAPY FORM: assigning a {Reaction: coefficient} dict to model.objective
    # goes through the optlang interface. The 2-D code set the deprecated
    # model.reactions[j].objective_coefficient attribute instead; that attribute has been removed in
    # recent cobrapy. Using the dict form is what makes cobrapy >= 0.13 the minimum version here.
    objective = {}
    for j in range(nrxns):
        if c[j] != 0.0:
            objective[model.reactions[j]] = float(c[j])
    if objective:
        model.objective = objective
        # [FIX-3D] Apply the requested direction. Must come AFTER the assignment above,
        # which replaces the objective object and resets its direction to the default.
        model.objective.direction = 'min' if sense > 0 else 'max'
    else:
        print("complab3d_cobrapy WARNING: no non-zero objective coefficient; the model has a null objective.")

    return [model]


def run_FBA(*args):
    """Solve one FBA per microbe present in the current voxel.

    See the FLUX-VECTOR CONTRACT block at the top of this file for the exact argument layout and for
    the length/ordering guarantee of the returned flux vectors.

    Only the reactions listed in loc[i] are re-bounded on every call; all other bounds persist from
    feed_cobraModel(). The model objects are reused across voxels and time steps by design.
    """
    num_of_microbes = int(args[0])
    total_length = len(args)

    num_of_subs = [int(args[1 + i]) for i in range(num_of_microbes)]
    total_num_of_subs = sum(num_of_subs)

    model = [None] * num_of_microbes
    loc = [None] * num_of_microbes
    lb = [None] * num_of_microbes
    ub = [None] * num_of_microbes

    off = 0
    base_loc = 1 + num_of_microbes
    base_lb = base_loc + total_num_of_subs
    base_ub = base_loc + 2 * total_num_of_subs
    for i in range(num_of_microbes):
        n_i = num_of_subs[i]
        model[i] = args[total_length - num_of_microbes + i][0]   # unwrap the [model] handle
        loc[i] = [int(args[base_loc + off + k]) for k in range(n_i)]
        lb[i] = [float(args[base_lb + off + k]) for k in range(n_i)]
        ub[i] = [float(args[base_ub + off + k]) for k in range(n_i)]
        off += n_i

    # Apply the voxel-local exchange bounds.
    for i in range(num_of_microbes):
        reactions = model[i].reactions
        for k in range(num_of_subs[i]):
            # Atomic assignment, and robust when the new lower bound is above the old upper bound.
            reactions[loc[i][k]].bounds = (lb[i][k], ub[i][k])

    grate = [0.0] * num_of_microbes
    flux = [None] * num_of_microbes
    status = [""] * num_of_microbes

    for i in range(num_of_microbes):
        n_i = num_of_subs[i]
        try:
            solution = model[i].optimize()
        except Exception as err:                      # solver blow-up must not kill the simulation
            print("complab3d_cobrapy ERROR: optimize() raised for microbe %d: %s" % (i, err))
            grate[i] = 0.0
            flux[i] = [0.0] * n_i
            status[i] = "solver_error"
            continue

        status[i] = str(solution.status)
        # An infeasible / unbounded LP gives objective_value = None (and NaN fluxes). Never let that
        # reach the C++ side: PyFloat_AsDouble(None) returns -1.0 and sets a Python exception, which
        # the 2-D bridge never checked -- an infeasible voxel silently grew at -1.0 h^-1.
        if status[i] != STATUS_OPTIMAL or solution.objective_value is None:
            grate[i] = 0.0
            flux[i] = [0.0] * n_i
            continue

        grate[i] = float(solution.objective_value)
        # CONTRACT: exactly one flux per entry of loc[i], in the same order.
        # solution.fluxes is a pandas Series indexed by reaction ID, so it is looked up by ID here.
        # The 2-D code did solution.fluxes[loc[i][j]] with an *integer* on a string-labelled index,
        # which only worked through pandas' positional fallback -- removed in pandas 2.x.
        fluxes = solution.fluxes
        reactions = model[i].reactions
        flux[i] = [float(fluxes[reactions[l].id]) for l in loc[i]]

    return [grate, flux, status]


if __name__ == "__main__":
    # ---------------------------------------------------------------------------------------------
    # Smoke test: builds a tiny 3-reaction network through the same wire format the C++ side uses and
    # solves it through run_FBA(), so this file can be checked without compiling CompLaB3D.
    #
    #   metabolites: A (index 0), B (index 1)
    #   reactions:   R0: -> A        S[A][R0] = +1
    #                R1: A -> B      S[A][R1] = -1, S[B][R1] = +1
    #                R2: B ->        S[B][R2] = -1     (objective)
    # ---------------------------------------------------------------------------------------------
    nmets, nrxns = 2, 3
    S_dense = [1.0, -1.0, 0.0,      # metabolite A over R0, R1, R2
               0.0, 1.0, -1.0]      # metabolite B over R0, R1, R2
    lb0 = [0.0, 0.0, 0.0]
    ub0 = [10.0, 1000.0, 1000.0]
    c0 = [0.0, 0.0, 1.0]            # maximize R2

    def build():
        return feed_cobraModel(*([float(nmets), float(nrxns)] + S_dense + lb0 + ub0 + c0))

    failures = 0

    # --- test 1: feasible, loc = [R0, R2], uptake capped at 7 --------------------------------------
    handle = build()
    args = [1.0,            # M
            2.0,            # n_0
            0.0, 2.0,       # loc[0] = R0, R2
            0.0, 0.0,       # lb[0]
            7.0, 1000.0,    # ub[0]
            handle]
    grate, flux, status = run_FBA(*args)
    print("test 1 (feasible): status = %s, growth = %.6f, flux = %s" % (status[0], grate[0], flux[0]))
    if status[0] != STATUS_OPTIMAL or abs(grate[0] - 7.0) > 1e-6:
        print("  FAILED: expected an optimal solution with objective 7.0")
        failures += 1
    if len(flux[0]) != 2:
        print("  FAILED: flux vector must have exactly len(loc) = 2 entries (contract violation)")
        failures += 1
    elif abs(flux[0][0] - 7.0) > 1e-6 or abs(flux[0][1] - 7.0) > 1e-6:
        print("  FAILED: expected flux [7.0, 7.0] in loc order")
        failures += 1

    # --- test 2: infeasible (forced influx of 5 with the outlet shut) ------------------------------
    handle = build()
    args = [1.0,
            2.0,
            0.0, 2.0,       # loc[0] = R0, R2
            5.0, 0.0,       # lb[0]: R0 forced to 5, R2 forced to 0
            5.0, 0.0,       # ub[0]
            handle]
    grate, flux, status = run_FBA(*args)
    print("test 2 (infeasible): status = %s, growth = %.6f, flux = %s" % (status[0], grate[0], flux[0]))
    if status[0] == STATUS_OPTIMAL:
        print("  FAILED: this LP has no feasible point, status should not be 'optimal'")
        failures += 1
    if grate[0] != 0.0 or any(f != 0.0 for f in flux[0]):
        print("  FAILED: an infeasible LP must be reported as zero growth and zero fluxes, never None/NaN")
        failures += 1

    print("cobrapy version: %s" % cobra.__version__)
    print("SELF-TEST %s" % ("PASSED" if failures == 0 else "FAILED (%d problem(s))" % failures))
