#!/usr/bin/env python3
"""
/* This file is a part of the CompLaB program.
 *
 * The CompLaB3D software (3-D pore-scale extension) is developed since 2024
 * by the University of Georgia (United States, Meile Lab, Department of
 * Marine Sciences). The original 2-D CompLaB v1.0 was a collaboration of
 * the University of Georgia and Chungnam National University (South Korea).
 *
 * Contact:
 * Shahram Asgari, Christof Meile
 * Department of Marine Sciences (Meile Lab)
 * University of Georgia, Athens, GA 30602, USA
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

===============================================================================
STEP 1 OF 3 -- BUILD THE TRAINING SET FOR A CompLaB SURROGATE MODEL
===============================================================================

WHAT A SURROGATE IS FOR

  Flux balance analysis solves one linear program per microbe, per active
  voxel, per time step.  In a 3-D domain of a few million voxels that cost
  dominates everything else in the run.  A surrogate replaces the linear
  program with a small fitted function:

        substrate uptake rates  ->  predicted specific growth rate

  evaluated in a few hundred multiplications.  Speed-ups of 100-1000x versus
  in-line FBA are typical.  The price is that the fit is only valid inside the
  range of uptake rates it was trained on -- a neural network extrapolates
  badly and will confidently return nonsense outside that box.  That is why
  this script records the sampled range and writes it into the CSV header,
  and why trainSurrogate.py/.m carries it through into the generated C++.

WHAT THIS SCRIPT DOES

  It runs the FBA offline, once, over a grid (or a random sample) of uptake
  bounds for the exchange reactions you nominate, and writes one row per
  sample:

        v_in_0, v_in_1, ..., growth_rate

  Uptake rates are written POSITIVE (mmol/gDW/h), which is CompLaB's Fin
  convention.  Internally the bound applied to the exchange reaction is
  NEGATIVE, because an exchange reaction is written "metabolite ->" and
  consuming the metabolite therefore runs it backwards.  This script does that
  sign flip for you; do not do it twice.

USAGE

  Sweep acetate and Fe(III) uptake for Geobacter, 60 x 60 logarithmic grid:

    python3 generateTrainingData.py iAF987.xml \\
        --exchange EX_ac_e   --range 1e-3 10   --log \\
        --exchange EX_fe3_e  --range 1e-5 0.5  --log \\
        --grid 60 \\
        -o geobacter_training.csv

  Same thing with 4000 random (Latin hypercube) samples instead of a grid,
  which scales far better past two inputs:

    python3 generateTrainingData.py iAF987.xml \\
        --exchange EX_ac_e  --range 1e-3 10  --log \\
        --exchange EX_fe3_e --range 1e-5 0.5 --log \\
        --samples 4000 --seed 1 \\
        -o geobacter_training.csv

  Exchanges may be given by reaction ID ("EX_ac_e") or by the zero-based index
  used in CompLaB.xml's <exchange_reaction_indices>.  Mixing the two is fine.

CHOOSING THE RANGE

  Use the range the SIMULATION will actually visit, not the range the organism
  can tolerate.  The Michaelis-Menten bound CompLaB applies is

        v = Vmax * C / (Kc + C)

  so the largest uptake rate any voxel can ever request is Vmax, i.e. your
  <fba_maximum_uptake_flux> entry -- that is the natural upper end.  The lower
  end should be small but non-zero; sampling exactly zero adds a row that the
  network then has to fit through, and growth is usually discontinuous there.

  Sample LOGARITHMICALLY (--log) whenever the range spans more than about one
  order of magnitude.  Concentrations in a pore network routinely span four,
  and a linear grid puts essentially all of its points in the top decade, so
  the fit is worst exactly where the biology is most sensitive.

REQUIREMENTS

    python3 -m pip install cobrapy numpy

===============================================================================
"""

import argparse
import os
import sys
import time

try:
    import numpy as np
except ImportError:
    sys.exit("ERROR: numpy is required.  python3 -m pip install numpy")


# ---------------------------------------------------------------------------
# Model loading.  Identical dispatch to extractMM.py, so the same model file
# feeds the FBA path and the surrogate path with no conversion in between.
# ---------------------------------------------------------------------------
def load_model(path):
    try:
        import cobra
    except ImportError:
        sys.exit("ERROR: cobrapy is required.  python3 -m pip install cobrapy")

    ext = os.path.splitext(path)[1].lower()
    if ext == ".xml" or ext == ".sbml":
        return cobra.io.read_sbml_model(path)
    if ext == ".mat":
        return cobra.io.load_matlab_model(path)
    if ext == ".json":
        return cobra.io.load_json_model(path)
    sys.exit("ERROR: '%s' has extension '%s'.  Only .xml (SBML), .mat and "
             ".json are supported." % (path, ext))


def resolve_exchange(model, token):
    """Accept either a reaction ID or a zero-based reaction index.

    The index form is deliberate: CompLaB.xml's <exchange_reaction_indices>
    stores indices, not names, because extractMM.py strips the names when it
    flattens the model to a matrix.  Being able to paste those same numbers in
    here removes a whole class of off-by-one mistakes."""
    ids = [r.id for r in model.reactions]
    if token in ids:
        return ids.index(token)
    try:
        idx = int(token)
    except ValueError:
        sys.exit("ERROR: '%s' is neither a reaction ID in this model nor an "
                 "integer index." % token)
    if idx < 0 or idx >= len(ids):
        sys.exit("ERROR: reaction index %d is out of range; the model has %d "
                 "reactions (0..%d)." % (idx, len(ids), len(ids) - 1))
    return idx


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------
def make_grid(ranges, logflags, n):
    """Full factorial grid, n points per axis.  Exact and reproducible, but
    n**d rows -- fine for the 1-3 inputs a surrogate should have, ruinous
    beyond that.  Past three inputs use --samples."""
    axes = []
    for (lo, hi), lg in zip(ranges, logflags):
        axes.append(np.logspace(np.log10(lo), np.log10(hi), n) if lg
                    else np.linspace(lo, hi, n))
    mesh = np.meshgrid(*axes, indexing="ij")
    return np.column_stack([m.ravel() for m in mesh])


def make_lhs(ranges, logflags, n, seed):
    """Latin hypercube: n rows regardless of dimension, and every 1-D
    projection is still uniformly covered, which a random uniform sample does
    not guarantee.  No scipy dependency -- a Latin hypercube is just a
    stratified permutation per axis."""
    rng = np.random.default_rng(seed)
    d = len(ranges)
    out = np.empty((n, d))
    for j, ((lo, hi), lg) in enumerate(zip(ranges, logflags)):
        strata = (rng.permutation(n) + rng.random(n)) / n     # one point per stratum
        if lg:
            out[:, j] = 10.0 ** (np.log10(lo) + strata * (np.log10(hi) - np.log10(lo)))
        else:
            out[:, j] = lo + strata * (hi - lo)
    return out


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        prog="generateTrainingData.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Sweep a genome-scale model with FBA and write a CompLaB "
                    "surrogate training set.",
        epilog="See the header of this file for a worked example.")

    ap.add_argument("model", help="metabolic model (.xml SBML / .mat / .json)")
    ap.add_argument("--exchange", action="append", required=True, metavar="ID_OR_INDEX",
                    help="exchange reaction to sweep; repeat once per input, "
                         "IN THE SAME ORDER as the substrates in CompLaB.xml")
    ap.add_argument("--range", action="append", nargs=2, type=float, required=True,
                    metavar=("LO", "HI"),
                    help="uptake range for the matching --exchange, POSITIVE "
                         "mmol/gDW/h; repeat once per input")
    ap.add_argument("--log", action="append_const", const=True, dest="logs",
                    help="sample the matching --exchange logarithmically")
    ap.add_argument("--linear", action="append_const", const=False, dest="logs",
                    help="sample the matching --exchange linearly (the default)")
    ap.add_argument("--grid", type=int, default=0, metavar="N",
                    help="full factorial grid with N points per axis")
    ap.add_argument("--samples", type=int, default=0, metavar="N",
                    help="N Latin hypercube samples instead of a grid")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed for --samples")
    ap.add_argument("--objective", default=None, metavar="ID",
                    help="objective reaction ID (default: whatever the model "
                         "already has set)")
    ap.add_argument("--fix-others", action="store_true",
                    help="close every OTHER exchange reaction's uptake before "
                         "sweeping, so the swept substrates are the only "
                         "supply.  Off by default: most curated models come "
                         "with a sensible medium already open, and closing it "
                         "usually makes every solve infeasible.")
    ap.add_argument("-o", "--output", default="training_data.csv")
    args = ap.parse_args()

    # ---- argument cross-checks --------------------------------------------
    nin = len(args.exchange)
    if len(args.range) != nin:
        sys.exit("ERROR: %d --exchange but %d --range.  Give one --range per "
                 "--exchange, in the same order." % (nin, len(args.range)))
    logs = args.logs if args.logs else []
    if len(logs) not in (0, nin):
        sys.exit("ERROR: %d --exchange but %d --log/--linear flags.  Either "
                 "give one per exchange or none at all." % (nin, len(logs)))
    if not logs:
        logs = [False] * nin
    if (args.grid > 0) == (args.samples > 0):
        sys.exit("ERROR: give exactly one of --grid or --samples.")
    for k, (lo, hi) in enumerate(args.range):
        if lo <= 0 and logs[k]:
            sys.exit("ERROR: --range for input %d starts at %g; a logarithmic "
                     "sweep needs a strictly positive lower bound." % (k, lo))
        if hi <= lo:
            sys.exit("ERROR: --range for input %d is (%g, %g); HI must exceed "
                     "LO." % (k, lo, hi))
    if nin > 3:
        print("WARNING: %d inputs.  Surrogates fit well at 1-3 inputs; past "
              "that the sample count needed for the same accuracy grows fast "
              "and the fit is usually not worth the trouble." % nin,
              file=sys.stderr)

    # ---- model -------------------------------------------------------------
    print("Loading %s ..." % args.model)
    model = load_model(args.model)
    print("  %d metabolites, %d reactions" % (len(model.metabolites), len(model.reactions)))

    idx = [resolve_exchange(model, tok) for tok in args.exchange]
    rxn = [model.reactions[i] for i in idx]
    for k, r in enumerate(rxn):
        print("  input %d: index %d  %-16s  %s" % (k, idx[k], r.id, r.name or ""))

    if args.objective:
        model.objective = args.objective
    obj_expr = str(model.objective.expression)
    print("  objective: %s" % (obj_expr[:70] + (" ..." if len(obj_expr) > 70 else "")))

    if args.fix_others:
        swept = set(idx)
        closed = 0
        for j, r in enumerate(model.reactions):
            if j in swept:
                continue
            if r.id.startswith("EX_") and r.lower_bound < 0:
                r.lower_bound = 0.0
                closed += 1
        print("  --fix-others: closed uptake on %d other exchange reactions" % closed)

    # remember the original bounds so each sample starts from a clean state
    orig = [(r.lower_bound, r.upper_bound) for r in rxn]

    # ---- samples -----------------------------------------------------------
    if args.grid:
        X = make_grid(args.range, logs, args.grid)
        print("Grid: %d points per axis, %d samples total" % (args.grid, len(X)))
    else:
        X = make_lhs(args.range, logs, args.samples, args.seed)
        print("Latin hypercube: %d samples, seed %d" % (len(X), args.seed))

    # ---- the sweep ---------------------------------------------------------
    y = np.zeros(len(X))
    status = []
    t0 = time.time()
    for i, row in enumerate(X):
        for k, r in enumerate(rxn):
            # POSITIVE uptake in, NEGATIVE lower bound out.  Reaction.bounds
            # sets both at once; assigning lower_bound alone can raise when the
            # new lower bound momentarily exceeds the stored upper bound.
            r.bounds = (-float(row[k]), orig[k][1])
        sol = model.optimize()
        st = sol.status
        status.append(st)
        # An infeasible solve is information, not an error: it means the
        # organism cannot balance at that supply.  Zero is the physically
        # right label and keeps the sample, so the network learns where the
        # feasible region ends instead of having a hole punched in its domain.
        y[i] = float(sol.objective_value) if (st == "optimal" and sol.objective_value is not None) else 0.0
        if y[i] < 0.0:
            y[i] = 0.0
        if (i + 1) % max(1, len(X) // 20) == 0:
            el = time.time() - t0
            print("  %6d / %d   %5.1f%%   %.1f s elapsed, ~%.1f s left"
                  % (i + 1, len(X), 100.0 * (i + 1) / len(X), el,
                     el * (len(X) - i - 1) / (i + 1)))

    for k, r in enumerate(rxn):
        r.bounds = orig[k]

    nopt = sum(1 for s in status if s == "optimal")
    ngrow = int((y > 1e-12).sum())
    print("Done in %.1f s.  %d/%d solves optimal, %d/%d samples with growth."
          % (time.time() - t0, nopt, len(X), ngrow, len(X)))

    if ngrow == 0:
        print("\nWARNING: NOTHING GREW ANYWHERE.  The training set is all "
              "zeros and a network fitted to it will predict zero forever.\n"
              "  Likely causes: an essential nutrient is closed (drop "
              "--fix-others, or open it explicitly); the objective is not the "
              "biomass reaction (--objective); the swept range is below the "
              "maintenance requirement.", file=sys.stderr)
    elif ngrow < len(X) // 20:
        print("\nWARNING: only %.1f%% of samples grew.  The fit will be "
              "dominated by the flat zero region.  Consider narrowing the "
              "range to where the organism is actually alive."
              % (100.0 * ngrow / len(X)), file=sys.stderr)

    # ---- write -------------------------------------------------------------
    # The "#" lines are metadata that trainSurrogate.py/.m reads back: the
    # training box is carried all the way into the generated C++, so a run
    # that wanders outside it can be detected rather than silently trusted.
    with open(args.output, "w") as f:
        f.write("# CompLaB surrogate training set\n")
        f.write("# model: %s\n" % os.path.basename(args.model))
        f.write("# objective: %s\n" % obj_expr.replace("\n", " "))
        f.write("# n_inputs: %d\n" % nin)
        for k in range(nin):
            f.write("# input %d: %s (index %d) range %.17g %.17g %s\n"
                    % (k, rxn[k].id, idx[k], args.range[k][0], args.range[k][1],
                       "log" if logs[k] else "linear"))
        f.write("# output: specific growth rate, 1/h, range %.17g %.17g\n"
                % (y.min(), y.max()))
        f.write("# samples: %d (%d optimal, %d growing)\n" % (len(X), nopt, ngrow))
        f.write(",".join(["v_in_%d" % k for k in range(nin)]) + ",growth\n")
        for i in range(len(X)):
            f.write(",".join("%.17g" % v for v in X[i]) + ",%.17g\n" % y[i])

    print("Wrote %s  (%d rows, %d input columns)" % (args.output, len(X), nin))
    print("\nNext:  python3 trainSurrogate.py %s -o surrogate_weights.hh"
          % args.output)


if __name__ == "__main__":
    main()
