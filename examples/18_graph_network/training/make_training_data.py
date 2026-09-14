#!/usr/bin/env python3
"""
make_training_data.py -- write the table the graph network is fitted to.

The shipped aom_samples.csv was drawn from the law written down here, so this
case has an answer key: tests/check_aom.py holds the trained network up against
the same law on 300 samples it never saw.  This script is what makes such a
table, and it records the law beside it so the key cannot drift from the data.

    python3 make_training_data.py                        # an equivalent table
    python3 make_training_data.py --n 2000 --noise 0.02  # more, with noise
    python3 make_training_data.py -o wider.csv --ch4-hi 2e-2

WHAT THE SHIPPED TABLE IS, EXACTLY.

Anaerobic oxidation of methane coupled to sulfate reduction:

    CH4 + SO4(2-)  ->  HS(-) + HCO3(-)

One turn of the reaction consumes one methane and one sulfate and makes one
sulfide and one bicarbonate, which is the aom_stoich.csv beside this file.  The
rate of that turn is dual Monod on the two reactants, and every species rate is
its stoichiometric number times that rate.  Growth is the rate times the yield.

The shipped table carries NO noise.  That is deliberate and worth knowing when
you read the accuracy figures: whatever error the network shows is its own
approximation error, not scatter it could not have fitted.  Pass --noise to add
some if you want to see how it copes.

REGENERATING MEANS RETRAINING.  aom_samples.csv and the shipped
pipelines/B_offline_models/B4_graph_network/expected/aom.gnn were made as a
pair.  If you overwrite the table, run ./offline.sh to fit a new network to it,
or the two stop belonging together.  This script therefore writes to
aom_samples.csv only if you ask for that name.
"""
import argparse
import json
import os
import sys

import numpy as np

# The law the data comes from.  Rate in mol/L/h, concentrations in mol/L.
VMAX = 4.0e-3
K_CH4 = 1.0e-3
K_SO4 = 5.0e-4
YIELD = 5.0            # biomass made per unit of methane turned over

# One turn of CH4 + SO4 -> HS + HCO3.  Must agree with aom_stoich.csv.
STOICH = [("CH4", -1.0), ("SO4", -1.0), ("HS", +1.0), ("HCO3", +1.0)]

# The range each species is sampled over.  CH4 reaches 5 times its constant and
# SO4 16 times its, so the saturating part of the curve is well covered and both
# constants can be recovered from the table.  HS and HCO3 are products: they do
# not enter the law, and they are sampled so the network has to learn that.
RANGES = {"CH4": (1.0e-5, 5.0e-3), "SO4": (1.0e-5, 8.0e-3),
          "HS": (0.0, 2.0e-3), "HCO3": (1.0e-4, 4.0e-3)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", default="aom_samples_new.csv")
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--noise", type=float, default=0.0,
                    help="relative noise on the rate columns; the shipped table has none")
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--vmax", type=float, default=VMAX)
    ap.add_argument("--k-ch4", type=float, default=K_CH4)
    ap.add_argument("--k-so4", type=float, default=K_SO4)
    ap.add_argument("--yield-coef", type=float, default=YIELD)
    ap.add_argument("--ch4-hi", type=float, default=RANGES["CH4"][1])
    ap.add_argument("--so4-hi", type=float, default=RANGES["SO4"][1])
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    hi = dict(RANGES)
    hi["CH4"] = (RANGES["CH4"][0], args.ch4_hi)
    hi["SO4"] = (RANGES["SO4"][0], args.so4_hi)
    C = {k: rng.uniform(lo, up, args.n) for k, (lo, up) in hi.items()}

    extent = (args.vmax * C["CH4"] / (args.k_ch4 + C["CH4"])
              * C["SO4"] / (args.k_so4 + C["SO4"]))
    if args.noise > 0:
        extent = np.maximum(extent * (1.0 + args.noise * rng.standard_normal(args.n)), 0.0)

    cols = ["CH4", "SO4", "HS", "HCO3"]
    with open(args.output, "w", newline="") as f:
        f.write(",".join(cols + [c + "_rate" for c, _ in STOICH] + ["growth"]) + "\n")
        for i in range(args.n):
            row = ["%.7e" % C[c][i] for c in cols]
            row += ["%.7e" % (nu * extent[i]) for _, nu in STOICH]
            row += ["%.7e" % (args.yield_coef * extent[i])]
            f.write(",".join(row) + "\n")

    truth = {
        "reaction": "CH4 + SO4 -> HS + HCO3",
        "extent_law": "%g * CH4 / (%g + CH4) * SO4 / (%g + SO4)"
                      % (args.vmax, args.k_ch4, args.k_so4),
        "units": "mol/L/h",
        "vmax": args.vmax, "k_ch4": args.k_ch4, "k_so4": args.k_so4,
        "yield": args.yield_coef, "stoichiometry": dict(STOICH),
        "samples": args.n, "relative_noise": args.noise, "seed": args.seed,
        "ranges": {k: list(v) for k, v in hi.items()},
    }
    side = os.path.splitext(args.output)[0] + "_truth.json"
    with open(side, "w") as f:
        json.dump(truth, f, indent=2)
        f.write("\n")

    print("wrote %s   %d rows, %.0f%% noise, seed %d"
          % (args.output, args.n, 100 * args.noise, args.seed))
    print("the law behind it, written to %s:" % os.path.basename(side))
    print("    extent = %s" % truth["extent_law"])
    print("    every species rate is its stoichiometric number times that")
    print("how far the sampling reaches past each half-saturation constant:")
    print("    CH4  to %.1f times its constant" % (hi["CH4"][1] / args.k_ch4))
    print("    SO4  to %.1f times its constant" % (hi["SO4"][1] / args.k_so4))
    if os.path.basename(args.output) == "aom_samples.csv":
        print("\n  You have overwritten the shipped table. Run ./offline.sh to fit a")
        print("  new network to it, or aom.gnn no longer belongs with this data.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
