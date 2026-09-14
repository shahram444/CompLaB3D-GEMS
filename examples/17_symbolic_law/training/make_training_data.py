#!/usr/bin/env python3
"""
make_training_data.py -- write the table the symbolic search is fitted to.

The point of this case is to ask whether a search told nothing but four
operators and two column names can recover a rate law.  That question only has
an answer if we know what the law was, so the data is drawn from a law written
down here, with noise added, and the law is written into a companion file as
the answer key.

    python3 make_training_data.py                       # the shipped table
    python3 make_training_data.py --n 2000 --noise 0.05 # more, noisier
    python3 make_training_data.py --seed 7 -o other.csv

WHERE THE SAMPLING RANGE HAS TO SIT, AND WHY THIS ONE MOVED.

A half-saturation constant can only be recovered from data that reaches it.
Below it the law is a straight line in that variable, and every value of the
constant fits a straight line equally well.  The first version of this case
sampled acetate to 0.0049 mol/L against a half-saturation constant of 0.05,
which is a tenth of the way there, so the curvature that identifies the
constant was barely present: a plain product of the two concentrations already
fitted to within 4 per cent, and the search had almost nothing to find.

The concentrations were not the problem.  1e-6 to 5e-3 mol/L is what this case
actually runs at, and what a pore at a seep actually holds.  The constants were
the problem: 0.05 mol/L for acetate and 0.01 for oxygen are one to two orders
of magnitude larger than any measured aerobic heterotroph.  Bringing them to
values that are defensible on their own terms also puts them inside the
sampled range, so the same table now spans from well below each constant to
ten and twenty times above it, and the curve is fully covered.
"""
import argparse
import json
import os
import sys

import numpy as np

# The law the data comes from.  Per hour, concentrations in mol/L.
MU_MAX = 0.35          # maximum specific growth rate
KS_DONOR = 5.0e-4      # half-saturation for the electron donor, acetate
KS_ACCEPTOR = 1.0e-4   # half-saturation for the electron acceptor, oxygen

# The range the case runs over: the inlet holds 4.0e-3 acetate and 1.5e-3 o2.
LO_DONOR, HI_DONOR = 1.0e-6, 5.0e-3
LO_ACCEPTOR, HI_ACCEPTOR = 1.0e-6, 2.0e-3


def law(a, o, mu=MU_MAX, ks=KS_DONOR, ko=KS_ACCEPTOR):
    """Dual Monod: whichever of the two is scarcer holds the rate down."""
    return mu * a / (ks + a) * o / (ko + o)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", default="growth_samples.csv")
    ap.add_argument("--n", type=int, default=400, help="how many rows")
    ap.add_argument("--noise", type=float, default=0.02,
                    help="relative noise on the rate column, 0.02 = 2 per cent")
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--mu-max", type=float, default=MU_MAX)
    ap.add_argument("--ks-donor", type=float, default=KS_DONOR)
    ap.add_argument("--ks-acceptor", type=float, default=KS_ACCEPTOR)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    # Log-uniform, because a rate law is used across orders of magnitude and a
    # uniform draw would put almost every sample in the saturated corner where
    # the law is flat and carries no information about the constants.
    a = np.exp(rng.uniform(np.log(LO_DONOR), np.log(HI_DONOR), args.n))
    o = np.exp(rng.uniform(np.log(LO_ACCEPTOR), np.log(HI_ACCEPTOR), args.n))
    g = law(a, o, args.mu_max, args.ks_donor, args.ks_acceptor)
    g = g * (1.0 + args.noise * rng.standard_normal(args.n))
    g = np.maximum(g, 0.0)

    with open(args.output, "w", newline="") as f:
        f.write("acetate,o2,growth\n")
        for i in range(args.n):
            f.write("%.7e,%.7e,%.7e\n" % (a[i], o[i], g[i]))

    truth = {
        "law": "%g * acetate / (%g + acetate) * o2 / (%g + o2)"
               % (args.mu_max, args.ks_donor, args.ks_acceptor),
        "units": "per_hour",
        "mu_max": args.mu_max,
        "ks_acetate": args.ks_donor,
        "ks_o2": args.ks_acceptor,
        "samples": args.n,
        "relative_noise": args.noise,
        "seed": args.seed,
        "acetate_range": [LO_DONOR, HI_DONOR],
        "o2_range": [LO_ACCEPTOR, HI_ACCEPTOR],
    }
    side = os.path.splitext(args.output)[0] + "_truth.json"
    with open(side, "w") as f:
        json.dump(truth, f, indent=2)
        f.write("\n")

    print("wrote %s   %d rows, %.0f%% noise, seed %d"
          % (args.output, args.n, 100 * args.noise, args.seed))
    print("the law behind it, written to %s:" % os.path.basename(side))
    print("    growth = %s" % truth["law"])
    print("how far the sampling reaches past each half-saturation constant,")
    print("which is what decides whether the constants can be recovered at all:")
    print("    acetate  %.4g .. %.4g   =  %.3f to %.1f times its constant"
          % (a.min(), a.max(), a.min() / args.ks_donor, a.max() / args.ks_donor))
    print("    o2       %.4g .. %.4g   =  %.3f to %.1f times its constant"
          % (o.min(), o.max(), o.min() / args.ks_acceptor, o.max() / args.ks_acceptor))
    if a.max() / args.ks_donor < 3 or o.max() / args.ks_acceptor < 3:
        print("\n  WARNING: the sampling stops short of saturation. Below its")
        print("  half-saturation constant the law is a straight line in that")
        print("  variable, so the constant cannot be identified from this table")
        print("  and the search will return a product rather than a Monod form.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
