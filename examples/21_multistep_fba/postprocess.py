#!/usr/bin/env python3
"""
21_multistep_fba  --  POST-PROCESSING

Growth from a linear program solved in STAGES, with a fraction of each stage's answer retained for the next.

  reads   output/summary.csv   one row per diagnostic interval
          output/*.log         the solver log, if it was captured
  writes  nothing; prints a verdict

Self-contained: standard library only, nothing imported from tools/. Run it
after the solver, or let pipeline.sh run it for you:

    python3 postprocess.py

The point of this case is the <multi_step> block. Instead of one linear program per
voxel, an early stage is solved, <retain_fraction> of its answer is kept, and a
later stage is solved against what remains. Acetate and formate are products, so
they start at zero and can only be made.

WHICH TOTALS CAN BE BALANCED AT ALL. A species held at a Dirichlet boundary is
supplied from outside the domain, so its total is not conserved and never will
be. Naming it in a conservation check guarantees a failure that means nothing.
In this case `glc` and `o2` are fed and `ac` and `for` are closed.
"""
import csv
import glob
import os
import sys

OUT = "output"


def rows():
    p = os.path.join(OUT, "summary.csv")
    if not os.path.exists(p):
        sys.exit("no %s. Run the solver first." % p)
    with open(p) as f:
        return [r for r in csv.DictReader(l for l in f if not l.startswith("#"))]


def series(rs, key):
    return [float(r[key]) for r in rs if key in r and r[key] not in ("", None)]


def main():
    rs = rows()
    if len(rs) < 2:
        sys.exit("summary.csv has %d row(s); at least two are needed to say anything "
                 "changed. Lower <interval> or raise <ade_max_iT>." % len(rs))
    print("Run      : %s" % os.path.abspath("."))
    print("Snapshots: %d" % len(rs))
    print()

    bad = 0
    for name in ['glc', 'o2', 'ac', 'for']:
        k = name + "_total"
        s = series(rs, k)
        if not s:
            print("  %-10s no column in summary.csv" % name)
            continue
        d = s[-1] - s[0]
        word = "rose" if d > 0 else ("fell" if d < 0 else "unchanged")
        print("  %-10s %.6e -> %.6e   %s" % (name, s[0], s[-1], word))

    for name in ['glc', 'o2', 'ac', 'for']:
        k = name + "_min"
        s = series(rs, k)
        if s and min(s) < 0:
            print("\n  NEGATIVE VALUES in: %s (worst %.3e)" % (name, min(s)))
            print("  [v1.3] The compiled kinetics and flux balance paths apply each")
            print("  increment as computed. A time step or a rate constant large enough to")
            print("  consume more than a voxel holds will drive it negative. Reduce")
            print("  <ade_dt> or the uptake bounds, then run again.")
            bad += 1

    for f in sorted(glob.glob(os.path.join(OUT, "*.log"))):
        txt = open(f, errors="replace").read()
        for marker in ("infeasible", "did not converge", "NEG!"):
            n = txt.lower().count(marker.lower())
            if n:
                print("\n  log reports %r %d time(s) in %s" % (marker, n, os.path.basename(f)))
                bad += 1

    print("\nWhat to check")
    print("  1. acetate and formate appear")
    print("     Both start at zero and are closed at both ends, so any amount present was made by the organism. If they stay at zero the later stage never ran.")
    print("  2. glucose and oxygen fall")
    print("     Both are fed from the left, so their totals need not fall; what must fall is the concentration away from the inlet.")
    print("  3. no field goes negative")
    print("     A staged solve draws twice from the same voxel, so this is the check that the retained fraction was subtracted from what the second stage was allowed.")
    print("  4. growth is not zero everywhere")
    print("     Zero everywhere usually means the model is closed: an exchange the configuration does not open.")

    print("\nVERDICT: %s" % ("LOOK AGAIN" if bad else "nothing obviously wrong"))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
