#!/usr/bin/env python3
"""
20_upscaling  --  POST-PROCESSING

Reads output/upscaling.csv and says whether the effectiveness factor this run
reports is worth quoting, and what it means next to the classical result.

  reads   output/upscaling.csv, and output/summary.csv if it is there
  writes  nothing

Standard library only.

    python3 postprocess.py

THREE QUESTIONS, IN ORDER. Each one can invalidate the answer to the next, so
they are asked in this order and the script says which one failed.

  1. Did the measurement start from a state it could measure?
     At iteration zero every voxel is still at the bulk composition, so the
     measured mean rate must equal the rate at the bulk and eta must be 1. If it
     is not, the chain from the increment lattices to the reported number is
     broken and nothing below it means anything.

  2. Did it reach steady state?
     eta is a property of the concentration profile inside the aggregate, which
     relaxes on R^2/D. Read off a run that has not got there, it is a number that
     is still moving.

  3. What did the energy limit do to it?
     The classical sphere result assumes first-order kinetics and no energy
     limit. This case has neither. The gap is the whole point.
"""
import csv
import math
import os
import sys

UPS = os.path.join("output", "upscaling.csv")


def rows(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return [{k: float(v) for k, v in r.items()}
                for r in csv.DictReader(l for l in f if not l.startswith("#"))]


def rule(title=""):
    print("-" * 78)
    if title:
        print(title)
        print("-" * 78)


def main():
    print("=" * 78)
    print("20_upscaling")
    print("=" * 78)

    r = rows(UPS)
    if not r:
        print("\nNo %s." % UPS)
        print("The <upscaling> block writes it. Check that <enabled> is true and that")
        print("the run got past its first diagnostic interval.")
        return 1

    first, last = r[0], r[-1]
    print("\n%d record(s), iteration %d to %d"
          % (len(r), int(first["iteration"]), int(last["iteration"])))
    print("aggregate %d voxel(s), bulk %d voxel(s)"
          % (int(last["aggregate_voxels"]), int(last["bulk_voxels"])))

    # ---- 1. the measurement checking itself ------------------------------
    rule("1. DOES THE MEASUREMENT CHAIN WORK AT ALL")
    e0 = first["eta"]
    print("   At iteration 0 every voxel is still at the bulk composition, so the")
    print("   mean rate inside the aggregate must equal the rate at the bulk, and")
    print("   eta must be exactly 1. This is the measurement checking itself.")
    print("      eta(0) = %.12g   (departure from 1: %.3g)" % (e0, abs(e0 - 1.0)))
    ok1 = abs(e0 - 1.0) < 1e-6
    print("      %s" % ("PASS" if ok1 else "FAIL -- nothing below this is meaningful"))
    if not ok1:
        print("\n   A non-unit value here means the numerator and the denominator are not")
        print("   the same quantity. The numerator is read from the dC increment")
        print("   lattices, which hold the increment as computeDensity() -- sum(f)+1,")
        print("   not sum(f). The denominator is defineRxnKinetics evaluated at the")
        print("   bulk. Check those two before anything else.")
        return 1

    # ---- 2. steady state -------------------------------------------------
    rule("2. HAS IT STOPPED MOVING")
    if len(r) >= 3:
        tail = [x["eta"] for x in r[-3:]]
        spread = (max(tail) - min(tail)) / abs(max(tail)) if max(tail) else 0.0
        print("   the last three intervals: " + ", ".join("%.6g" % t for t in tail))
        print("   spread %.3g%%" % (100.0 * spread))
        steady = spread < 0.01
        print("      %s" % ("PASS -- steady" if steady
                            else "FAIL -- still a transient. Raise <ade_max_iT>."))
    else:
        steady = False
        print("   fewer than three intervals; raise <ade_max_iT> or lower <interval>")

    # ---- 3. what the energy limit did ------------------------------------
    rule("3. THE TWO NUMBERS A CONTINUUM MODEL NEEDS")
    print("   bulk concentration        %.6g mol/L" % last["bulk_conc"])
    print("   rate at the bulk          %.6g mol/L/s" % last["rate_at_bulk"])
    print("   measured mean rate        %.6g mol/L/s" % last["rate_mean"])
    print()
    print("   effectiveness factor eta  %.6g" % last["eta"])
    print("   Thiele modulus phi        %.6g" % last["thiele"])
    print("   classical eta at that phi %.6g" % last["eta_classical"])
    print("   gate at the bulk F_T      %.6g" % last["F_T_bulk"])

    if last["eta_classical"] > 0:
        ratio = last["eta"] / last["eta_classical"]
        print()
        print("   measured / classical      %.4g" % ratio)
        print()
        if ratio > 1.05:
            print("   ABOVE the classical curve. The classical result is first order in the")
            print("   substrate; a Monod law saturates, so the interior keeps reacting at")
            print("   nearly the full rate on a smaller concentration than first order")
            print("   would allow. An aggregate can be MORE effective than the textbook")
            print("   curve says, and this is when.")
        elif ratio < 0.95:
            print("   BELOW the classical curve. This is the energy limit. The core of the")
            print("   aggregate has not merely slowed down -- it has stopped, at the depth")
            print("   where its own sulfide raised dG past what the organism can use. That")
            print("   is a moving internal boundary, and no first-order Thiele analysis")
            print("   contains it.")
        else:
            print("   The two agree here. That happens where the aggregate is small enough,")
            print("   or the bulk far enough above the threshold, that neither the Monod")
            print("   saturation nor the gate has taken hold.")

    # ---- what a reader should do next ------------------------------------
    rule("WHAT TO DO WITH THIS")
    print("Evaluating F_T at the bulk is NOT a shortcut past any of the above.")
    print("r(C_bulk) already includes F_T(C_bulk) -- it is in the denominator -- so")
    print("whatever eta departs from 1 is exactly the error that shortcut makes.")
    print("Here that error is %.0f%%." % (100.0 * abs(1.0 - last["eta"])))
    print()
    print("ONE RUN IS ONE POINT. A continuum model needs the curve, over the sizes")
    print("and compositions it will actually be asked about:")
    print()
    print("    python3 offline/upscale.py --radii 3,4,5 --bulk 1e-5,3e-4,1e-3")
    print()
    if not steady:
        print("Do not quote the numbers above until step 2 passes.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
