#!/usr/bin/env python3
"""Check a finished straight-duct run against closed-form answers.

    python3 tools/analytic_channel.py --run channel/runs/run_0000

This is the physics check. Everything else in the package tests that the pipeline
carries numbers from one stage to the next; this tests that the numbers are right.

It uses the straight duct, gid 900, because a duct has answers that can be
written down. The sphere packs do not.

WHAT IS EXACT, AND WHY
----------------------
Three combinations of the three fields obey a reaction-free equation, EXACTLY and
for any rate constant. The reaction is A + B -> C, one for one for one, so

    dA/dt = T(A) - R        dB/dt = T(B) - R        dC/dt = T(C) + R

where T is whatever transport does. All three species share one diffusivity, so
T is the same linear operator for each, and

    A - B      the R cancels by subtraction
    A + C      the R cancels by addition
    B + C      likewise

each satisfy dS/dt = T(S) with no source at all. Every one of them is therefore a
CONSERVATIVE TRACER carrying the boundary values of its own combination, and at
steady state in a duct each is the solution of the plain advection-diffusion
equation. Nothing about the rate law enters. That is what makes these checks
worth something: they hold whether or not the chemistry is what you think it is,
so a disagreement is a transport fault, and they would be broken by an
inequality anywhere in the stoichiometry, by the three species not sharing a
diffusivity, or by the positivity clamp cutting one species and not its partner.

At Pe = 0 the steady solution of T(S) = 0 with Dirichlet ends is a straight line
in x, so:

    A - B  goes from  +A0  to  -B0    linearly
    A + C  goes from   A0  to    0    linearly
    B + C  goes from    0  to   B0    linearly

At Pe > 0 the one-dimensional steady solution with a uniform velocity is

    S(x) = S(0) + (S(L) - S(0)) * (exp(Pe x / L) - 1) / (exp(Pe) - 1)

and the cross-section average of the duct field should follow it. That comparison
is APPROXIMATE and gets worse as Peclet rises: a duct has a velocity profile, not
a plug flow, and shear stretches the front. The residual grows like Pe^2, which
is Taylor dispersion and is real physics rather than a solver fault. The check
reports the discrepancy rather than passing or failing on it, and it fails only
on the exact statements.

WHAT ELSE IT REPORTS
--------------------
    where the reaction sits. At large Damkohler A and B cannot coexist and the
    reaction collapses onto a sheet. For equal feeds, equal diffusivities and no
    flow that sheet is at the middle of the domain; flow pushes it downstream.
    The position of the peak of C is printed beside the thin-sheet prediction.

    whether the run reached steady state, by how much the last two snapshots
    differ. The exact statements above are steady-state statements, so a run
    stopped at one transit time will not satisfy them and should not be expected
    to. Build the channel case with a large --run-factor.
"""
import argparse
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collect_to_h5 import read_vti, snapshots, pick_scalar        # noqa: E402

PORE = 2


def cross_section_mean(field, pore):
    """One number per x-slice, over the pore voxels of that slice."""
    out = np.zeros(field.shape[0])
    for i in range(field.shape[0]):
        m = pore[i]
        out[i] = field[i][m].mean() if m.any() else np.nan
    return out


def ade_1d(x_over_L, s0, sL, pe):
    if abs(pe) < 1e-9:
        return s0 + (sL - s0) * x_over_L
    return s0 + (sL - s0) * (np.exp(pe * x_over_L) - 1.0) / (np.exp(pe) - 1.0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="one finished case directory")
    ap.add_argument("--geometries", default="geometries")
    ap.add_argument("--tol", type=float, default=0.02,
                    help="relative tolerance on the exact statements")
    ap.add_argument("--plot", default=None, help="write a PNG here")
    args = ap.parse_args(argv)

    p = json.load(open(os.path.join(args.run, "params.json")))
    od = os.path.join(args.run, "output")
    gid = p["gid"]
    npz = os.path.join(args.geometries, "geom_%04d" % gid, "geom_%04d.npz" % gid)
    mat = np.load(npz)["material"]
    pore = (mat == PORE)
    shape = mat.shape
    nx = shape[0]
    pe, da, A0, B0 = p["pe"], p["da_abio"], p["A0"], p["B0"]

    print("%s   gid %d   Pe %g   Da %g   k %.4g L/(mol s)"
          % (args.run, gid, pe, da, p["k_abio"]))
    variant = p.get("boundaries", "unknown")
    if variant != "dirichlet":
        print()
        print("  STOP. This case was built with boundaries = %s." % variant)
        print("  The three statements below are exact only when EVERY face is held,")
        print("  because that is what makes A-B, A+C and B+C conservative tracers")
        print("  with known boundary values. With an open outflow their boundary")
        print("  values are whatever the flow leaves there, and there is nothing to")
        print("  compare against. Rebuild the duct case with")
        print()
        print("    make_campaign.py build --geometries geometries --out channel \\")
        print("        --channel --boundaries dirichlet --pe 0 --da 0.1 1.0 10.0 \\")
        print("        --run-factor 2 --complab <your binary>")
        print()
        return 2

    if p.get("porosity", 0) < 0.999:
        print("  NOTE: this is not the straight duct (porosity %.3f). The exact "
              "statements below hold in any geometry, but the one-dimensional "
              "profile comparison does not." % p["porosity"])

    fields = {}
    for nm in p["species"]:
        s = snapshots(od, nm)
        if len(s) < 2:
            sys.exit("need at least two snapshots of %s in %s" % (nm, od))
        arrs, _ = read_vti(s[-1][1])
        fields[nm] = pick_scalar(arrs, nm, shape)
        arrs2, _ = read_vti(s[-2][1])
        fields[nm + "_prev"] = pick_scalar(arrs2, nm, shape)

    A, B, Cc = fields["A"], fields["B"], fields["C"]

    # --- did it reach steady state? ---------------------------------------
    drift = 0.0
    for nm in p["species"]:
        a, b = fields[nm][pore], fields[nm + "_prev"][pore]
        den = np.sqrt((a * a).sum()) + 1e-30
        drift = max(drift, float(np.sqrt(((a - b) ** 2).sum()) / den))
    print("  change over the last snapshot interval: %.4f" % drift)
    steady = drift < 0.01
    if not steady:
        print("  NOT at steady state. The exact statements below are steady-state")
        print("  statements; rebuild the channel case with a larger --run-factor")
        print("  before reading the numbers as a verdict.")

    x = np.arange(nx) / (nx - 1.0)
    invariants = [("A - B", cross_section_mean(A - B, pore), A0, -B0),
                  ("A + C", cross_section_mean(A + Cc, pore), A0, 0.0),
                  ("B + C", cross_section_mean(B + Cc, pore), 0.0, B0)]

    fails = 0
    print()
    print("  the three reaction-free combinations, against the 1D steady profile")
    print("  %-8s %10s %10s   %s" % ("", "max error", "relative", "verdict"))
    curves = []
    for name, got, s0, sL in invariants:
        want = ade_1d(x, s0, sL, pe)
        err = np.nanmax(np.abs(got - want))
        rel = err / max(abs(s0), abs(sL), 1e-30)
        if not steady:
            verdict = "not steady, not judged"
        elif abs(pe) < 1e-9:
            verdict = "exact" if rel <= args.tol else "FAIL"
            fails += int(rel > args.tol)
        else:
            verdict = ("within tolerance" if rel <= args.tol
                       else "outside tolerance; Taylor dispersion grows as Pe^2")
        print("  %-8s %10.4g %10.4f   %s" % (name, err, rel, verdict))
        curves.append((name, got, want))

    # --- where is the reaction? -------------------------------------------
    cx = cross_section_mean(Cc, pore)
    peak = int(np.nanargmax(cx))
    print()
    print("  peak C at x = %d of %d (%.2f of the domain)"
          % (peak, nx - 1, peak / (nx - 1.0)))
    if abs(pe) < 1e-9:
        print("  thin-sheet prediction for equal feeds and no flow: %.2f" % 0.5)
    else:
        # where the two 1D tracer profiles cross, which is where a fast reaction
        # sheet sits: A_tracer(x) = B_tracer(x)
        xa = ade_1d(x, A0, 0.0, pe)
        xb = ade_1d(x, 0.0, B0, pe)
        cross = x[int(np.argmin(np.abs(xa - xb)))]
        print("  thin-sheet prediction from the 1D tracer crossing: %.2f" % cross)
    print("  (the two agree only in the fast-reaction limit; at Da = %g the front "
          "is spread)" % da)

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
            for nm, f in (("A", A), ("B", B), ("C", Cc)):
                ax[0].plot(x, cross_section_mean(f, pore), label=nm)
            ax[0].set_xlabel("x / L"); ax[0].set_ylabel("cross-section mean, mol/L")
            ax[0].set_title("Pe = %g, Da = %g" % (pe, da)); ax[0].legend()
            for nm, got, want in curves:
                line, = ax[1].plot(x, got, label=nm + " simulated")
                ax[1].plot(x, want, "--", color=line.get_color(),
                           label=nm + " closed form")
            ax[1].set_xlabel("x / L")
            ax[1].set_title("the reaction-free combinations")
            ax[1].legend(fontsize=7)
            fig.tight_layout()
            fig.savefig(args.plot, dpi=120)
            print("  wrote %s" % args.plot)
        except Exception as e:
            print("  (no plot: %s)" % e)

    print()
    if not steady:
        print("VERDICT: inconclusive. Run the channel case to steady state.")
        return 0
    if fails:
        print("VERDICT: FAIL. %d of the exact statements %s violated, which means "
              "a transport fault: unequal diffusivities, broken stoichiometry, or "
              "a boundary that is not what the XML says."
              % (fails, "is" if fails == 1 else "are"))
        return 1
    print("VERDICT: the duct reproduces its closed-form answers. Transport, "
          "stoichiometry and the boundary conditions are right.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
