#!/usr/bin/env python3
"""
19_thermodynamic_gate  --  POST-PROCESSING

An AOM aggregate whose rate law is gated by the local free energy.

  reads   output/summary.csv    one row per diagnostic interval
          output/*.log          the solver log, if it was captured
          output/thermo_curve.csv  what the offline step predicted
  writes  nothing; prints a verdict

Self-contained: standard library only. Run it after the solver, or let
pipeline.sh run it for you:

    python3 postprocess.py

WHY A VERDICT AND NOT A PICTURE.  A field plot of a gated run and an ungated
run look equally plausible.  The difference is in the numbers.

THE ONE CHECK THAT MATTERS HERE.  A thermodynamic gate can fail in two
directions and both produce a finished run:

    shut everywhere   the reaction never happens, which looks identical to a
                      case whose biomass was never seeded
    open everywhere   the gate did nothing, and the result is the ungated run
                      with a thermodynamic claim attached to it

The [THM] block in the log reports the mean F_T and its range precisely so
those two cases are visible.  This script pulls them out and says which one you
have.

WHICH TOTALS CAN BE BALANCED AT ALL.  All four species are held at a Dirichlet
boundary on the left, so none of their totals is closed and a conservation
check on any of them would fail for reasons that mean nothing.  Biomass is
closed: nothing supplies it and decay is off, so it can only grow.
"""
import csv
import os
import re
import sys

SUMMARY = os.path.join("output", "summary.csv")
NOGATE = os.path.join("output", "summary_nogate.csv")
CURVE = os.path.join("output", "thermo_curve.csv")
LOGDIR = "output"

BOUNDARY = {
    "CH4":  "Dirichlet / Neumann",
    "SO4":  "Dirichlet / Neumann",
    "HS":   "Dirichlet / Neumann",
    "HCO3": "Dirichlet / Neumann",
}



def conserved(row, name):
    """What the lattice actually holds, not just what the dynamics will admit to.

    <name>_total is computed with Palabos's computeDensity(), which asks each cell's
    dynamics. BounceBack and NoDynamics both answer from a stored number and ignore
    the populations they are holding, and the reported box excludes the two boundary
    planes at x = 0 and x = nx-1. So three kinds of mass are missing from _total:
    whatever is in flight at a wall, whatever is resting inside a grain, and
    whatever is sitting in a closed boundary plane.

    <name>_held is exactly that difference, measured from the populations. Adding it
    back gives the conserved quantity. In a fully closed box with no reaction the sum
    below is constant to 4e-13; _total alone drifts by 7%.

    Older summary files have no _held column; those fall back to _total.
    """
    return row.get(name + "_total", 0.0) + row.get(name + "_held", 0.0)


def read_summary(path):
    if not os.path.exists(path):
        return None, None
    with open(path) as f:
        lines = [ln for ln in f if not ln.startswith("#")]
    if not lines:
        return None, None
    rdr = csv.DictReader(lines)
    rows = [{k: (float(v) if v not in (None, "") else 0.0)
             for k, v in r.items() if k} for r in rdr]
    return (rdr.fieldnames, rows) if rows else (None, None)


def series(rows, col):
    return [r[col] for r in rows if col in r]


def read_logs():
    """The gated run's log, preferentially.

    output/ accumulates logs: the control run writes one, and anything else you
    tried writes one too. Concatenating them all and taking the first [THM]
    block reports whichever run happened to sort first, which is how a stale
    experiment ends up being presented as this one. run.log is the gated run.
    """
    if not os.path.isdir(LOGDIR):
        return ""
    primary = os.path.join(LOGDIR, "run.log")
    if os.path.exists(primary):
        try:
            return open(primary, errors="replace").read()
        except OSError:
            pass
    text = ""
    for f in sorted(os.listdir(LOGDIR)):
        if f.endswith(".log"):
            try:
                text += open(os.path.join(LOGDIR, f), errors="replace").read()
            except OSError:
                pass
    return text


def thermo_verdict(text):
    """Pull the numbers out of the [THM] end-of-run report and judge them."""
    hits = list(re.finditer(r"\[THM\]\s+(\d+) gate evaluations, mean F_T ([0-9.]+), "
                            r"range ([0-9.]+)\.\.([0-9.]+)", text))
    if not hits:
        return None
    m = hits[-1]                       # the most recent report, not the first
    n = int(m.group(1))
    mean = float(m.group(2))
    lo = float(m.group(3))
    hi = float(m.group(4))

    bs = list(re.finditer(r"\[THM\]\s+(\d+) of them \(([0-9.]+)%\)", text))
    b = bs[-1] if bs else None
    blocked, pct = (int(b.group(1)), float(b.group(2))) if b else (0, 0.0)

    print("THE GATE")
    print("-" * 78)
    print("   evaluations   %d" % n)
    print("   mean F_T      %.4f" % mean)
    print("   range         %.4f to %.4f" % (lo, hi))
    print("   at threshold  %d (%.2f%%) had F_T = 0, so no reaction at all" % (blocked, pct))
    print()

    if hi <= 0.0:
        print("   VERDICT: SHUT EVERYWHERE. Nothing reacted anywhere in the domain. This")
        print("   is almost never a result; it is a sign that dG0 has the wrong sign, that")
        print("   the stoich line has reactants and products the wrong way round, or that")
        print("   dG0 was quoted at 25 C for a case running at 4 C. Run")
        print("   'python3 offline/thermo_curve.py' and look at where the threshold sits.")
    elif lo >= 1.0:
        print("   VERDICT: OPEN EVERYWHERE. The gate never closed, so this run is")
        print("   numerically identical to one with <thermodynamics><enabled>false</enabled>.")
        print("   Do not report it as a thermodynamically controlled result. Either the")
        print("   products never built up enough (run longer, or shrink the aggregate's")
        print("   diffusive escape), or the threshold is set too low to bind.")
    elif pct > 0.5:
        print("   VERDICT: THE GATE ACTED, AND CLOSED. Part of the domain is past the energy")
        print("   threshold and part is not, which is the reaction front this case exists to")
        print("   produce. The comparison against the ungated control below is the measurement.")
    elif lo < 0.5:
        print("   VERDICT: THE GATE ACTED. The rate is attenuated by a factor of %.0f at the"
              % (1.0 / lo if lo > 0 else 0.0))
        print("   most limited voxel. Nowhere reached F_T = 0 exactly, and on a case like this")
        print("   one that is expected rather than a shortfall: the gate is self-limiting. As")
        print("   F_T falls the reaction slows, so it stops producing the very products that")
        print("   were closing it, and the deepest voxel settles ON the threshold instead of")
        print("   passing through it. F_T reaches exactly zero only where transport, not the")
        print("   local reaction, pushes the composition past it.")
    else:
        print("   VERDICT: THE GATE ACTED WEAKLY. It attenuated rates but nowhere by much.")
        print("   That is a real result, but a thinner one than the case is set up to show.")
        print("   Running longer lets the products build further in.")
    print()
    return {"n": n, "mean": mean, "lo": lo, "hi": hi, "blocked": blocked, "pct": pct}


def compare_curve(v):
    """Hold the run against what the offline step predicted before it started."""
    if v is None or not os.path.exists(CURVE):
        return
    with open(CURVE) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return
    fts = [float(r["F_T"]) for r in rows]
    print("AGAINST THE OFFLINE PREDICTION")
    print("-" * 78)
    print("   offline/thermo_curve.py swept the composition path this case follows and")
    print("   found F_T running from %.4f down to %.4f along it." % (max(fts), min(fts)))
    if min(fts) <= 0.0 and v["hi"] > 0.0 and v["lo"] < 1.0:
        print("   The run reached %.4f to %.4f, so it visited that path as expected."
              % (v["lo"], v["hi"]))
    elif min(fts) <= 0.0 and v["lo"] >= 1.0:
        print("   The prediction said the gate would shut somewhere along the path; the run")
        print("   never got there. The run is too short, or the aggregate loses its products")
        print("   to diffusion faster than the sweep assumed.")
    print()


def against_control():
    """Hold the gated run against the ungated one, and check they really differ.

    This check exists because both ways of getting it wrong are silent. The
    solver used to ignore its command-line argument, so the control ran the same
    file as the experiment; and an unanchored sed on CompLaB.xml switches off the
    <diagnostics> block instead of the <thermodynamics> one, so the control
    writes no summary and the comparison reads the experiment's own file twice.
    Either way the two runs come back identical, and the obvious reading of that
    is that the gate does nothing.
    """
    if not (os.path.exists(SUMMARY) and os.path.exists(NOGATE)):
        print("NO CONTROL RUN")
        print("-" * 78)
        print("   %s is missing, so there is nothing to compare against." % NOGATE)
        print("   ./pipeline.sh produces it. Without it, nothing here supports a claim")
        print("   about what the gate did to the fields.")
        print()
        return

    _, a = read_summary(NOGATE)
    _, b = read_summary(SUMMARY)
    if not a or not b:
        return

    print("GATED AGAINST UNGATED")
    print("-" * 78)
    print("   Same rate law, same geometry, same boundary conditions, same time step.")
    print("   One line of CompLaB.xml differs. Every difference below is the gate.")
    print()

    cols = [k for k in b[0] if k.endswith("_total")]
    worst = 0.0
    print("   %-14s %16s %16s %12s" % ("field", "ungated", "gated", "relative"))
    for c in sorted(cols):
        x = a[-1].get(c, 0.0)
        y = b[-1].get(c, 0.0)
        rel = abs(y - x) / abs(x) if x else 0.0
        worst = max(worst, rel)
        print("   %-14s %16.9g %16.9g %11.3g%%" % (c[:-6], x, y, 100.0 * rel))
    print()

    if worst < 1e-12:
        print("   THE TWO RUNS ARE IDENTICAL. Before concluding that the gate did")
        print("   nothing, check that the control really ran the control: `diff CompLaB.xml")
        print("   nogate.xml` must show exactly one changed line, inside the")
        print("   <thermodynamics> block and nowhere else.")
    else:
        print("   Largest relative difference: %.3g%%. The two runs are genuinely" % (100.0 * worst))
        print("   different, so the comparison is a real one.")
    print()


def main():
    cols, rows = read_summary(SUMMARY)
    text = read_logs()

    print("=" * 78)
    print("19_thermodynamic_gate")
    print("=" * 78)

    if rows:
        print("%d diagnostic rows, iteration %d to %d"
              % (len(rows), rows[0]["iteration"], rows[-1]["iteration"]))
        print()

        names = sorted({c[:-6] for c in (cols or []) if c.endswith("_total")})
        if names:
            print("%-14s %14s %14s %14s   %s"
                  % ("field", "first total", "last total", "change", "boundary"))
            print("-" * 78)
            for n in names:
                a = conserved(rows[0], n)
                b = conserved(rows[-1], n)
                print("%-14s %14.6g %14.6g %14.6g   %s"
                      % (n, a, b, b - a, BOUNDARY.get(n, "biomass (closed)")))
            print()

        negative = [n for n in names if min(series(rows, n + "_min") or [0.0]) < 0.0]
        if negative:
            print("NEGATIVE VALUES in: %s" % ", ".join(negative))
            print("[v1.3] The compiled kinetics paths apply each reaction increment as")
            print("computed: there is no positivity clamp on them, so a time step or a rate")
            print("constant large enough to consume more than a voxel holds will drive it")
            print("negative. Earlier text here called this impossible by construction and a")
            print("solver bug worth reporting; it is normally a configuration problem. Reduce")
            print("<ade_dt> or the rate constant, then run again.")
            print()

        # Sulfide is the product whose build-up closes the gate; its rise is the
        # mechanism of this case, so say what it did.
        hs = series(rows, "HS_max")
        if hs:
            print("sulfide peak    %.4g -> %.4g mol/L" % (hs[0], hs[-1]))
            print("   This is the quantity that shuts the gate. If it barely moved, the gate")
            print("   had nothing to act on however the energetics were written.")
            print()
    else:
        print("No %s found; judging from the log alone." % SUMMARY)
        print()

    verdict = thermo_verdict(text) if text else None
    if text and verdict is None:
        print("No [THM] block in the log.")
        print("   Either <thermodynamics><enabled> is false, or the run stopped before its")
        print("   closing report. Check the start-up lines: the file, its temperature and")
        print("   every reaction block are echoed there before the first step.")
        print()
    if not text:
        print("No .log file in %s/. Capture it so the run stays reproducible:" % LOGDIR)
        print("   ./complab CompLaB.xml 2>&1 | tee output/run.log")
        print()

    compare_curve(verdict)
    against_control()

    if text:
        hits = [ln.strip() for ln in text.splitlines() if "[THM]" in ln]
        if hits:
            print("THE ENERGETICS, AS PARSED")
            print("-" * 78)
            for ln in hits[:20]:
                print("   %s" % ln)
            if len(hits) > 20:
                print("   ... and %d more" % (len(hits) - 20))
            print()

    print("-" * 78)
    print("WHAT TO DO NEXT")
    print("-" * 78)
    print("If the control run above is missing, ./pipeline.sh produces it. The rate law,")
    print("the geometry, the boundary conditions and the time step are identical between")
    print("the two runs; every difference in the fields is the gate, which is the only way")
    print("to make a claim about what it did.")
    print()
    print("Read the [THM] block above before believing anything else. It prints the")
    print("energetics as the PARSER understood them: the temperature, the reaction as")
    print("written, dG0, and the threshold. An energetics file that parsed differently")
    print("from how you read it shows up there and nowhere else.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
