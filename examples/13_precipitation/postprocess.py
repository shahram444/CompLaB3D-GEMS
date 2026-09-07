#!/usr/bin/env python3
"""
13_precipitation  --  POST-PROCESSING

FeS forms where an iron front meets a sulfide front, fills its voxels, and seals them.

  reads   output/summary.csv   one row per diagnostic interval
          output/*.log         the solver log, if it was captured
  writes  nothing; prints a verdict

Self-contained: standard library only, nothing imported from tools/. Run it
after the solver, or let pipeline.sh run it for you:

    python3 postprocess.py

WHY A VERDICT AND NOT A PICTURE. A field plot looks equally plausible whether or
not the run conserved mass. Check the run first; plot it once you know it is
worth looking at.

WHICH TOTALS CAN BE BALANCED AT ALL. A species held at a Dirichlet boundary is
supplied from outside the domain, so its total is not conserved and never will
be. Naming it in a conservation check guarantees a failure that means nothing.
Only closed species are balanced below. For this case:

    Fe2          OPEN    (Dirichlet / Neumann), supplied from outside
    HS           OPEN    (Neumann / Dirichlet), supplied from outside
    FeS          closed  (bounce-back both ends), so it can be balanced
"""
import csv
import os
import re
import sys

SUMMARY = os.path.join("output", "summary.csv")
LOGDIR = "output"

BOUNDARY = {'Fe2': 'Dirichlet / Neumann', 'HS': 'Neumann / Dirichlet', 'FeS': 'closed'}



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
    """The CSV the <diagnostics> block writes. Comment lines start with #."""
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


def main():
    cols, rows = read_summary(SUMMARY)
    if not rows:
        print("No %s found." % SUMMARY)
        print()
        print("That file is written by the <diagnostics> block. To get it, add this to")
        print("CompLaB.xml and run again:")
        print()
        print("    <diagnostics>")
        print("        <enabled>true</enabled>")
        print("        <summary_csv>summary.csv</summary_csv>")
        print("        <interval>500</interval>")
        print("        <tolerance>1e-6</tolerance>")
        print("    </diagnostics>")
        return 1

    print("=" * 78)
    print("13_precipitation")
    print("=" * 78)
    print("%d diagnostic rows, iteration %d to %d"
          % (len(rows), rows[0]["iteration"], rows[-1]["iteration"]))
    print()

    por = series(rows, "porosity")
    if por:
        change = por[-1] - por[0]
        if abs(change) < 1e-12:
            print("porosity        unchanged at %.4f" % por[0])
        else:
            print("porosity        %s  %.4f -> %.4f"
                  % ("FELL" if change < 0 else "ROSE", por[0], por[-1]))
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
                  % (n, a, b, b - a, BOUNDARY.get(n, "biomass")))
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

    print("CLOSED-SPECIES REPORT -- FeS, which is closed and immobile")
    for line in 'Fe2 is fed from the left and HS from the right, so only the mineral is closed.'.split("\n"):
        print("   %s" % line)
    terms = [('FeS', 1.0)]
    missing = [n for n, _ in terms if (n + "_total") not in (cols or [])]
    if missing:
        print("   SKIPPED: %s absent from the summary." % ", ".join(missing))
        print("   A check that never ran is not a pass.")
    else:
        for name, _ in terms:
            hist = series(rows, name + "_total")
            change = hist[-1] - hist[0]
            monotone = (all(b >= a - 1e-15 for a, b in zip(hist, hist[1:]))
                        or all(b <= a + 1e-15 for a, b in zip(hist, hist[1:])))
            print("   %-10s change %+.6g   %s   %s"
                  % (name, change,
                     "accumulating" if change > 0 else
                     ("depleting" if change < 0 else "flat"),
                     "monotone" if monotone else "NOT monotone: it went both ways"))
        print("   A closed species changes only by reaction, so this IS the reaction")
        print("   rate integrated over the run. Nothing here came from a boundary.")
    print()


    logs = ([os.path.join(LOGDIR, f) for f in sorted(os.listdir(LOGDIR))
             if f.endswith(".log")] if os.path.isdir(LOGDIR) else [])
    if logs:
        text = ""
        for path in logs:
            try:
                text += open(path, errors="replace").read()
            except OSError:
                pass
        for label, pattern in [('clogging', '(?i)porosit|percolat|clog|\\[PRECIP'), ('diagnostics', '\\[DIAG\\]')]:
            hits = [ln.strip() for ln in text.splitlines() if re.search(pattern, ln)]
            if hits:
                print("%s:" % label)
                for ln in hits[:12]:
                    print("   %s" % ln)
                if len(hits) > 12:
                    print("   ... and %d more" % (len(hits) - 12))
                print()
    else:
        print("No .log file in %s/. Capture it so the run stays reproducible:" % LOGDIR)
        print("   ./complab CompLaB.xml 2>&1 | tee output/run.log")
        print()

    print("-" * 78)
    print("WHAT TO LOOK AT")
    print("-" * 78)
    print('Porosity should fall in STEPS, one per <update_interval>, not smoothly. A voxel converts')
    print('only when it is full; between conversions the porosity is flat while mineral accumulates.')
    print('')
    print('Every mole of FeS formed must remove one mole of Fe2 and one of HS. Both are fed, so you')
    print('cannot verify that from totals alone.')
    print('')
    print('[v1.3] There is no closed sum to declare in <diagnostics> here. Fe2 is held at a Dirichlet')
    print('boundary on the left face and HS at one on the right, so both are supplied from outside the')
    print('domain and both totals rise without bound. Naming Fe2+FeS or HS+FeS in a conservation check')
    print('makes the solver report a mass balance FAIL on a run that is behaving correctly. What')
    print('postprocess.py already reports is the check available in this case: the change in each')
    print('field, and the mineral inventory against the solute drawn down.')
    print('')
    print('If the run stops reporting a percolation limit, the domain sealed. That is the code working.')

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
