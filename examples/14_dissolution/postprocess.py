#!/usr/bin/env python3
"""
14_dissolution  --  POST-PROCESSING

Calcite is consumed by the water touching it and the voxel reopens once enough has gone.

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

    H            OPEN    (Dirichlet / Neumann), supplied from outside
    Ca           closed  (Neumann both ends), so it can be balanced
    calcite      closed  (Neumann both ends), so it can be balanced
"""
import csv
import os
import re
import sys

SUMMARY = os.path.join("output", "summary.csv")
LOGDIR = "output"

BOUNDARY = {'H': 'Dirichlet / Neumann', 'Ca': 'closed', 'calcite': 'closed'}


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
    print("14_dissolution")
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
            a = rows[0].get(n + "_total", 0.0)
            b = rows[-1].get(n + "_total", 0.0)
            print("%-14s %14.6g %14.6g %14.6g   %s"
                  % (n, a, b, b - a, BOUNDARY.get(n, "biomass")))
        print()

    negative = [n for n in names if min(series(rows, n + "_min") or [0.0]) < 0.0]
    if negative:
        print("NEGATIVE VALUES in: %s" % ", ".join(negative))
        print("The mass-budget clamp makes this impossible by construction, so this is")
        print("a solver bug rather than a configuration problem. Worth reporting.")
        print()

    print("BALANCE CHECK -- calcite lost against calcium released")
    for line in 'Both are closed: only H is fed. Every mole of calcite dissolved must appear as one mole of Ca,\nso the SUM of the two changes should be zero to the diagnostics tolerance.'.split("\n"):
        print("   %s" % line)
    terms = [('calcite', 1.0), ('Ca', 1.0)]
    missing = [n for n, _ in terms if (n + "_total") not in (cols or [])]
    if missing:
        print("   SKIPPED: %s absent from the summary." % ", ".join(missing))
        print("   A check that never ran is not a pass.")
    else:
        net = sum(c * (rows[-1][n + "_total"] - rows[0][n + "_total"]) for n, c in terms)
        scale = max(abs(rows[-1][n + "_total"] - rows[0][n + "_total"])
                    for n, _ in terms) or 1.0
        rel = abs(net) / scale
        print("   residual %+.6g   relative to the largest change %.3g   %s"
              % (net, rel, "PASS" if rel <= 1e-3 else "FAIL"))
        if rel > 1e-3:
            print("   The two are not tracking each other. Suspect the stoichiometry")
            print("   before the numerics.")
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
        for label, pattern in [('reopening', '(?i)porosit|reopen|\\[DISSOL'), ('diagnostics', '\\[DIAG\\]')]:
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
        print("   ./build/complab CompLaB.xml 2>&1 | tee output/run.log")
        print()

    print("-" * 78)
    print("WHAT TO LOOK AT")
    print("-" * 78)
    print('RUN THIS ON ONE PROCESS. A dissolving voxel deposits its products into its open face')
    print('neighbours, and where a neighbour lies in a different MPI block the deposit is written into')
    print('the local envelope and never communicated back. The result therefore depends on the processor')
    print('count, with the error concentrated at block interfaces. If the balance above fails on several')
    print('ranks and passes on one, that is this limitation and not your chemistry.')
    print('')
    print('Porosity should RISE, in steps, for the same reason 13 falls in steps. A voxel reopens below')
    print('<reopen_fraction> of full rather than at zero, so a voxel sitting on the threshold cannot')
    print('flicker between open and closed.')

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
