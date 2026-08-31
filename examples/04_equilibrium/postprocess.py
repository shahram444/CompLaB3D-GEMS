#!/usr/bin/env python3
"""
04_equilibrium  --  POST-PROCESSING

Aqueous speciation solved in every pore voxel, every step.

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

    HCO3         OPEN    (Dirichlet / Neumann), supplied from outside
    H            OPEN    (Dirichlet / Neumann), supplied from outside
    CO3          closed  (Neumann both ends), so it can be balanced
    H2CO3        closed  (Neumann both ends), so it can be balanced
"""
import csv
import os
import re
import sys

SUMMARY = os.path.join("output", "summary.csv")
LOGDIR = "output"

BOUNDARY = {'HCO3': 'Dirichlet / Neumann', 'H': 'Dirichlet / Neumann', 'CO3': 'closed', 'H2CO3': 'closed'}


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
    print("04_equilibrium")
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

    print("CLOSED-SPECIES REPORT -- the complexes, which are formed in place")
    for line in 'HCO3 and H are fed; the two complexes are made from them inside the domain.'.split("\n"):
        print("   %s" % line)
    terms = [('CO3', 1.0), ('H2CO3', 1.0)]
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
        for label, pattern in [('diagnostics', '\\[DIAG\\]')]:
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
    print("Compare a component's FREE concentration with its TOTAL. They are equal only when no")
    print('complex carries it; once CO3 and H2CO3 appear the two differ, and that difference is exactly')
    print('what <fba_concentration_basis>total</fba_concentration_basis> exists to account for.')
    print('')
    print('Speciation here is ideal and isothermal: no activity correction, no temperature dependence,')
    print('no redox couple, no saturation index. Your log K values are conditional constants at your own')
    print('ionic strength and 25 C, and nothing in the output will remind you of that.')

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
