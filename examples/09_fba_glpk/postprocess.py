#!/usr/bin/env python3
"""
09_fba_glpk  --  POST-PROCESSING

Growth predicted by a genome-scale linear program, solved in process by GLPK.

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

    S            OPEN    (Dirichlet / Neumann), supplied from outside
    O            OPEN    (Dirichlet / Neumann), supplied from outside
    P            closed  (Neumann both ends), so it can be balanced
"""
import csv
import os
import re
import sys

SUMMARY = os.path.join("output", "summary.csv")
LOGDIR = "output"

BOUNDARY = {'S': 'Dirichlet / Neumann', 'O': 'Dirichlet / Neumann', 'P': 'closed'}


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
    print("09_fba_glpk")
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

    print("CLOSED-SPECIES REPORT -- P, the product")
    for line in 'S and O are both fed from the left.'.split("\n"):
        print("   %s" % line)
    terms = [('P', 1.0)]
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
        for label, pattern in [('diagnostics', '\\[DIAG\\]'), ('the linear program', '(?i)\\[FBA\\]|infeasible|no growth|warm.start')]:
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
    print('Growth should appear only where BOTH substrates are present. Where either runs out the')
    print('linear program has no feasible optimum and the organism falls back to decaying at')
    print('<decay_coefficient>.')
    print('')
    print('A biomass total that falls everywhere usually means the uptake bounds are zero, not that the')
    print('model cannot grow. Check <maximum_uptake_flux> before you suspect the model.')

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
