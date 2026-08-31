#!/usr/bin/env python3
"""
17_symbolic_law  --  POST-PROCESSING

The rate law read from input/growth.sym at start-up, with no rebuild.

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

    acetate      OPEN    (Dirichlet / Neumann), supplied from outside
    o2           OPEN    (Dirichlet / Neumann), supplied from outside
"""
import csv
import os
import re
import sys

SUMMARY = os.path.join("output", "summary.csv")
LOGDIR = "output"

BOUNDARY = {'acetate': 'Dirichlet / Neumann', 'o2': 'Dirichlet / Neumann'}


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
    print("17_symbolic_law")
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

    print("BALANCE CHECK -- none is possible for this case.")
    print("   Every species here is held at a boundary, so no total is closed and")
    print("   any conservation check would fail for reasons that mean nothing.")
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
        for label, pattern in [('the law, AS PARSED', '\\[SYM\\]'), ('diagnostics', '\\[DIAG\\]')]:
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
    print('Read the [SYM] block above before believing anything else. It prints every expression as')
    print('the PARSER understood it, the resolved units, and the valid range of each variable. An')
    print('expression that parsed differently from how you read it shows up there and nowhere else.')
    print('')
    print('Check the units line first. Omit it and the file is read as per_second; a law fitted from an')
    print('hourly sweep is then 3600 times too fast, and nothing warns you.')
    print('')
    print('Then the clamp count. A large fraction of evaluations clamped means the run is spending its')
    print('time outside the box the law was fitted in, and the answer is an extrapolation whatever the')
    print('numbers look like.')
    print('')
    print('Both substrates are fed here, so no total is closed. Compare against example 05, which is the')
    print('same colony driven by compiled kinetics.')

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
