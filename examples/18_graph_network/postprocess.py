#!/usr/bin/env python3
"""
18_graph_network  --  POST-PROCESSING

All four rates plus growth from one evaluation of a message-passing network.

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

    CH4          OPEN    (Dirichlet / Neumann), supplied from outside
    SO4          OPEN    (Dirichlet / Neumann), supplied from outside
    HS           closed  (Neumann both ends), so it can be balanced
    HCO3         closed  (Neumann both ends), so it can be balanced
"""
import csv
import os
import re
import sys

SUMMARY = os.path.join("output", "summary.csv")
LOGDIR = "output"

BOUNDARY = {'CH4': 'Dirichlet / Neumann', 'SO4': 'Dirichlet / Neumann', 'HS': 'closed', 'HCO3': 'closed'}


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
    print("18_graph_network")
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

    print("BALANCE CHECK -- HS against HCO3: the stoichiometric ratio")
    for line in 'AOM makes one HS and one HCO3 per methane, so the DIFFERENCE of the two changes should be zero.\nBoth are closed; CH4 and SO4 are fed from the left.'.split("\n"):
        print("   %s" % line)
    terms = [('HS', 1.0), ('HCO3', -1.0)]
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
        for label, pattern in [('the network provenance', '\\[GNN\\]'), ('diagnostics', '\\[DIAG\\]')]:
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
    print('THE NUMBER THAT MATTERS IS THE RATIO, not any single total. AOM is')
    print('')
    print('    CH4 + SO4  ->  HS + HCO3')
    print('')
    print('so HS and HCO3 must be produced one for one, and the check above tests exactly that. If they')
    print('diverge, the stoichiometric matrix given to the trainer is wrong -- NOT the fit. The network')
    print('cannot break a ratio it was handed as structure, which is the whole argument for this path')
    print('over a symbolic law with the species written out by hand.')
    print('')
    print('Then the clamp count in the [GNN] block. trainmin and trainmax travel inside the .gnn file and')
    print('are enforced at every evaluation. A count that rises through the run means the colony has')
    print('walked the concentrations out of the box the network was trained in.')

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
