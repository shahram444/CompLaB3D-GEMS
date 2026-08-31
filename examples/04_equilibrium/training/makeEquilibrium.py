#!/usr/bin/env python3
"""
Build the <equilibrium> tableau from chemical formulas.

WHY THIS EXISTS

The speciation solver wants a stoichiometry matrix: one row per species, one
column per component, saying how many of each component that species contains.
Assembling that by hand is error-prone in a way that does not announce itself
-- a wrong row gives a solver that converges happily to the wrong speciation.

You already know the chemistry as formulas. This works the matrix out from
them, by balancing elements and charge, and refuses to guess when the answer is
not unique.

    python3 makeEquilibrium.py carbonate.spc

THE FILE FORMAT

    # components: the master species, one per line, in any order
    component HCO3-1
    component H+1

    # every species, IN CompLaB.xml SUBSTRATE ORDER, with its formation
    # constant. Components appear here too, with logK 0.
    species HCO3-1    logK 0
    species H+1       logK 0
    species CO3-2     logK -10.33
    species H2CO3     logK 6.35

    # optional: a species the solver should ignore, such as a mineral or
    # biomass. It gets a row of zeros.
    inert FeS

FORMULAS

    H2CO3     plain
    CO3-2     charge -2
    Fe+2      charge +2
    Ca+2
    CH3COO-1  organic ligands are fine; the parser only needs elements

Written as element symbols with optional counts, then an optional signed
charge. Parentheses are not supported: write Ca(OH)2 as CaO2H2.

WHAT IT CHECKS

  * that the components are chemically independent. If they are not, no unique
    tableau exists, and the tool says so instead of picking one.
  * that every species can be written in terms of the components at all.
  * that each row balances exactly, elements and charge.

KNOWN LIMIT, stated plainly: the solver these numbers feed uses ideal
activities, with no temperature dependence, no gas phase, no redox couple and
no mineral saturation index. The log K values you supply are conditional
constants at your own ionic strength and 25 C.

Part of the CompLaB program.  GNU Affero General Public License v3 or later.
"""

import argparse
import re
import sys

try:
    import numpy as np
except ImportError:
    sys.exit("ERROR: numpy is required.  python3 -m pip install numpy")

FORMULA = re.compile(r"([A-Z][a-z]?)(\d*)")
CHARGE = re.compile(r"([+-])(\d*)$")


def parse_formula(f, where):
    """'CO3-2' into ({'C':1,'O':3}, -2)."""
    charge = 0
    m = CHARGE.search(f)
    body = f
    if m:
        body = f[:m.start()]
        charge = int(m.group(2) or 1) * (1 if m.group(1) == "+" else -1)
    if "(" in body or ")" in body:
        raise SystemExit("ERROR: %s: '%s' has parentheses, which this parser "
                         "does not handle. Write Ca(OH)2 as CaO2H2." % (where, f))
    counts = {}
    pos = 0
    for mm in FORMULA.finditer(body):
        if mm.start() != pos:
            raise SystemExit("ERROR: %s: cannot read the formula '%s' at "
                             "character %d" % (where, f, pos))
        counts[mm.group(1)] = counts.get(mm.group(1), 0) + int(mm.group(2) or 1)
        pos = mm.end()
    if pos != len(body):
        raise SystemExit("ERROR: %s: cannot read the formula '%s'" % (where, f))
    if not counts:
        raise SystemExit("ERROR: %s: '%s' contains no elements" % (where, f))
    return counts, charge


def name_of(f):
    """A tag-safe name: CO3-2 becomes CO3, Fe+2 becomes Fe.

    The tableau needs names that match <name_of_substrates>, and those cannot
    carry a charge suffix comfortably."""
    m = CHARGE.search(f)
    return f[:m.start()] if m else f


def main():
    ap = argparse.ArgumentParser(
        prog="makeEquilibrium.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Derive the CompLB3D <equilibrium> tableau from formulas.",
        epilog="See the header of this file for the file format.")
    ap.add_argument("spec", help="the .spc species file")
    ap.add_argument("-o", "--output", default=None,
                    help="write the XML block here (default: print it)")
    args = ap.parse_args()

    comps, species, inert = [], [], []
    with open(args.spec) as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.split("#")[0].strip()
            if not line:
                continue
            parts = line.split()
            where = "line %d" % lineno
            if parts[0] == "component":
                comps.append(parts[1])
            elif parts[0] == "species":
                logk = 0.0
                if "logK" in parts:
                    logk = float(parts[parts.index("logK") + 1])
                species.append((parts[1], logk, where))
            elif parts[0] == "inert":
                species.append((parts[1], 0.0, where))
                inert.append(parts[1])
            else:
                raise SystemExit("ERROR: %s: '%s' should start with component, "
                                 "species or inert" % (where, parts[0]))

    if not comps:
        raise SystemExit("ERROR: no 'component' lines.")
    if not species:
        raise SystemExit("ERROR: no 'species' lines.")

    # ---- element basis ----------------------------------------------------
    parsed = {}
    elements = []
    for f in comps + [s[0] for s in species if s[0] not in inert]:
        c, q = parse_formula(f, f)
        parsed[f] = (c, q)
        for e in c:
            if e not in elements:
                elements.append(e)
    elements.sort()
    rows = elements + ["charge"]

    def vec(f):
        c, q = parsed[f]
        return np.array([c.get(e, 0) for e in elements] + [q], float)

    A = np.column_stack([vec(c) for c in comps])          # (nrow, ncomp)

    rank = np.linalg.matrix_rank(A)
    if rank < len(comps):
        raise SystemExit(
            "ERROR: the %d components are not chemically independent (rank %d).\n"
            "       No unique tableau exists, so this tool will not guess one.\n"
            "       Drop a component, or pick a different set."
            % (len(comps), rank))

    print("Components (%d): %s" % (len(comps), " ".join(comps)))
    print("Elements balanced: %s" % " ".join(rows))
    print()

    table = []
    for f, logk, where in species:
        if f in inert:
            table.append((name_of(f), [0.0] * len(comps), 0.0, "inert"))
            continue
        b = vec(f)
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        resid = A @ x - b
        if np.abs(resid).max() > 1e-8:
            bad = [rows[i] for i in range(len(rows)) if abs(resid[i]) > 1e-8]
            raise SystemExit(
                "ERROR: %s: '%s' cannot be written in terms of the components.\n"
                "       Unbalanced: %s\n"
                "       Either a component is missing, or the formula is wrong."
                % (where, f, ", ".join(bad)))
        x = np.where(np.abs(x) < 1e-10, 0.0, x)
        table.append((name_of(f), list(x), logk, ""))

    # ---- report -----------------------------------------------------------
    w = max(len(t[0]) for t in table) + 2
    print("%-*s %s   logK" % (w, "species", "  ".join("%6s" % name_of(c) for c in comps)))
    for nm, row, logk, tag in table:
        print("%-*s %s   %-8g %s"
              % (w, nm, "  ".join("%6g" % v for v in row), logk, tag))

    # ---- XML --------------------------------------------------------------
    L = ["    <equilibrium>", "        <enabled>true</enabled>", "",
         "        <components>%s</components>"
         % " ".join(name_of(c) for c in comps), "",
         "        <stoichiometry>"]
    for i, (nm, row, logk, tag) in enumerate(table):
        L.append("            <species%d>%s</species%d>   <!-- %s -->"
                 % (i, " ".join("%g" % v for v in row), i, nm))
    L += ["        </stoichiometry>", "", "        <logK>"]
    for i, (nm, row, logk, tag) in enumerate(table):
        L.append("            <species%d>%g</species%d>" % (i, logk, i))
    L += ["        </logK>", "    </equilibrium>"]
    xml = "\n".join(L) + "\n"

    print()
    if args.output:
        with open(args.output, "w") as f:
            f.write(xml)
        print("Wrote %s" % args.output)
    else:
        print(xml)

    print("Paste that into CompLaB.xml, and make sure <name_of_substrates> "
          "reads, in this order:")
    print("    " + " ".join(t[0] for t in table))
    print("\nA species with a row of zeros is ignored by the solver. That is "
          "how minerals and biomass stay out of the tableau.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
