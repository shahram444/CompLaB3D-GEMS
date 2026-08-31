#!/usr/bin/env python3
"""
Write defineKinetics.hh and defineAbioticKinetics.hh from a reaction list.

WHY THIS EXISTS

CompLB3D compiles its rate laws in, which means that until now the only way to
describe your chemistry was to write C++ by hand: get the substrate indices
right, get the signs right, get the per-second conversion right, remember to
zero the output arrays, remember to floor the negatives, and remember the guard
on array size. Every one of those is a silent failure if you get it wrong.

This generates all of it from a plain reaction file, and generates a self-test
alongside so you can prove mass balance before the code goes anywhere near
Palabos.

    python3 makeKinetics.py mychemistry.rxn -o .

WHAT IT COVERS

    abiotic mass action     a A + b B -> c C,   rate = k [A]^a [B]^b
    microbial growth        multi-Monod on any number of substrates,
                            with yields, product release and first-order decay
    dissolution             first order in one solute and in remaining mineral

Anything outside that, write by hand. The generated file is ordinary readable
C++ and is a good starting point to edit.

THE FILE FORMAT

    # comments start with a hash

    substrates: A B C           # in CompLaB.xml order. Mandatory.
    microbes:   Bug             # in CompLaB.xml order. Omit if abiotic only.

    abiotic R1:                 # a mass-action reaction
        equation: A + B -> C
        k: 5.0e-2               # per second, in whatever order the equation is

    microbe Bug:                # a growing organism
        mu_max: 2.0e-4          # 1/s
        monod: A 0.05           # substrate, half-saturation. Repeatable.
        yield: A 0.4            # gDW of cells per mol of A consumed
        release: C 2.0          # mol of C per gDW of cells made. Repeatable.
        decay: 1.0e-6           # 1/s, first order

    dissolve calcite:           # a dissolution phase, numbered in <phaseN> order
        equation: calcite + H -> Ca
        k: 1.0e-2               # per second per mol/L of remaining mineral
        driver: H               # which solute drives it

UNITS, everywhere: rates PER SECOND, concentrations per litre, biomass gDW/L.
A rate constant published per day must be divided by 86400 before it goes in.

Part of the CompLaB program.  GNU Affero General Public License v3 or later.
"""

import argparse
import os
import re
import sys


class SpecError(Exception):
    pass


# ===========================================================================
def parse(path):
    """Read the reaction file into a dict. Errors name the line number."""
    spec = {"substrates": [], "microbes": [], "abiotic": [], "microbe": [],
            "dissolve": []}
    cur = None
    with open(path) as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.split("#")[0].rstrip()
            if not line.strip():
                continue
            indented = line[0] in " \t"
            line = line.strip()

            if not indented:
                if line.startswith("substrates:"):
                    spec["substrates"] = line.split(":", 1)[1].split()
                    cur = None
                elif line.startswith("microbes:"):
                    spec["microbes"] = line.split(":", 1)[1].split()
                    cur = None
                elif line.startswith(("abiotic ", "microbe ", "dissolve ")):
                    kind, name = line.split(None, 1)
                    name = name.rstrip(":").strip()
                    cur = {"name": name, "line": lineno}
                    spec[kind].append(cur)
                else:
                    raise SpecError(
                        "line %d: '%s' is not a block header.\n"
                        "  Expected one of: substrates:, microbes:, "
                        "'abiotic <name>:', 'microbe <name>:', "
                        "'dissolve <name>:'" % (lineno, line))
                continue

            if cur is None:
                raise SpecError("line %d: '%s' is indented but no block is "
                                "open above it" % (lineno, line))
            if ":" not in line:
                raise SpecError("line %d: '%s' has no colon; every setting is "
                                "'key: value'" % (lineno, line))
            key, val = line.split(":", 1)
            key, val = key.strip(), val.strip()
            if key in ("monod", "yield", "release", "inhibit"):
                cur.setdefault(key, []).append(val)
            else:
                cur[key] = val
    if not spec["substrates"]:
        raise SpecError("no 'substrates:' line. It is mandatory and must list "
                        "the substrate names in CompLaB.xml order.")
    return spec


def idx_of(spec, name, where):
    if name not in spec["substrates"]:
        raise SpecError("%s names '%s', which is not in the substrates list "
                        "(%s)" % (where, name, " ".join(spec["substrates"])))
    return spec["substrates"].index(name)


def parse_equation(eq, spec, where):
    """'2 A + B -> C' into [(coeff, index, name), ...] with reactants negative."""
    if "->" not in eq:
        raise SpecError("%s: equation '%s' has no '->'" % (where, eq))
    lhs, rhs = eq.split("->", 1)
    terms = []
    for side, sign in ((lhs, -1.0), (rhs, +1.0)):
        for part in side.split("+"):
            part = part.strip()
            if not part:
                continue
            m = re.match(r"^(?:([0-9.eE+-]+)\s+)?([A-Za-z_][A-Za-z_0-9]*)$", part)
            if not m:
                raise SpecError("%s: cannot read '%s' in '%s'. Write terms as "
                                "'A' or '2 A'." % (where, part, eq))
            coeff = float(m.group(1)) if m.group(1) else 1.0
            name = m.group(2)
            terms.append((sign * coeff, idx_of(spec, name, where), name))
    return terms


# ===========================================================================
LICENCE = """/* This file is a part of the CompLaB program.
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/
"""


def banner(spec, source, what):
    L = [LICENCE, "", "/* " + "=" * 74,
         " * %s  --  GENERATED FILE, DO NOT EDIT BY HAND." % what,
         " *",
         " * Produced by makeKinetics.py from %s" % os.path.basename(source),
         " *",
         " * SUBSTRATE INDICES, in CompLaB.xml order. Reordering the substrates",
         " * in the XML changes what every line below means, so regenerate this",
         " * file whenever you touch that list."]
    for i, s in enumerate(spec["substrates"]):
        L.append(" *     C[%d] = %s" % (i, s))
    if spec["microbes"]:
        L.append(" *")
        for i, m in enumerate(spec["microbes"]):
            L.append(" *     B[%d] = %s" % (i, m))
    L += [" *",
          " * UNITS: every rate is PER SECOND. Concentrations are per litre and",
          " * biomass is gDW/L. A constant published per day was divided by 86400",
          " * before it got here -- check that it was.",
          " * " + "=" * 74 + " */"]
    return "\n".join(L)


def gen_abiotic(spec, source):
    out = [banner(spec, source, "defineAbioticKinetics.hh"),
           "#ifndef DEFINE_ABIOTIC_KINETICS_HH",
           "#define DEFINE_ABIOTIC_KINETICS_HH", "",
           "#include <vector>", "#include <cstddef>", "#include <cmath>",
           "#include <algorithm>", ""]
    ns = len(spec["substrates"])

    if spec["abiotic"]:
        out.append("namespace GeneratedAbiotic {")
        for r in spec["abiotic"]:
            if "k" not in r:
                raise SpecError("abiotic %s: no 'k'" % r["name"])
            out.append("    const double k_%s = %s;" % (r["name"], r["k"]))
        out.append("}")
        out.append("")

    out += ["void defineAbioticRxnKinetics(std::vector<double> C,",
            "                              std::vector<double>& subsR,",
            "                              plb::plint mask)",
            "{",
            "    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;",
            "    if (C.size() < %d || subsR.size() < %d) return;" % (ns, ns),
            "    if (mask < 2) return;   // nothing happens in solid or wall voxels",
            ""]

    if not spec["abiotic"]:
        out.append("    // no abiotic reactions were declared")
    for r in spec["abiotic"]:
        where = "abiotic %s (line %d)" % (r["name"], r["line"])
        terms = parse_equation(r["equation"], spec, where)
        out.append("    // %s :  %s" % (r["name"], r["equation"]))
        reactants = [(c, i, n) for (c, i, n) in terms if c < 0]
        rate = "GeneratedAbiotic::k_%s" % r["name"]
        for c, i, n in reactants:
            order = int(round(-c))
            for _ in range(order):
                rate += " * std::max(C[%d], 0.0)" % i
        out.append("    {")
        out.append("        const double R = %s;" % rate)
        for c, i, n in terms:
            out.append("        subsR[%d] += %s%g * R;   // %s"
                       % (i, "+" if c > 0 else "-", abs(c), n))
        out.append("    }")
        out.append("")

    out += ["    for (std::size_t i = 0; i < subsR.size(); ++i)",
            "        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;",
            "}", ""]

    # ---- dissolution -------------------------------------------------------
    if spec["dissolve"]:
        out.append("namespace GeneratedDissolution {")
        for i, d in enumerate(spec["dissolve"], start=1):
            out.append("    const double k_%s = %s;   // phase %d"
                       % (d["name"], d.get("k", "0.0"), i))
        out.append("}")
        out.append("")

    out += ["inline void defineDissolutionRate(plb::plint phaseId,",
            "                                  const std::vector<double>& C,",
            "                                  double mineral,",
            "                                  std::vector<double>& subsR,",
            "                                  double& mineralR,",
            "                                  plb::plint mask)",
            "{",
            "    (void) mask;",
            "    mineralR = 0.0;",
            "    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;",
            "    if (C.size() < %d || subsR.size() < %d) return;" % (ns, ns),
            ""]
    if not spec["dissolve"]:
        out.append("    (void) phaseId; (void) mineral;")
        out.append("    // no dissolving phases were declared")
    for i, d in enumerate(spec["dissolve"], start=1):
        where = "dissolve %s (line %d)" % (d["name"], d["line"])
        terms = parse_equation(d["equation"], spec, where)
        driver = d.get("driver")
        if driver is None:
            raise SpecError("%s: no 'driver'. Name the solute whose "
                            "concentration drives the rate." % where)
        di = idx_of(spec, driver, where)
        mineral_names = [n for (c, ix, n) in terms if c < 0 and n != driver]
        if not mineral_names:
            raise SpecError("%s: the left side must contain the mineral as well "
                            "as the driver" % where)
        out += ["    if (phaseId == %d) {          // %s :  %s"
                % (i, d["name"], d["equation"]),
                "        const double drive = std::max(C[%d], 0.0);   // %s" % (di, driver),
                "        const double rate  = GeneratedDissolution::k_%s * drive"
                " * std::max(mineral, 0.0);" % d["name"],
                "        mineralR = -rate;              // the mineral disappears"]
        for c, ix, n in terms:
            if n in mineral_names:
                continue        # the solver handles the mineral from mineralR
            out.append("        subsR[%d] = %s%g * rate;   // %s"
                       % (ix, "+" if c > 0 else "-", abs(c), n))
        out.append("    }")
    out += ["",
            "    if (std::isnan(mineralR) || std::isinf(mineralR)) mineralR = 0.0;",
            "    for (std::size_t i = 0; i < subsR.size(); ++i)",
            "        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;",
            "}", "", "#endif", ""]
    return "\n".join(out)


def gen_biotic(spec, source):
    ns, nm = len(spec["substrates"]), len(spec["microbes"])
    out = [banner(spec, source, "defineKinetics.hh"),
           "#ifndef DEFINE_KINETICS_HH", "#define DEFINE_KINETICS_HH", "",
           "#include <vector>", "#include <cstddef>", "#include <cmath>",
           "#include <algorithm>", ""]

    if spec["microbe"]:
        out.append("namespace GeneratedMonod {")
        for m in spec["microbe"]:
            n = m["name"]
            out.append("    // %s" % n)
            out.append("    const double mu_max_%s = %s;   // 1/s" % (n, m.get("mu_max", "0.0")))
            out.append("    const double decay_%s  = %s;   // 1/s" % (n, m.get("decay", "0.0")))
        out.append("}")
        out.append("")

    out += ["inline double gm_monod(double S, double K) { return S / (K + S); }",
            "",
            "void defineRxnKinetics(std::vector<double> B, std::vector<double> C,",
            "                       std::vector<double>& subsR, std::vector<double>& bioR,",
            "                       plb::plint mask)",
            "{",
            "    for (std::size_t i = 0; i < subsR.size(); ++i) subsR[i] = 0.0;",
            "    for (std::size_t i = 0; i < bioR.size();  ++i) bioR[i]  = 0.0;",
            "    if (C.size() < %d || subsR.size() < %d) return;" % (ns, ns),
            "    if (mask < 2) return;   // no biology in solid or wall voxels",
            ""]
    if not spec["microbe"]:
        out.append("    (void) B;")
        out.append("    // no microbes were declared")

    for m in spec["microbe"]:
        n = m["name"]
        where = "microbe %s (line %d)" % (n, m["line"])
        if n not in spec["microbes"]:
            raise SpecError("%s: '%s' is not in the microbes list (%s)"
                            % (where, n, " ".join(spec["microbes"]) or "empty"))
        mi = spec["microbes"].index(n)
        out.append("    // ---- %s ----" % n)
        out.append("    if (B.size() > %d && bioR.size() > %d) {" % (mi, mi))
        out.append("        const double Bm = std::max(B[%d], 0.0);" % mi)
        out.append("        if (Bm > 0.0) {")

        limit = ["GeneratedMonod::mu_max_%s" % n]
        for entry in m.get("monod", []):
            parts = entry.split()
            if len(parts) != 2:
                raise SpecError("%s: 'monod: %s' should be 'monod: <substrate> "
                                "<half saturation>'" % (where, entry))
            s, K = parts
            limit.append("gm_monod(std::max(C[%d], 0.0), %s)"
                         % (idx_of(spec, s, where), K))
        for entry in m.get("inhibit", []):
            parts = entry.split()
            if len(parts) != 2:
                raise SpecError("%s: 'inhibit: %s' should be 'inhibit: "
                                "<substrate> <constant>'" % (where, entry))
            s, K = parts
            limit.append("(%s / (%s + std::max(C[%d], 0.0)))"
                         % (K, K, idx_of(spec, s, where)))
        out.append("            const double growth = " +
                   ("\n                                * ".join(limit)) +
                   "\n                                * Bm;   // gDW/L/s")
        out.append("")
        out.append("            bioR[%d] += growth - GeneratedMonod::decay_%s * Bm;"
                   % (mi, n))
        for entry in m.get("yield", []):
            parts = entry.split()
            if len(parts) != 2:
                raise SpecError("%s: 'yield: %s' should be 'yield: <substrate> "
                                "<gDW per mol>'" % (where, entry))
            s, Y = parts
            if float(Y) == 0.0:
                raise SpecError("%s: a yield of zero would divide by zero" % where)
            out.append("            subsR[%d] -= growth / %s;   // %s consumed"
                       % (idx_of(spec, s, where), Y, s))
        for entry in m.get("release", []):
            parts = entry.split()
            if len(parts) != 2:
                raise SpecError("%s: 'release: %s' should be 'release: "
                                "<substrate> <mol per gDW>'" % (where, entry))
            s, f = parts
            out.append("            subsR[%d] += growth * %s;   // %s released"
                       % (idx_of(spec, s, where), f, s))
        out.append("        }")
        out.append("    }")
        out.append("")

    out += ["    for (std::size_t i = 0; i < subsR.size(); ++i)",
            "        if (std::isnan(subsR[i]) || std::isinf(subsR[i])) subsR[i] = 0.0;",
            "    for (std::size_t i = 0; i < bioR.size(); ++i)",
            "        if (std::isnan(bioR[i]) || std::isinf(bioR[i])) bioR[i] = 0.0;",
            "}", "", "#endif", ""]
    return "\n".join(out)


def _left_nullspace(N, tol=1e-10):
    """Vectors w with w . N = 0, i.e. the conserved moieties.

    N is (n_substrates, n_reactions). Any w orthogonal to every column is a
    quantity no reaction in the set can create or destroy, whatever the rate
    constants are. Checking w . subsR == 0 is therefore the correct, general
    statement of mass balance for a reaction network -- and unlike a
    per-reaction ratio check it stays valid when several reactions share a
    substrate, which was the bug in the first version of this function.
    """
    import numpy as np
    N = np.asarray(N, float)
    if N.size == 0 or N.shape[1] == 0:
        return np.eye(N.shape[0])
    u, s, vt = np.linalg.svd(N)
    rank = int((s > tol * max(1.0, s[0])).sum())
    W = u[:, rank:].T                      # rows span the left null space
    # Tidy the basis so the generated C++ reads like chemistry rather than
    # like an SVD: scale each row so its smallest non-zero entry is 1.
    out = []
    for w in W:
        w = np.where(np.abs(w) < tol, 0.0, w)
        nz = np.abs(w[np.abs(w) > 0])
        if nz.size:
            w = w / nz.min()
        out.append(w)
    return np.array(out) if out else np.zeros((0, N.shape[0]))


def gen_selftest(spec, source):
    """A driver that proves the generated abiotic rates conserve mass.

    Compiles and runs with no Palabos and no CompLB3D: a closed test of the
    chemistry alone, which is the part a generator can get wrong."""
    try:
        import numpy as np                                   # noqa: F401
    except ImportError:
        return None
    ns = len(spec["substrates"])

    N = [[0.0] * len(spec["abiotic"]) for _ in range(ns)]
    for j, r in enumerate(spec["abiotic"]):
        for c, i, _n in parse_equation(r["equation"], spec, "selftest"):
            N[i][j] += c
    W = _left_nullspace(N)

    lines = [LICENCE, "",
             "/* Self-test for the generated abiotic kinetics.",
             " *",
             " * Feeds random states to defineAbioticRxnKinetics and checks every",
             " * CONSERVED MOIETY of the declared reaction network: a combination of",
             " * substrates that no reaction in the set can create or destroy, whatever",
             " * the rate constants are. Those combinations are the left null space of",
             " * the stoichiometry matrix, computed when this file was generated.",
             " *",
             " * This is the right check even when several reactions share a substrate,",
             " * where a per-reaction ratio test would be meaningless because the rates",
             " * are sums.",
             " *",
             " * It needs neither Palabos nor a working build of the solver:",
             " *     g++ -O2 -std=c++11 -o selftest selftest_kinetics.cpp && ./selftest",
             " */",
             "#include <vector>", "#include <cstdio>", "#include <cmath>",
             "namespace plb { typedef long long plint; }",
             '#include "defineAbioticKinetics.hh"', "",
             "int main() {",
             "    const int NS = %d;" % ns,
             "    const int NW = %d;   // conserved moieties" % len(W)]

    if len(W):
        rows = []
        for w in W:
            rows.append("        { " + ", ".join("%.17g" % v for v in w) + " }")
        lines.append("    static const double W[%d][%d] = {" % (len(W), ns))
        lines.append(",\n".join(rows))
        lines.append("    };")
        lines.append("    static const char* WNAME[%d] = {" % len(W))
        names = []
        for w in W:
            terms = ["%+g %s" % (v, spec["substrates"][i])
                     for i, v in enumerate(w) if abs(v) > 1e-12]
            names.append('        "%s"' % " ".join(terms).replace('"', "'"))
        lines.append(",\n".join(names))
        lines.append("    };")
    else:
        lines.append("    static const double (*W)[%d] = 0; (void) W;" % ns)
        lines.append("    static const char** WNAME = 0; (void) WNAME;")

    lines += ["",
              "    unsigned seed = 12345u;",
              "    double worst = 0.0; int trials = 20000; int worstw = -1;",
              "    for (int t = 0; t < trials; ++t) {",
              "        std::vector<double> C(NS), R(NS);",
              "        for (int i = 0; i < NS; ++i) {",
              "            seed = seed * 1103515245u + 12345u;",
              "            C[i] = ((seed >> 16) & 0x7fff) / 32767.0 * 10.0;",
              "        }",
              "        defineAbioticRxnKinetics(C, R, 2);",
              "        double scale = 0.0;",
              "        for (int i = 0; i < NS; ++i) scale += std::fabs(R[i]);",
              "        if (scale < 1e-300) continue;",
              "        for (int k = 0; k < NW; ++k) {",
              "            double dot = 0.0;",
              "            for (int i = 0; i < NS; ++i) dot += W[k][i] * R[i];",
              "            double d = std::fabs(dot) / scale;",
              "            if (d > worst) { worst = d; worstw = k; }",
              "        }",
              "    }",
              '    std::printf("%d trials, %d conserved moiety(ies), worst relative "',
              '                "imbalance %.3e\\n", trials, NW, worst);',
              "    if (NW == 0)",
              '        std::printf("NOTE: this network has no conserved moiety, so there "',
              '                    "is nothing for this test to check.\\n");',
              "    if (worst > 1e-12) {",
              '        std::printf("FAIL: %s is not conserved\\n",',
              "                    worstw >= 0 ? WNAME[worstw] : \"?\");",
              "        return 1;",
              "    }",
              '    std::printf("PASS: every conserved moiety is conserved\\n");',
              "    return 0;", "}", ""]
    return "\n".join(lines)


# ===========================================================================
def main():
    ap = argparse.ArgumentParser(
        prog="makeKinetics.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Generate CompLB3D kinetics headers from a reaction file.",
        epilog="See the header of this file for the file format.")
    ap.add_argument("spec", help="the .rxn reaction file")
    ap.add_argument("-o", "--output", default=".", help="folder to write into")
    ap.add_argument("--no-selftest", action="store_true")
    args = ap.parse_args()

    try:
        spec = parse(args.spec)
        biotic = gen_biotic(spec, args.spec)
        abiotic = gen_abiotic(spec, args.spec)
        test = None if args.no_selftest else gen_selftest(spec, args.spec)
    except SpecError as e:
        sys.exit("ERROR in %s\n  %s" % (args.spec, e))

    os.makedirs(args.output, exist_ok=True)
    for name, text in (("defineKinetics.hh", biotic),
                       ("defineAbioticKinetics.hh", abiotic)):
        p = os.path.join(args.output, name)
        with open(p, "w") as f:
            f.write(text)
        print("Wrote %s   (%d lines)" % (p, text.count("\n")))
    if test:
        p = os.path.join(args.output, "selftest_kinetics.cpp")
        with open(p, "w") as f:
            f.write(test)
        print("Wrote %s" % p)
        print("\nProve the chemistry conserves mass before you build anything:")
        print("    g++ -O2 -std=c++11 -o selftest %s && ./selftest" % p)

    print("\nSummary")
    print("  substrates : %s" % " ".join(spec["substrates"]))
    print("  microbes   : %s" % (" ".join(spec["microbes"]) or "none"))
    print("  abiotic    : %d reaction(s)" % len(spec["abiotic"]))
    print("  microbial  : %d organism(s)" % len(spec["microbe"]))
    print("  dissolving : %d phase(s)" % len(spec["dissolve"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
