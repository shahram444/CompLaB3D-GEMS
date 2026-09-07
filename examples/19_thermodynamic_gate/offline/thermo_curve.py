#!/usr/bin/env python3
"""
19_thermodynamic_gate  --  OFFLINE STEP

Reads input/aom.thm and prints the gate it describes, before the solver runs.

  reads   input/aom.thm
  writes  output/thermo_curve.csv   and a table on the terminal

Self-contained: nothing beyond the standard library.

WHY THIS STEP EXISTS.  Every other rate path in CompLaB3D has an offline stage:
a metabolic model to fetch, a network to train, an expression to fit.  This one
has a file of measured energetics, and the equivalent question is whether those
energetics put the threshold anywhere near the concentrations the run will
actually visit.  A dG0 quoted at 25 C instead of 4 C, or a threshold set to 40
kJ/mol by copying a number from an aerobic paper, produces a gate that is shut
everywhere -- which looks in the output exactly like a reaction that does not
happen, and is impossible to tell apart after the fact.

So: read the file, sweep the product concentrations across the range the run
will cover, and look at where F_T falls through 0.5 and where it reaches zero.
If both are outside the range the case visits, fix the file, not the solver.

    python3 offline/thermo_curve.py
"""
import math
import os
import sys

RGAS = 8.314462618e-3            # kJ / mol / K


# -----------------------------------------------------------------------------
# The .thm parser, in the small subset this script needs.  Deliberately a second
# implementation rather than a binding to the C++ one: if the two disagree about
# what the file says, that disagreement is worth seeing here rather than after a
# week of simulation.
# -----------------------------------------------------------------------------
def read_thm(path):
    model = {"T": 298.15, "scale": 1.0, "conc_scale": 1.0, "reactions": []}
    cur = None
    with open(path) as f:
        for raw in f:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            tok = line.split()
            key = tok[0]

            if key == "temperature":
                model["T"] = float(tok[1])
            elif key == "energy_units":
                model["scale"] = 1.0 if tok[1].lower() == "kj" else 1e-3
            elif key == "concentration_scale":
                model["conc_scale"] = float(tok[1])
            elif key == "gas_constant":
                model["R"] = float(tok[1])
            elif key == "reaction":
                cur = {"name": tok[1], "species": [], "nu": [], "dG0": 0.0, "atp": 0.0,
                       "dGatp": 50.0, "chi": 1.0, "cmin": 1e-9, "dGloss": 0.0,
                       "yield": False, "dGana": 0.0, "carbon": 1.0, "reduction": 4.0,
                       "yieldmax": 1.0, "microbe": "all"}
                model["reactions"].append(cur)
            elif cur is None:
                continue
            elif key == "stoich":
                rest = tok[1:]
                for i in range(0, len(rest) - 1, 2):
                    cur["species"].append(rest[i])
                    cur["nu"].append(float(rest[i + 1]))
            elif key == "microbe":
                cur["microbe"] = tok[1]
            elif key == "yield":
                cur["yield"] = tok[1].lower() in ("on", "yes", "true", "1")
            elif key in ("dG0", "atp", "dGatp", "chi", "cmin", "dGloss",
                         "dGana", "carbon", "reduction", "dGdis", "yieldmax"):
                cur[key] = float(tok[1])

    for r in model["reactions"]:
        for k in ("dG0", "dGatp", "dGloss", "dGana"):
            r[k] *= model["scale"]
    return model


def dG(rxn, conc, T, R):
    """dG = dG0 + RT ln Q, with activities floored at cmin."""
    lnQ = 0.0
    for name, nu in zip(rxn["species"], rxn["nu"]):
        a = max(conc.get(name, 0.0), rxn["cmin"])
        lnQ += nu * math.log(a)
    return rxn["dG0"] + R * T * lnQ + rxn["dGloss"]


def F_T(rxn, g, T, R):
    """F_T = max(0, 1 - exp((dG + m dG_ATP)/(chi R T)))"""
    f = g + rxn["atp"] * rxn["dGatp"]
    x = f / (rxn["chi"] * R * T)
    if x >= 0.0:
        return 0.0
    if x < -700.0:
        return 1.0
    return max(0.0, min(1.0, 1.0 - math.exp(x)))


def dissipation(NoC, red):
    """Heijnen & van Dijken (1992), kJ per C-mol biomass."""
    return (200.0 + 18.0 * (6.0 - NoC) ** 1.8
            + math.exp(((3.8 - red) ** 2) ** 0.16 * (3.6 + 0.4 * NoC)))


def yield_of(rxn, g):
    if not rxn["yield"]:
        return 0.0
    dis = rxn.get("dGdis", 0.0) or dissipation(rxn["carbon"], rxn["reduction"])
    den = -(rxn["dGana"] + dis)
    if den >= 0.0:
        return 0.0
    return min(max(g / den, 0.0), rxn["yieldmax"])


# -----------------------------------------------------------------------------
# The sweep.  The case starts at the inlet composition and, inside the
# aggregate, moves towards the depleted-reactant / accumulated-product end.  The
# path below interpolates between the two, which is the path a voxel actually
# takes as the front moves inwards.
# -----------------------------------------------------------------------------
# The rim is the case's own boundary condition. The core is where the interior
# ends up once the reaction has run: methane drawn down by an order of magnitude,
# and one mole of sulfide and bicarbonate for every mole of methane consumed.
INLET = {"CH4": 1.0e-2, "SO4": 2.8e-2, "HS": 1.0e-5, "HCO3": 2.3e-3}
DEEP = {"CH4": 2.0e-3, "SO4": 2.0e-2, "HS": 1.2e-2, "HCO3": 1.45e-2}
STEPS = 21


def main():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    thm = os.path.join(here, "input", "aom.thm")
    if not os.path.exists(thm):
        print("offline: %s is missing." % thm)
        return 1

    model = read_thm(thm)
    R = model.get("R", RGAS)
    T = model["T"]

    print("== 19_thermodynamic_gate offline: the gate this file describes")
    print("   file        %s" % os.path.relpath(thm, here))
    print("   T           %.2f K   RT = %.4f kJ/mol" % (T, R * T))

    for rxn in model["reactions"]:
        eq = "  ".join("%+g %s" % (n, s) for s, n in zip(rxn["species"], rxn["nu"]))
        print()
        print("   reaction    %s   (%s)" % (rxn["name"], rxn["microbe"]))
        print("   %s" % eq)
        print("   dG0         %.4g kJ/mol" % rxn["dG0"])
        print("   threshold   %.4g kJ/mol   (%g ATP x %g kJ/mol)"
              % (rxn["atp"] * rxn["dGatp"], rxn["atp"], rxn["dGatp"]))
        if rxn["yield"]:
            print("   dG_dis      %.4g kJ per C-mol biomass (Heijnen)"
                  % (rxn.get("dGdis", 0.0) or dissipation(rxn["carbon"], rxn["reduction"])))
        print()
        print("   %-6s %10s %10s %10s %10s %10s %8s %9s"
              % ("depth", "CH4", "SO4", "HS", "HCO3", "dG", "F_T", "yield"))

        rows = []
        half = None
        shut = None
        for k in range(STEPS):
            t = k / float(STEPS - 1)
            conc = {}
            for s in INLET:
                # geometric interpolation: concentrations move over orders of
                # magnitude, not linearly
                conc[s] = INLET[s] * (DEEP[s] / INLET[s]) ** t
                conc[s] *= model["conc_scale"]
            g = dG(rxn, conc, T, R)
            f = F_T(rxn, g, T, R)
            y = yield_of(rxn, g)
            rows.append((t, conc, g, f, y))
            if half is None and f < 0.5:
                half = t
            if shut is None and f <= 0.0:
                shut = t
            print("   %-6.2f %10.3e %10.3e %10.3e %10.3e %10.2f %8.4f %9.5f"
                  % (t, conc["CH4"], conc["SO4"], conc["HS"], conc["HCO3"], g, f, y))

        print()
        if shut is None:
            print("   The gate never shuts along this path. Either the reaction really is")
            print("   favourable everywhere the run will visit, or dG0 or the threshold is")
            print("   wrong. Check dG0 is quoted at %.2f K." % T)
        else:
            # The transition is narrower than one row of the table above, so
            # locate it properly rather than reporting the row it fell between.
            def at(t):
                c = {}
                for s2 in INLET:
                    c[s2] = INLET[s2] * (DEEP[s2] / INLET[s2]) ** t * model["conc_scale"]
                return dG(rxn, c, T, R)

            def bisect(target):
                lo, hi = 0.0, 1.0
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    if F_T(rxn, at(mid), T, R) > target:
                        lo = mid
                    else:
                        hi = mid
                return 0.5 * (lo + hi)

            d50 = bisect(0.5)
            d00 = bisect(1e-12)
            print("   F_T = 0.5  at depth %.4f, dG = %+7.2f kJ/mol" % (d50, at(d50)))
            print("   F_T = 0    at depth %.4f, dG = %+7.2f kJ/mol" % (d00, at(d00)))
            print("   transition width %.2f kJ/mol, which is chi x RT x ln 2 = %.2f"
                  % (abs(at(d50) - at(d00)), rxn["chi"] * R * T * math.log(2.0)))
            print()
            print("   THE TRANSITION IS SHARP, AND THAT IS THE EQUATION, NOT A BUG. With")
            print("   chi = %g the whole fall from open to shut happens over about %.0f kJ/mol,"
                  % (rxn["chi"], 2.0 * rxn["chi"] * R * T))
            print("   which is a few percent of dG0. In a real aggregate that shows up as a")
            print("   reaction front rather than a gradual fade. A larger chi widens it.")
            print()
            print("   Both depths are inside the path the run covers, so the gate will act.")

        out = os.path.join(here, "output", "thermo_curve.csv")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, "w") as f:
            f.write("depth,CH4,SO4,HS,HCO3,dG_kJ_per_mol,F_T,yield_Cmol_per_mol\n")
            for t, conc, g, ft, y in rows:
                f.write("%.4f,%.6e,%.6e,%.6e,%.6e,%.6f,%.6f,%.6f\n"
                        % (t, conc["CH4"], conc["SO4"], conc["HS"], conc["HCO3"], g, ft, y))
        print("   wrote       output/thermo_curve.csv")

    return 0


if __name__ == "__main__":
    sys.exit(main())
