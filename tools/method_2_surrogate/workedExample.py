#!/usr/bin/env python3
"""
/* This file is a part of the CompLaB program.
 *
 * The CompLaB3D software (3-D pore-scale extension) is developed since 2024
 * by the University of Georgia (United States, Meile Lab, Department of
 * Marine Sciences). The original 2-D CompLaB v1.0 was a collaboration of
 * the University of Georgia and Chungnam National University (South Korea).
 *
 * Contact:
 * Shahram Asgari, Christof Meile
 * Department of Marine Sciences (Meile Lab)
 * University of Georgia, Athens, GA 30602, USA
 *
 * CompLaB is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>.
*/

===============================================================================
THE WORKED NUMERICAL EXAMPLE OF APPENDIX E, SECTION E9
===============================================================================

WHAT IT DOES

  Carries two voxels through every equation of the surrogate path, in both
  the single-output (MISO) and the multi-output (MIMO) form, and compares
  both with the linear program the networks were fitted to:

      Step A   equation (E2)   concentration  ->  uptake bound
               equation (E3)   uptake bound   ->  rescaled input
      Step B   equations (E4), (E5)           hidden layer
      Step C   equations (E6), (E7)           output layer, back to units
               equation (E8)                  what a MISO run has to assume
      Step D   equations (E9), (E10), (E11)   source term and mass budget

  Every intermediate is printed, so the trace can be checked by hand.

WHAT YOU NEED FIRST

  Two fitted networks in the run-time .srg format, from the same sweep:

      python3 generateTrainingData.py e_coli_core.xml \
          --exchange EX_glc__D_e --range 2 20 --log \
          --exchange EX_o2_e     --range 2 30 --log \
          --also EX_ac_e --objective Biomass_Ecoli_core \
          --grid 21 -o train_mimo.csv
      python3 trainSurrogate.py train_mimo.csv --layers 4 --restarts 6 \
          --seed 7 --name mimo_demo --srg mimo.srg -o mimo_weights.hh

  and the same two commands with --growth-only and --name miso_demo for the
  single-output network.  cobrapy is needed only for the reference answers.

CHANGING THE EXAMPLE

  MONOD, BIOMASS and POINTS below are the whole of the physical setup: the
  Michaelis-Menten parameters, the biomass and time step, and the voxel
  concentrations.  A different voxel is a two-line edit.

    python3 workedExample.py                    # print the trace
    python3 workedExample.py --json out.json    # and machine-readable
"""
import argparse, json, math, sys

# ---------------------------------------------------------------- the network

def read_srg(path):
    """Read the run-time .srg format written by trainSurrogate.py."""
    net = {"W": [], "B": []}
    key = None
    for raw in open(path):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        tok = line.split()
        if tok[0] in ("W", "B"):
            key = (tok[0], int(tok[1]))
            net[tok[0]].append([])
            continue
        if key and all(_is_num(t) for t in tok) and tok[0] not in ("inputs",):
            net[key[0]][key[1]].append([float(t) for t in tok])
            continue
        head, rest = tok[0], tok[1:]
        if head in ("inputnames", "outputnames"):
            net[head] = rest
        elif head in ("version", "inputs", "outputs", "logoutput"):
            net[head] = int(rest[0])
        elif head == "layers":
            net["layers"] = [int(t) for t in rest]
        else:
            net[head] = [float(t) for t in rest]
    net["W"] = [[row for row in m] for m in net["W"]]
    return net


def _is_num(t):
    try:
        float(t); return True
    except ValueError:
        return False


def tansig(n):
    return 2.0 / (1.0 + math.exp(-2.0 * n)) - 1.0


def forward(net, x, trace=None):
    """Equations (E3) to (E7), with every intermediate recorded."""
    xo, xg = net["xoffset"], net["xgain"]
    a = [(x[j] - xo[j]) * xg[j] - 1.0 for j in range(len(x))]
    if trace is not None:
        trace["a0"] = a[:]
    for L, (Wm, bv) in enumerate(zip(net["W"], net["B"])):
        b = bv[0]
        n = [sum(Wm[i][j] * a[j] for j in range(len(a))) + b[i]
             for i in range(len(Wm))]
        last = (L == len(net["W"]) - 1)
        a = n[:] if last else [tansig(v) for v in n]
        if trace is not None:
            trace.setdefault("n", []).append(n[:])
            trace.setdefault("a", []).append(a[:])
    yo, yg = net["yoffset"], net["ygain"]
    o = [(a[k] + 1.0) / yg[k] + yo[k] for k in range(len(a))]
    if trace is not None:
        trace["o"] = o[:]
    return o


# ------------------------------------------------------------------ the truth

def lp_truth(model_path, objective, bounds):
    import cobra
    m = cobra.io.read_sbml_model(model_path)
    m.objective = objective
    m.reactions.get_by_id("EX_glc__D_e").lower_bound = -bounds[0]
    m.reactions.get_by_id("EX_o2_e").lower_bound = -bounds[1]
    s = m.optimize()
    return {"status": s.status,
            "growth": float(s.objective_value),
            "glc": float(-s.fluxes["EX_glc__D_e"]),
            "o2": float(-s.fluxes["EX_o2_e"]),
            "ac": float(s.fluxes["EX_ac_e"])}


# ------------------------------------------------------- the two worked cases

MONOD = {                       # V_max in mmol gDW-1 h-1, K in mol/L
    "glucose": {"Vmax": 15.0, "K": 5.0e-4},
    "oxygen":  {"Vmax": 30.0, "K": 2.0e-4},
}
BIOMASS = {"B": 2.0e-2, "M_B": 24.6, "dt": 1.0}   # mol/L, g/mol, s

POINTS = [
    {"tag": "P1", "name": "oxygen scarce, carbon plentiful",
     "C": {"glucose": 1.00e-3, "oxygen": 1.00e-4}},
    {"tag": "P2", "name": "carbon scarce, oxygen plentiful",
     "C": {"glucose": 1.25e-4, "oxygen": 1.00e-3}},
]


def monod(name, C):
    p = MONOD[name]
    return p["Vmax"] * C / (p["K"] + C)


def work(point, miso, mimo, model_path, objective):
    out = {"tag": point["tag"], "name": point["name"], "C": point["C"]}

    # ---- Step A ---------------------------------------------------------
    x = [monod("glucose", point["C"]["glucose"]),
         monod("oxygen",  point["C"]["oxygen"])]
    out["x"] = x

    # ---- Steps B and C, both forms --------------------------------------
    tm, tM = {}, {}
    o_miso = forward(miso, x, tm)
    o_mimo = forward(mimo, x, tM)
    out["miso"] = {"trace": tm, "o": o_miso}
    out["mimo"] = {"trace": tM, "o": o_mimo}

    # ---- the linear program the networks were fitted to ------------------
    out["lp"] = lp_truth(model_path, objective, x)

    # ---- what each form hands the transport step ------------------------
    #   MISO: growth from the network, consumption assumed = the bound (E8)
    out["miso"]["fluxes"] = {"growth": o_miso[0], "glc": x[0], "o2": x[1],
                             "ac": 0.0}
    out["mimo"]["fluxes"] = {"growth": o_mimo[0], "glc": o_mimo[1],
                             "o2": o_mimo[2], "ac": o_mimo[3]}

    #   the linear program, expressed the same way, as the reference
    out["lp"]["fluxes"] = {"growth": out["lp"]["growth"], "glc": out["lp"]["glc"],
                           "o2": out["lp"]["o2"], "ac": out["lp"]["ac"]}

    # ---- Step D ---------------------------------------------------------
    B, MB, dt = BIOMASS["B"], BIOMASS["M_B"], BIOMASS["dt"]
    bm = B * MB                                     # g dry weight per litre
    out["BM"] = bm
    for form in ("miso", "mimo", "lp"):
        f = out[form]["fluxes"]
        R = {k: -f[k] * bm / 3600.0 * 1e-3 for k in ("glc", "o2")}
        R["ac"] = +f["ac"] * bm / 3600.0 * 1e-3
        dC = {"glc": max(R["glc"] * dt, -point["C"]["glucose"]),
              "o2":  max(R["o2"] * dt, -point["C"]["oxygen"]),
              "ac":  R["ac"] * dt}
        out[form]["R"] = R
        out[form]["dC"] = dC
        out[form]["dB"] = f["growth"] * B * dt / 3600.0
        out[form]["clamped"] = {k: (R[k] * dt < -point["C"].get(
            {"glc": "glucose", "o2": "oxygen"}.get(k, ""), 1e30))
            for k in ("glc", "o2")}
    # ---- how far each form is from the linear program --------------------
    out["err"] = {}
    for form in ("miso", "mimo"):
        e = {}
        for k in ("growth", "glc", "o2", "ac"):
            t = out["lp"]["fluxes"][k]
            v = out[form]["fluxes"][k]
            e[k] = {"value": v, "truth": t, "abs": v - t,
                    "rel": (v - t) / t if abs(t) > 1e-12 else None}
        out["err"][form] = e
    return out


def fmt(v, n=6):
    return ("%." + str(n) + "g") % v


def report(res, miso, mimo):
    L = []
    A = L.append
    for r in res:
        A("=" * 74)
        A("POINT %s -- %s" % (r["tag"], r["name"]))
        A("=" * 74)
        A("Step A  concentrations -> uptake bounds, equation (E2)")
        for i, s in enumerate(("glucose", "oxygen")):
            p = MONOD[s]
            A("   %-8s C = %-10s Vmax = %-6s K = %-8s ->  x = %s"
              % (s, fmt(r["C"][s]), fmt(p["Vmax"]), fmt(p["K"]), fmt(r["x"][i])))
        A("")
        for form, net in (("miso", miso), ("mimo", mimo)):
            t = r[form]["trace"]
            A("--- %s ---" % form.upper())
            A("Step A  rescale, equation (E3)")
            A("   a(0) = %s" % [fmt(v) for v in t["a0"]])
            A("Step B  hidden layer, equations (E4) and (E5)")
            A("   n    = %s" % [fmt(v) for v in t["n"][0]])
            A("   a(1) = %s" % [fmt(v) for v in t["a"][0]])
            A("Step C  output layer, equation (E6), then (E7)")
            A("   a(out) = %s" % [fmt(v) for v in t["n"][-1]])
            A("   o      = %s" % [fmt(v) for v in t["o"]])
            f = r[form]["fluxes"]
            A("   growth %s 1/h   glc %s   o2 %s   ac %s"
              % (fmt(f["growth"]), fmt(f["glc"]), fmt(f["o2"]), fmt(f["ac"])))
            A("Step D  equations (E9), (E10), (E11)")
            for k in ("glc", "o2", "ac"):
                A("   R_%-4s = %-12s mol/L/s    dC = %s"
                  % (k, fmt(r[form]["R"][k]), fmt(r[form]["dC"][k])))
            A("   dB     = %s mol/L over dt = %s s"
              % (fmt(r[form]["dB"]), BIOMASS["dt"]))
            A("")
        lp = r["lp"]
        A("--- the linear program, for comparison ---")
        A("   growth %s 1/h   glc %s   o2 %s   ac %s"
          % (fmt(lp["growth"]), fmt(lp["glc"]), fmt(lp["o2"]), fmt(lp["ac"])))
        A("")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--miso", default="miso.srg")
    ap.add_argument("--mimo", default="mimo.srg")
    ap.add_argument("--model", default="e_coli_core.xml")
    ap.add_argument("--objective", default="Biomass_Ecoli_core")
    ap.add_argument("--json")
    a = ap.parse_args()

    miso, mimo = read_srg(a.miso), read_srg(a.mimo)
    res = [work(p, miso, mimo, a.model, a.objective) for p in POINTS]
    print(report(res, miso, mimo))
    if a.json:
        json.dump({"points": res, "monod": MONOD, "biomass": BIOMASS,
                   "miso": miso, "mimo": mimo},
                  open(a.json, "w"), indent=1)
        print("wrote %s" % a.json)


if __name__ == "__main__":
    main()
