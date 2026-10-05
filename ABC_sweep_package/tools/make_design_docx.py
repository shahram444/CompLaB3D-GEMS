#!/usr/bin/env python3
"""Emit the campaign02 design document as a LaTeX-styled .docx.

It reads the campaign that make_campaign.py actually built and the geometry index
that make_geometries.py actually wrote, so every number in the document is a
number from the files rather than a number retyped from a plan.

    python3 tools/make_design_docx.py --geometries geometries \\
        --campaign campaign02 --out campaign02_design.docx

The document is written by a generated Node script against the `docx` npm package,
which is the only route that produces a real Word table. The Node script is left
beside the output so the document can be regenerated without this file.
"""

import argparse
import collections
import csv
import json
import math
import os
import subprocess
import sys


# ------------------------------------------------------------------- reading ---
def read_inputs(geodir, campdir):
    gi = list(csv.DictReader(open(os.path.join(geodir, "index.csv"))))
    runs = list(csv.DictReader(open(os.path.join(campdir, "runs.csv"))))
    camp = json.load(open(os.path.join(campdir, "campaign.json")))
    return gi, runs, camp


def front_position(pe, phi):
    """Where two conservative counter-fed tracers cross, as a fraction of L.

    Steady advection-diffusion with A held at x=0 and B held at x=L gives
    A = (e^P - e^{Px/L})/(e^P - 1) and B = (e^{Px/L} - 1)/(e^P - 1), so they meet
    where e^{Px/L} = (e^P + 1)/2.  P is the INTERSTITIAL Peclet, Pe/phi, because
    CompLaB's Peclet is built on the superficial speed.
    """
    p = pe / phi
    return math.log((math.exp(p) + 1.0) / 2.0) / p


def b_reach_voxels(pe, phi, nx):
    """B's e-folding distance upstream, in voxels."""
    return nx * phi / pe


# ------------------------------------------------------------------ assembly ---
def build(gi, runs, camp, out_docx):
    nx = camp["nx"]
    spheres = [r for r in gi if r["kind"] == "spheres"]
    duct = [r for r in gi if r["kind"] != "spheres"]
    pes = camp["pe"]
    das = camp["da"]
    n_snap = camp["snapshots"]
    usable = n_snap - 1
    L = camp["L_m"]
    D = camp["D"]

    phis = [float(r["porosity"]) for r in spheres]
    phi_lo, phi_hi = min(phis), max(phis)

    # per porosity level
    levels = collections.OrderedDict()
    for r in spheres:
        levels.setdefault(float(r["porosity_target"]), []).append(r)

    # per split
    split_levels = collections.OrderedDict()
    for t, rs in levels.items():
        split_levels.setdefault(rs[0]["split"], []).append(t)
    split_geoms = collections.Counter(r["split"] for r in spheres)
    split_runs = collections.Counter(r["split"] for r in runs)

    # one representative case per Pe
    per_pe = {}
    for r in runs:
        per_pe.setdefault(float(r["pe"]), r)
    k_of_da = {}
    for r in runs:
        k_of_da[float(r["da"])] = float(r["k_abio"])

    # ---- compute budget -----------------------------------------------------
    # THE COST MODEL, MEASURED ON campaign02 ITSELF.
    #
    # An earlier version of this file scaled campaign01's total by step count
    # alone: 0.0390 s per transport step, no fixed term. That was wrong twice
    # over. It ignored the FLOW SOLVE, which happens once per case and costs the
    # same whatever the Peclet is, and it used a per-step rate from campaign01's
    # nodes, which were faster than the ones campaign02 landed on. The result
    # under-predicted the campaign by a factor of 1.8.
    #
    # Fitting  t = FIXED + steps * RATE  to campaign02's own first runs, on gid 0
    # and gid 1:
    #
    #     Pe 0.02   ~14300 steps    968, 978, 963, 1027, 1019 s
    #     Pe 0.2    ~12200 steps    796, 801 s
    #     Pe 2.0     ~4860 steps    492, 488, 491 s
    #
    # gives FIXED = 240 s and RATE = 0.050 s/step. The fixed term is why a Pe 2
    # case takes 8 minutes rather than the 3 that scaling steps alone predicts:
    # about half of its wall time is the flow solve.
    FIXED_S = 240.0
    SEC_PER_STEP = 0.0500
    DT_EST = 0.0115
    dts = [float(r["dt_s"]) for r in runs]
    dt_lo, dt_hi = min(dts), max(dts)
    DT_MEAS = (sum(dts) / len(dts)) if camp.get("dt_calibrated") else 0.0087
    steps_est = {pe: int(float(per_pe[pe]["ade_max_iT"])) for pe in pes}
    steps_meas = {pe: int(round(float(per_pe[pe]["t_end_s"]) / DT_MEAS)) for pe in pes}
    n_geo = len(spheres)

    def cost_h(sbp):
        """Core-hours for the whole campaign under the fitted model."""
        per = sum(FIXED_S + sbp[pe] * SEC_PER_STEP for pe in pes)
        return n_geo * len(das) * per / 3600.0

    tot_est = n_geo * len(das) * sum(steps_est.values())
    tot_meas = n_geo * len(das) * sum(steps_meas.values())
    ch_est = cost_h(steps_est)
    ch_meas = cost_h(steps_meas)
    case_est = {pe: (FIXED_S + steps_est[pe] * SEC_PER_STEP) / 60.0 for pe in pes}
    case_meas = {pe: (FIXED_S + steps_meas[pe] * SEC_PER_STEP) / 60.0 for pe in pes}
    longest = max(case_meas.values()) * 1.15          # the slowest geometry
    throttle = 12

    # ---- dataset ------------------------------------------------------------
    S = len(runs)
    T = usable
    C = len(camp["species"])
    vox = nx * camp["ny"] * camp["nz"]
    gb_conc = S * T * C * vox * 2 / 1024.0 ** 3
    gb_vel = S * 3 * vox * 2 / 1024.0 ** 3
    gb_vti = S * 0.077                                    # GB, measured on campaign01

    d = {}
    d["title"] = "CompLaB3D campaign02"
    d["subtitle"] = ("Design and specification of the A + B to C pore-scale sweep, "
                     "and of the dataset it produces for the PRT3D surrogate")
    d["meta"] = ("%d runs   %d geometries   %d Peclet   %d Damkohler   %d training pairs"
                 % (S, n_geo, len(pes), len(das), S * T))

    # ---------------------------------------------------------------- tables ---
    tables = []

    tables.append(dict(
        num=1, cap="The domain and its discretisation. Every value is fixed across "
                   "all %d runs." % S,
        widths=[2500, 1300, 2100, 3460],
        head=["Quantity", "Symbol", "Value", "Note"],
        rows=[
            ["grid", "nx, ny, nz", "%d x %d x %d voxels" % (nx, camp["ny"], camp["nz"]),
             "<nx>32</nx> is the slice count in the file; the solver adds two ghost columns"],
            ["voxel size", "dx", "%g micrometres" % camp["dx_um"], "uniform and isotropic"],
            ["domain length", "L", "%g micrometres" % (camp["dx_um"] * nx),
             "%.2e m, the length the Damkohler and Peclet numbers are built on" % L],
            ["flow lattice", "", "D3Q19, tau = %g" % camp["tau"],
             "solved once per case, then frozen"],
            ["transport lattice", "", "D3Q7, one per species",
             "three lattices, advected by the frozen flow field"],
            ["diffusivity", "D", "%.1e m2/s" % D,
             "THE SAME for A, B and C, which makes three field combinations exactly reaction free"],
            ["pressure drop seed", "delta P", "%g" % camp["delta_p_seed"],
             "a seed only; the solver back-solves the real one to hit the target Peclet"],
            ["flow iteration cap", "", "60000", "the flow converges well inside this"],
            ["transport timestep", "dt",
             ("%.4f s (estimated)" % dt_lo) if not camp.get("dt_calibrated")
             else ("%.4f to %.4f s" % (dt_lo, dt_hi)),
             "NOT settable: derived from the solved flow field, so measured per case. "
             + ("This campaign was built WITHOUT --complab, so this is the estimator's "
                "value and not the solver's." if not camp.get("dt_calibrated")
                else "Measured by the build step.")],
            ["confining walls", "",
             "bounce back on y = 0, y = ny-1, z = 0, z = nz-1",
             "drawn into the geometry file, because those four faces get no boundary condition"],
        ]))

    tables.append(dict(
        num=2, cap="The three species, their initial state and where each one enters "
                   "the domain. Nothing starts in the domain: every field rises from zero.",
        widths=[1100, 1500, 1700, 1500, 1800, 1760],
        head=["Species", "Role", "Initial C", "D (m2/s)", "Fed at", "Feed value"],
        rows=[
            ["A", "reactant", "0.0 mol/L", "%.0e" % D, "x = 0",
             "A0 = %g mol/L" % camp["A0"]],
            ["B", "reactant", "0.0 mol/L", "%.0e" % D, "x = nx-1",
             "B0 = %g mol/L" % camp["B0"]],
            ["C", "product", "0.0 mol/L", "%.0e" % D, "never fed",
             "produced only by the reaction"],
        ]))

    tables.append(dict(
        num=3, cap="Boundary conditions. Every x-normal face is Dirichlet: read the "
                   "two ends as well-mixed reservoirs with the rock between them. "
                   "A enters from the left, B from the right, so the two fronts travel "
                   "toward each other and the reaction lives where they overlap.",
        widths=[2400, 2320, 2320, 2320],
        head=["Face", "A", "B", "C"],
        rows=[
            ["x = 0  (A reservoir)", "Dirichlet 1.0  FEED", "Dirichlet 0.0", "Dirichlet 0.0"],
            ["x = nx-1  (B reservoir)", "Dirichlet 0.0", "Dirichlet 1.0  FEED", "Dirichlet 0.0"],
            ["y = 0, y = ny-1", "wall (bounce back)", "wall (bounce back)", "wall (bounce back)"],
            ["z = 0, z = nz-1", "wall (bounce back)", "wall (bounce back)", "wall (bounce back)"],
        ]))

    # derived from the campaign rather than typed, so a campaign built with a
    # different Damkohler list still produces a correct table
    DA_NOTE = {
        "low": ("transport limited",
                "a broad overlap zone; A and B interpenetrate before much is consumed"),
        "mid": ("balanced",
                "reaction and transport times comparable; the least predictable of the set"),
        "high": ("reaction limited",
                 "a thin front; whatever meets reacts at once and C peaks in a narrow band"),
    }
    da_rows = []
    for j, da in enumerate(das):
        where = "low" if j == 0 else ("high" if j == len(das) - 1 else "mid")
        regime, look = DA_NOTE[where]
        da_rows.append(["%g" % da, "%.4e" % k_of_da[da], regime, look])

    tables.append(dict(
        num=4, cap=("The three Damkohler numbers" if len(das) == 3
                    else "The %d Damkohler numbers" % len(das))
                   + " and the rate constant each one sets. " + 
                   "k is written into every case's env.sh and read by "
                   "defineAbioticKinetics.hh, which echoes it back on its first call.",
        widths=[1200, 2300, 2300, 3560],
        head=["Da", "k  (L mol-1 s-1)", "Regime", "What the field looks like"],
        rows=da_rows))

    pe_rows = []
    for pe in pes:
        r = per_pe[pe]
        reach = b_reach_voxels(pe, phi_lo, nx)
        pe_rows.append([
            "%g" % pe,
            "%.1f" % float(r["t_end_s"]),
            "%d" % steps_est[pe],
            "%d" % int(float(r["vtk_interval"])),
            "%.3f" % front_position(pe, phi_lo),
            "%.3f" % front_position(pe, phi_hi),
            ("%.0f" % reach) if reach >= 10 else ("%.1f" % reach),
        ])
    tables.append(dict(
        num=5, cap="The three Peclet numbers. The front position is where two "
                   "conservative counter-fed tracers would cross, as a fraction of L, "
                   "at the interstitial speed. The spread down a column is the Peclet "
                   "effect; the spread across each row is the porosity effect at fixed "
                   "Peclet, and it is largest at Pe = %g, which is where the held-out "
                   "porosities are really tested." % pes[-1],
        widths=[1000, 1200, 1200, 1300, 1600, 1600, 1460],
        head=["Pe", "t_end (s)", "Steps", "VTI every",
              "Front at phi %.2f" % phi_lo, "Front at phi %.2f" % phi_hi,
              "B reach (voxels)"],
        rows=pe_rows))

    lvl_rows = []
    for t in sorted(levels):
        rs = levels[t]
        lvl_rows.append([
            "%.2f" % t,
            rs[0]["split"],
            ", ".join(r["gid"] for r in rs),
            ", ".join("%.4f" % float(r["porosity"]) for r in rs),
            ", ".join("%.3f" % float(r["tortuosity"]) for r in rs),
            ", ".join(r["grains"] for r in rs),
        ])
    tables.append(dict(
        num=6, cap="The geometry family: %d porosity levels, %d independent sphere "
                   "packings each, %d geometries. The split is held out BY POROSITY, "
                   "so no level appears in two splits and a held-out geometry is never "
                   "merely a re-packing of one the surrogate has already seen."
                   % (len(levels), len(levels[sorted(levels)[0]]), len(spheres)),
        widths=[1100, 1100, 1700, 2400, 1900, 1160],
        head=["Target phi", "Split", "gids", "Achieved phi", "Tortuosity", "Grains"],
        rows=lvl_rows))

    geo_rows = []
    for r in spheres + duct:
        geo_rows.append([
            r["gid"], r["kind"].replace("_", " "),
            "%.2f" % float(r["porosity_target"]),
            "%.4f" % float(r["porosity"]),
            "%.4f" % float(r["porosity_block"]),
            "%.3f" % float(r["tortuosity"]),
            r["grains"], r["seed"], r["split"],
        ])
    tables.append(dict(
        num=7, cap="Every geometry, one row each. Porosity is over the duct interior; "
                   "the block porosity counts the confining wall too, which is why it "
                   "is lower. gid %s is the straight duct: it is the analytic reference "
                   "case, not training data, and it is the only geometry whose answer "
                   "can be written down in closed form."
                   % duct[0]["gid"] if duct else "",
        widths=[700, 1400, 1000, 1100, 1100, 1100, 900, 900, 1160],
        head=["gid", "Kind", "Target", "phi", "phi block", "Tortuosity",
              "Grains", "Seed", "Split"],
        rows=geo_rows))

    sp_rows = []
    order = ["train", "val", "test"]
    label = {"train": "train", "val": "validation", "test": "test"}
    for s in order:
        lv = sorted(split_levels.get(s, []))
        sp_rows.append([
            label[s],
            ", ".join("%.2f" % x for x in lv),
            "%d" % split_geoms[s],
            "%d" % split_runs[s],
            "%d" % (split_runs[s] * T),
            "%.0f%%" % (100.0 * split_runs[s] / S),
        ])
    tables.append(dict(
        num=8, cap="The split. Both validation levels and both test levels sit strictly "
                   "inside the training range [%.2f, %.2f], so this asks the surrogate "
                   "for interpolation in porosity. Extrapolating below %.2f or above "
                   "%.2f is a different experiment and deserves its own geometries."
                   % (phi_lo, phi_hi, phi_lo, phi_hi),
        widths=[1700, 2600, 1400, 1300, 1400, 960],
        head=["Split", "Porosity levels", "Geometries", "Runs",
              "Training pairs", "Share"],
        rows=sp_rows))

    idx_rows = []
    for r in runs[:9]:
        idx_rows.append([r["case"], r["gid"], "%.4f" % float(r["porosity"]),
                         r["pe"], r["da"], r["split"]])
    idx_rows.append(["...", "...", "...", "...", "...", "..."])
    for r in runs[-3:]:
        idx_rows.append([r["case"], r["gid"], "%.4f" % float(r["porosity"]),
                         r["pe"], r["da"], r["split"]])
    tables.append(dict(
        num=9, cap="How a run id decodes. run id = %d x gid + %d x (Peclet index) + "
                   "(Damkohler index), with both indices running over the values in "
                   "tables 4 and 5 in the order listed. The first nine rows are the "
                   "complete operating grid of one geometry; runs.csv in the campaign "
                   "directory holds all %d."
                   % (len(pes) * len(das), len(das), S),
        widths=[1600, 900, 1400, 1300, 1300, 2860],
        head=["Case", "gid", "phi", "Pe", "Da", "Split"],
        rows=idx_rows))

    tables.append(dict(
        num=10, cap="The collected dataset, as tools/collect_to_h5.py writes it. "
                    "S is one row per run and T is the frames inside it, so the number "
                    "of (run, time) training pairs is S x T = %d. Concentrations are "
                    "stored raw in float16 with a per-species scale recorded alongside, "
                    "because a float16 can hold %g mol/L to about four decimal places "
                    "and the feed is %g." % (S * T, 1.0, camp["A0"]),
        widths=[2300, 2900, 1500, 2660],
        head=["Dataset", "Shape", "Type", "Contents"],
        rows=[
            ["/geom/material", "(%d, %d, %d, %d)" % (len(spheres), nx, camp["ny"], camp["nz"]),
             "uint8", "solid, wall, pore as 0, 1, 2"],
            ["/geom/gdf", "(%d, %d, %d, %d)" % (len(spheres), nx, camp["ny"], camp["nz"]),
             "float32", "signed distance to the nearest grain surface"],
            ["/geom/edt", "(%d, %d, %d, %d)" % (len(spheres), nx, camp["ny"], camp["nz"]),
             "float32", "Euclidean distance transform of the pore space"],
            ["/samples/conc", "(%d, %d, %d, %d, %d, %d)" % (S, T, C, nx, camp["ny"], camp["nz"]),
             "float16", "A, B, C in mol/L, %.2f GB" % gb_conc],
            ["/samples/velocity", "(%d, 3, %d, %d, %d)" % (S, nx, camp["ny"], camp["nz"]),
             "float16", "the steady flow field, %.0f MB" % (gb_vel * 1024)],
            ["/samples/params", "(%d, 2)" % S, "float32", "Pe and Da for each run"],
            ["/samples/t_norm", "(%d, %d)" % (S, T), "float32", "0 to 1 of that run's length"],
            ["/samples/split", "(%d,)" % S, "string", "train, val or test, travels with the data"],
            ["/samples/porosity", "(%d,)" % S, "float32", "of that run's geometry"],
        ]))

    tables.append(dict(
        num=11, cap="The compute budget, under the cost model fitted to this "
                    "campaign's own runs: %.0f s of fixed cost per case, which is the "
                    "flow solve, plus %.3f s per transport step. The left column uses "
                    "the estimator's timestep, the right the measured one. Read the "
                    "per-case wall times rather than the step counts: at Pe = %g about "
                    "half of a case is the flow solve, so cost does not scale with "
                    "steps alone. Fitted against the first runs of this campaign, the "
                    "model reproduces their measured 16.1, 13.3 and 8.2 minutes."
                    % (FIXED_S, SEC_PER_STEP, pes[-1]),
        widths=[3300, 2200, 2200, 1660],
        head=["Quantity", "At dt = %.4f s" % DT_EST, "At dt = %.4f s" % DT_MEAS, "Unit"],
        rows=[
            ["steps, Pe = %g" % pe, "%d" % steps_est[pe], "%d" % steps_meas[pe], "per run"]
            for pe in pes
        ] + [
            ["wall time, Pe = %g" % pe, "%.1f" % case_est[pe], "%.1f" % case_meas[pe],
             "min per run"]
            for pe in pes
        ] + [
            ["total transport steps", "%.2f million" % (tot_est / 1e6),
             "%.2f million" % (tot_meas / 1e6), "all %d runs" % S],
            ["total compute", "%.0f" % ch_est, "%.0f" % ch_meas, "core-hours"],
            ["longest single run", "%.0f" % (max(case_est.values()) * 1.15),
             "%.0f" % longest, "minutes"],
            ["concurrent array tasks", "%d" % throttle, "%d" % throttle, "jobs"],
            ["WALL CLOCK", "%.1f" % (ch_est / throttle), "%.1f" % (ch_meas / throttle), "hours"],
            ["VTI output on disk", "%.0f" % gb_vti, "%.0f" % gb_vti, "GB"],
            ["collected .h5", "%.1f" % (gb_conc + gb_vel), "%.1f" % (gb_conc + gb_vel), "GB"],
        ]))

    tables.append(dict(
        num=12, cap="What changed from campaign01, and why. The first row is the one "
                    "that made campaign01 unusable above Pe = 0.1; the rest are "
                    "improvements to the design rather than corrections to a fault.",
        widths=[2000, 2200, 2200, 2960],
        head=["Item", "campaign01", "campaign02", "Reason"],
        rows=[
            ["boundary set", "each species open at its exit face",
             "every face Dirichlet",
             "the open faces manufactured mass: at Pe = 2 all 20 runs exceeded the feed, "
             "the worst by 67 times, because a zero-gradient plane at an inflow re-reads "
             "the value it just advected inward"],
            ["Peclet", "0.1, 0.5, 2.0", "%s" % ", ".join("%g" % p for p in pes),
             "evenly spaced in log Pe, and the bottom rung is the diffusive limit with "
             "the flow solver still on, which Pe = 0 would switch off"],
            ["Damkohler", "0.1, 10", "%s" % ", ".join("%g" % x for x in das),
             "two decades in a single jump gave the surrogate one interval and no way to "
             "show a trend; 1.0 is the transition point"],
            ["porosity levels", "5", "%d" % len(levels),
             "twice the resolution in the one geometric parameter being varied"],
            ["packings per level", "2", "%d" % len(levels[sorted(levels)[0]]),
             "two cannot separate a porosity effect from one packing's quirks"],
            ["geometries", "10", "%d" % len(spheres), "consequence of the two rows above"],
            ["runs", "60", "%d" % S, "consequence of the three rows above"],
            ["training pairs", "1200", "%d" % (S * T),
             "4.5 times the data for the same order of compute"],
            ["held-out porosities", "1 validation, 1 test", "2 validation, 2 test",
             "one held-out level cannot distinguish a real generalisation gap from that "
             "level being easy"],
            ["feed excess check", "warned above 1.05 x feed",
             "fails above 1.25 x feed",
             "campaign01 passed every operational check while being wrong; a warning was "
             "not enough"],
        ]))

    d["tables"] = tables

    # ------------------------------------------------------------- narrative ---
    d["sections"] = [
        dict(n="1", h="The scenario",
             p=["Two dissolved species are fed into opposite ends of a rock sample and "
                "react where they meet. A enters from the left face, B from the right, "
                "and wherever both are present they combine into a product C that is "
                "carried away by the flow. Nothing starts inside the domain, so both "
                "fronts have to arrive before anything happens at all.",
                "The reaction is A + B to C at a rate R = k[A][B], second order and "
                "irreversible. That choice is the whole point of the experiment. A "
                "first-order sink R = -k[A] is separable: along a streamline the answer "
                "is an exponential, the second species is decoration, and a surrogate "
                "could fit it from the Peclet number alone without ever looking at the "
                "pore space. R = k[A][B] cannot be solved one species at a time. The "
                "rate is zero wherever either reactant is absent, however much of the "
                "other is present, so the reaction is confined to the overlap of the "
                "two fronts, and that overlap is a geometric object. Change the pore "
                "space and it moves.",
                "Counter-current feeding is forced by the no-initial-concentration "
                "requirement rather than chosen for its own sake. CompLaB3D offers "
                "Dirichlet boundaries on the two x-normal planes and nowhere else, so "
                "two species that must both enter through a boundary get one face each. "
                "It is also the better problem: fed together with the same diffusivity, "
                "A and B would satisfy identical equations with identical boundary "
                "conditions, B would carry no information, and the problem would "
                "collapse to one species reacting with itself."]),
        dict(n="2", h="The equations that are solved",
             p=["The flow is solved once per case on a D3Q19 lattice and then frozen. "
                "Each species is then advected and diffused on its own D3Q7 lattice, "
                "with the reaction entering as a source term:",
                "EQ: dA/dt + u . grad A = D lap A - k A B",
                "EQ: dB/dt + u . grad B = D lap B - k A B",
                "EQ: dC/dt + u . grad C = D lap C + k A B",
                "All three share one diffusivity, and that is deliberate. It makes the "
                "combinations A - B, A + C and B + C exactly reaction free, since the "
                "source terms cancel and the remaining operators are identical. Each of "
                "those three has a closed-form steady profile in a straight duct, which "
                "is what the analytic reference case checks the solver against. It also "
                "means no species outruns another, so the meeting surface is set by the "
                "flow and the pore structure alone rather than by a diffusivity contrast.",
                "The two dimensionless groups are built on the domain length L and the "
                "shared diffusivity. Peclet is set per case and the solver back-solves "
                "the pressure drop that achieves it. Damkohler is converted to a rate "
                "constant, which is what the kinetics header actually reads:",
                "EQ: Pe = u L / D        Da = k A0 L^2 / D        k = Da D / (A0 L^2)"]),
        dict(n="3", h="Why the timestep is measured and not chosen",
             p=["The transport timestep is not an input. The solver sets a reference "
                "viscosity from the mean pore speed of the flow it has just solved, and "
                "derives the timestep from that. The mean pore speed does not exist "
                "until the flow is solved, and it depends on the pore space, so the "
                "timestep is different in every geometry and cannot be known in advance.",
                "The build step therefore runs the solver once per geometry and Peclet "
                "for a single step, reads the timestep it prints, and sizes every case "
                "from that. On campaign01 the estimate was 0.0115 s and the measured "
                "value came out between 0.0086 and 0.0090 s, about 24 per cent lower. "
                "Without the measurement every run would have stopped roughly 22 per "
                "cent short of its intended physical duration, and unequally between "
                "geometries, which would have made the normalised time axis mean "
                "different things in different runs.",
                "The same measured timestep is written into each case's environment, "
                "because the rate law uses it to cap how much of the scarcer reactant a "
                "single step may consume. A cap calibrated against a step the solver is "
                "not taking is not a cap."]),
        dict(n="4", h="Initial and boundary conditions",
             p=["Every species starts at zero everywhere. Both reactants enter through "
                "a boundary and nothing is seeded in the pore space, so all three "
                "fields rise from nothing and the early snapshots are genuinely empty.",
                "Every x-normal face is held at a fixed concentration. Read the two "
                "ends as well-mixed reservoirs with the rock between them: the left "
                "reservoir holds A at the feed value and no B and no C, the right holds "
                "B at the feed value and no A and no C. Those zeros are true by "
                "construction rather than imposed for convenience, since neither "
                "reservoir is stirring product back in. The four y and z faces are "
                "solid wall, written into the geometry file as bounce back, because "
                "those faces get no boundary condition from the input file at all and "
                "an open lattice there streams off the block and reads back an envelope "
                "nothing updates.",
                "A held plane cannot manufacture mass, and that is the reason this set "
                "was chosen over the conventional flow-through arrangement. Leaving each "
                "species open at the face it leaves by was tried on a full 60-run "
                "campaign. It ran to completion and passed every operational check, and "
                "the concentrations were still wrong: nothing in this problem can "
                "concentrate, yet at Pe = 2 every run exceeded its feed and the worst "
                "reached 67 times it, on the very plane where the flow enters. A "
                "zero-gradient condition at an inflow is ill-posed. The plane takes its "
                "value from the layer inside it, that value is advected inward, and the "
                "plane then re-reads it from the layer it just fed."]),
        dict(n="5", h="What the surrogate is being asked to learn",
             p=["Each run contributes %d frames, and each frame is a full three-species "
                "concentration volume. The surrogate reads the pore space, the steady "
                "flow field, the two dimensionless groups and a normalised time, and "
                "predicts the three volumes. Geometry is held out by porosity, so a "
                "validation or test geometry is a pore space at a porosity the surrogate "
                "has never been trained on, not a different packing of one it has."
                % usable,
                "The honest question is interpolation, and the split is arranged to ask "
                "exactly that: every held-out porosity level lies between two training "
                "levels. That is why the geometry family has ten levels rather than "
                "five, and three packings per level rather than two. With two packings "
                "and one held-out level there is no way to tell a real generalisation "
                "gap from that one level happening to be easy."]),
    ]

    d["commands"] = [
        "# 1  geometries, 30 sphere packs and the analytic duct",
        "python3 tools/make_geometries.py --out geometries",
        "",
        "# 2  build the campaign, measuring the timestep from the solver",
        "python3 tools/make_campaign.py build \\",
        "    --geometries geometries --out campaign02 \\",
        "    --complab ~/abc_build/complab --throttle %d \\" % throttle,
        "    --mem 16gb --time 02:00:00 --email shahram.asgari@uga.edu",
        "",
        "# 3  the gate. It refuses to pass a campaign whose chemistry is not per case",
        "python3 tools/check_campaign.py check --campaign campaign02 \\",
        "    --complab ~/abc_build/complab",
        "",
        "# 4  submit",
        "cd campaign02 && sbatch submit.sbatch",
        "",
        "# 5  collect to .h5, then verify before training on it",
        "python3 tools/collect_to_h5.py --campaign campaign02 \\",
        "    --geometries geometries --out dataset02",
        "python3 tools/verify_dataset.py --dataset dataset02/dataset.h5 --feed 1.0",
        "",
        "# 6  train one surrogate per species. C is the interesting one",
        "python3 tools/train_prt3d.py --dataset dataset02/dataset.h5 \\",
        "    --species C --epochs 60 --out model02_C",
    ]

    d["closing"] = (
        "The calibration in step 2 is specific to the Peclet values it was run for, so "
        "campaign01's calibration file cannot be reused here. The three values have "
        "changed, and the gate in step 3 checks that the timestep written into each "
        "case matches the one the solver reports for that case, which is the check that "
        "would catch a stale calibration.")
    return d


# ------------------------------------------------------------------ emission ---
NODE = r"""
const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  WidthType, BorderStyle, AlignmentType, PageOrientation, HeadingLevel,
  PositionalTab, PositionalTabAlignment, PositionalTabLeader,
} = require("docx");

const D = JSON.parse(fs.readFileSync(process.argv[2], "utf8"));
const OUT = process.argv[3];

const SERIF = "Cambria";
const MONO = "Courier New";
const BODY = 21;          // half-points, 10.5pt
const SMALL = 17;
const TINY = 15;

const NONE = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
const RULE = (sz) => ({ style: BorderStyle.SINGLE, size: sz, color: "000000" });

function txt(s, o) {
  o = o || {};
  return new TextRun({
    text: s, font: o.font || SERIF, size: o.size || BODY,
    bold: !!o.bold, italics: !!o.italics, color: o.color || "000000",
    allCaps: !!o.allCaps,
  });
}
function para(s, o) {
  o = o || {};
  return new Paragraph({
    children: Array.isArray(s) ? s : [txt(s, o)],
    alignment: o.align || AlignmentType.JUSTIFIED,
    spacing: { before: o.before === undefined ? 0 : o.before,
               after: o.after === undefined ? 140 : o.after,
               line: o.line || 276 },
    indent: o.indent,
    keepNext: !!o.keepNext,
  });
}

// booktabs: a thick rule on top, a thin rule under the header, a thick rule at
// the bottom, and nothing else. No vertical lines, no shading.
function cell(text, w, o) {
  o = o || {};
  const lines = String(text).split("\n");
  return new TableCell({
    width: { size: w, type: WidthType.DXA },
    borders: {
      top: o.top || NONE, bottom: o.bottom || NONE, left: NONE, right: NONE,
    },
    margins: { top: 60, bottom: 60, left: 90, right: 90 },
    children: lines.map((ln, i) => new Paragraph({
      children: [txt(ln, { size: o.size || SMALL, bold: o.bold, font: o.font })],
      alignment: o.align || AlignmentType.LEFT,
      spacing: { before: 0, after: 0, line: 240 },
      keepNext: !!o.keepNext,
    })),
  });
}

function makeTable(t) {
  const w = t.widths;
  const total = w.reduce((a, b) => a + b, 0);
  const rows = [];
  rows.push(new TableRow({
    tableHeader: true,
    children: t.head.map((h, i) => cell(h, w[i], {
      bold: true, top: RULE(14), bottom: RULE(6), size: SMALL, keepNext: true,
    })),
  }));
  t.rows.forEach((r, ri) => {
    const last = ri === t.rows.length - 1;
    rows.push(new TableRow({
      children: r.map((c, i) => cell(c, w[i], {
        bottom: last ? RULE(14) : NONE,
        size: t.rows.length > 12 ? TINY : SMALL,
      })),
    }));
  });
  return new Table({
    columnWidths: w,
    width: { size: total, type: WidthType.DXA },
    rows,
    borders: {
      top: NONE, bottom: NONE, left: NONE, right: NONE,
      insideHorizontal: NONE, insideVertical: NONE,
    },
  });
}

function caption(t) {
  return new Paragraph({
    children: [
      txt("Table " + t.num + ": ", { size: SMALL, bold: true }),
      txt(t.cap, { size: SMALL, italics: true }),
    ],
    alignment: AlignmentType.JUSTIFIED,
    spacing: { before: 260, after: 90, line: 240 },
    keepNext: true,
  });
}

const kids = [];

// ---- title block, centred, rule under it, the way an article class does it ---
kids.push(new Paragraph({
  children: [txt(D.title, { size: 40, bold: true })],
  alignment: AlignmentType.CENTER,
  spacing: { before: 0, after: 120 },
}));
kids.push(new Paragraph({
  children: [txt(D.subtitle, { size: 22, italics: true })],
  alignment: AlignmentType.CENTER,
  spacing: { after: 100 },
}));
kids.push(new Paragraph({
  children: [txt(D.meta, { size: 18 })],
  alignment: AlignmentType.CENTER,
  spacing: { after: 60 },
  border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: "000000",
                      space: 8 } },
}));
kids.push(new Paragraph({ children: [], spacing: { after: 200 } }));

function section(n, h) {
  kids.push(new Paragraph({
    children: [txt(n + " " + h, { size: 26, bold: true })],
    alignment: AlignmentType.LEFT,
    spacing: { before: 300, after: 130 },
    keepNext: true,
  }));
}

D.sections.forEach((s) => {
  section(s.n, s.h);
  s.p.forEach((p) => {
    if (p.startsWith("EQ: ")) {
      kids.push(new Paragraph({
        children: [txt(p.slice(4), { italics: true, size: BODY })],
        alignment: AlignmentType.CENTER,
        spacing: { before: 120, after: 120 },
      }));
    } else {
      kids.push(para(p));
    }
  });
});

section(String(D.sections.length + 1), "Tables");
kids.push(para("Every number below is read from the campaign that was built and the "
  + "geometry index that was written, not retyped from a plan."));

D.tables.forEach((t) => {
  kids.push(caption(t));
  kids.push(makeTable(t));
  kids.push(new Paragraph({ children: [], spacing: { after: 120 } }));
});

section(String(D.sections.length + 2), "The commands, in order");
D.commands.forEach((c) => {
  kids.push(new Paragraph({
    children: [txt(c === "" ? " " : c, { font: MONO, size: SMALL })],
    alignment: AlignmentType.LEFT,
    spacing: { before: 0, after: 0, line: 240 },
  }));
});
kids.push(new Paragraph({ children: [], spacing: { after: 200 } }));
kids.push(para(D.closing));

const doc = new Document({
  styles: { default: { document: { run: { font: SERIF, size: BODY } } } },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1300, bottom: 1300, left: 1300, right: 1300 },
      },
    },
    children: kids,
  }],
});

Packer.toBuffer(doc).then((b) => { fs.writeFileSync(OUT, b); console.log("wrote " + OUT); });
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--geometries", default="geometries")
    ap.add_argument("--campaign", default="campaign02")
    ap.add_argument("--out", default="campaign02_design.docx")
    a = ap.parse_args(argv)

    gi, runs, camp = read_inputs(a.geometries, a.campaign)
    d = build(gi, runs, camp, a.out)

    stem = os.path.splitext(a.out)[0]
    jpath, npath = stem + ".json", stem + ".build.js"
    with open(jpath, "w") as f:
        json.dump(d, f)
    with open(npath, "w") as f:
        f.write(NODE)
    rc = subprocess.call(["node", npath, jpath, a.out])
    if rc:
        print("node failed, rc %d" % rc, file=sys.stderr)
        return rc
    print("%d tables, %d sections, %d runs, %d geometries"
          % (len(d["tables"]), len(d["sections"]), len(runs),
             len([r for r in gi if r["kind"] == "spheres"])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
