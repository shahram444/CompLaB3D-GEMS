#!/usr/bin/env python3
"""Build, calibrate and manage the ABC sweep campaign.

    python3 tools/make_campaign.py build    --geometries geometries --out campaign \\
                                            --complab ./complab
    python3 tools/make_campaign.py calibrate --campaign campaign --complab ./complab
    python3 tools/make_campaign.py status   --campaign campaign
    python3 tools/make_campaign.py retry    --campaign campaign

WHAT IT MAKES

    campaign/
      campaign.json          every constant this build used, so a run is reproducible
      runs.csv               one row per case
      calibration.json       the measured timestep per geometry and Peclet
      submit.sbatch          a Slurm ARRAY, one task per case
      run_one.sh             what one task does
      runs/
        run_0000/
          CompLaB.xml        rendered from xml/CompLaB.xml.template
          env.sh             PRT_KABIO, PRT_DT, PRT_MAXFRAC
          params.json        what the collector reads back
          input/geometry.dat
          output/            empty; the run fills it
        run_0001/ ...

THE SWEEP

    30 geometries x 3 Peclet x 3 Damkohler = 270 runs, 21 snapshots each.
    The collector drops the frame at iteration 0, which is identically zero in
    every case because nothing starts in the domain, leaving 5400 training pairs.
    plus the straight duct, built separately by --channel for the analytic check

WHY THERE IS A CALIBRATE STEP
-----------------------------
The transport timestep is not something this script can compute. complab.cpp does

    refNu  = PoreMeanU * characteristic_length / Pe          (line 1276)
    ade_dt = refNu * dx^2 / D[0]                             (line 1296)

and PoreMeanU is the mean speed of the flow field IN THAT GEOMETRY, which only
exists after the flow has been solved. So dt depends on the pore space, and a
campaign that assumed dt = ((tau-0.5)/3) dx^2 / D would be wrong by whatever
factor that geometry's mean speed differs from the nominal one -- typically 5 to
30 per cent, quietly, in the physical duration of every run.

Two things depend on getting it right:

    * the run length. Every case of a given Peclet is meant to end at the SAME
      physical time, so that t_norm = 0.5 in the dataset means the same thing in
      every sample. The number of steps that takes is dt-dependent.
    * PRT_DT, which defineAbioticKinetics.hh uses for its per-step reactant cap.
      A PRT_DT smaller than the real dt makes the cap too loose and lets the
      clamp fire instead; larger, and it throttles a reaction that was fine.

So `calibrate` runs the solver once per (geometry, Peclet) with <ade_max_iT>1,
which solves the flow and prints

    [ADE] dt=... s/iter

and nothing else worth the time, then reads that number back. Thirty short runs
for 270 exact ones. `build --complab PATH` does it for you in one command.

Without a binary, `build` falls back to an estimate from the porosity and says so
in every place a reader might look: on the terminal, in campaign.json, and in
each params.json. check_campaign.py reports an uncalibrated campaign as a warning,
not a failure -- it will still run, the physical durations will just be
approximate and slightly unequal across geometries.

HOW LONG EACH RUN IS
--------------------
    t_end = 1.25 * (L^2 / D) / (Pe + 1)

One crossing of the domain, then a quarter again. The two terms are the diffusive
time L^2/D and the advective time L^2/(Pe D); combining them as a harmonic sum is
what makes the expression right at both ends, advective at high Pe and diffusive
at Pe = 0.3, with no discontinuity between.

It is deliberately CONSERVATIVE. CompLaB's Peclet is built on the mean speed over
the whole box, which is a superficial (Darcy) speed, so the interstitial front
actually arrives sooner than L^2/(Pe D) by roughly the porosity. Running to the
superficial estimate means the front has comfortably crossed in every case rather
than only just.

It is also geometry-INDEPENDENT, which is the point: every case at one Peclet
covers the same physical interval, so time means one thing across the dataset.
"""
import argparse
import csv
import json
import math
import os
import re
import shutil
import stat
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
TEMPLATE = os.path.join(PKG, "xml", "CompLaB.xml.template")

# ---------------------------------------------------------------- constants ---
NX = NY = NZ = 32
DX_UM = 10.0                     # voxel size, micrometres
TAU = 0.8                        # flow relaxation time
D = 1.0e-9                       # m2/s, the same for A, B and C
# THE FEED CONCENTRATION, AND WHY IT IS ONE.
#
# The dataset stores concentration as float16 (that is what
# 3D/tools/dataset_reader.py reads and what every PRT-DeepONet dataset holds),
# and float16's smallest NORMAL value is 6.1e-5. A feed of 1e-3 mol/L would put
# the leading edge of the reaction front -- the part at a thousandth of the feed
# and below -- into the subnormal range, where the relative error is several per
# cent and rising. That edge is exactly where the interesting physics is.
#
# A feed of 1.0 puts five decades of the field inside the normal range, and it
# makes the stored number its own dimensionless concentration: C/A0 = C. The
# chemistry is a model chemistry and Pe and Da are what define the problem, so
# nothing physical is given up by choosing the unit this way.
A0 = 1.0                         # mol/L, the Dirichlet feed of A at x = 0
B0 = 1.0                         # mol/L, the Dirichlet feed of B at x = nx-1
# THE PECLET VALUES, AND WHY THEY STOP AT 2.
#
# B is fed at the OUTLET and has to work upstream against the flow. Its steady
# profile is B0 (e^{Pe x/L} - 1)/(e^{Pe} - 1), so it falls by a factor e over a
# distance L/Pe, and the reaction lives where that tail still overlaps A.
#
# CompLaB's Peclet is built on the mean speed over the whole box, a superficial
# speed, so the speed the solute actually feels is larger by about 1/porosity:
# at porosity 0.32 a nominal Pe of 2 is an interstitial Pe near 6.3, and B's
# reach is about six voxels. That is thin but resolved. A nominal Pe of 30 would
# put it under half a voxel, and the whole reaction would collapse onto the
# outlet plane, where the geometry has nothing left to say about it.
#
# WHY 0.02 RATHER THAN 0. Pe = 0 is a different experiment, not the small end of
# this one: complab.cpp:1156 turns the flow solver off entirely, so there is no
# velocity field, nothing for the surrogate to read as a flow input, and no
# PoreMeanU for the timestep to be derived from (see the calibration note below).
# 0.02 keeps the solver on with advection contributing about two per cent of the
# transport, which is the diffusive limit in every way that matters while leaving
# the machinery identical to the other two cases. Each step up is a factor of ten,
# so the three values are evenly spaced in log Pe, which is the axis the surrogate
# is conditioned on.
#
# WHERE THE TWO FRONTS MEET. Solving A and B as conservative tracers and asking
# where they cross, at the interstitial speed (CompLaB's Pe is superficial, so the
# solute feels about Pe/porosity):
#
#              porosity 0.32   porosity 0.59
#     Pe 0.02      0.508           0.504      symmetric, diffusion sets it
#     Pe 0.2       0.577           0.542
#     Pe 2.0       0.889           0.805      B's reach about 5 voxels
#
# Two things to read off that. The spread down the column is the Peclet effect the
# surrogate is being conditioned on. The spread ACROSS each row is the porosity
# effect at fixed Peclet, and it grows from under a voxel at Pe 0.02 to 2.7 voxels
# at Pe 2.0, so the geometry dependence is strongest at the high end. That is the
# regime the held-out porosities are really testing.
PE_LIST = (0.02, 0.2, 2.0)

# THREE DAMKOHLER VALUES, NOT TWO. campaign01 used 0.1 and 10, a single jump of two
# decades with nothing in between, so a surrogate conditioned on log Da had one
# interval and no way to show whether it had learned a trend or memorised two
# labels. 1.0 is the transition point where reaction and transport times are
# comparable and the front is neither kinetically nor transport limited, which is
# also the least predictable of the three. Evenly spaced in log Da, like Pe.
DA_LIST = (0.1, 1.0, 10.0)
N_SNAP = 21                      # snapshots per run, counting the one at iT = 0
RUN_FACTOR = 1.25                # how far past one crossing each run goes
DELTA_P_SEED = 1.0e-5            # a seed only; Peclet back-solves the real one
NS_MAX = 60000                   # flow iteration cap
MAX_RATE_FRACTION = 0.25         # PRT_MAXFRAC
SPECIES = ["A", "B", "C"]

# THE TWO BOUNDARY SETS.
#
#   mixed      each species has its feed face held and its exit face open. This is
#              the ordinary flow-through arrangement and it is what the sweep runs.
#              B's feed face is the RIGHT one, so B's Neumann goes on the left; a
#              Neumann on B's right would mean B never enters and nothing reacts.
#   dirichlet  every face held: two well-mixed reservoirs with the rock between
#              them. Nothing can be manufactured at a boundary and A-B, A+C and
#              B+C have closed-form steady profiles, which is the only way the
#              duct case can be checked against an answer written down in advance.
#              tools/analytic_channel.py requires this variant.
# WHAT A CAMPAIGN OF EACH ONE ACTUALLY PRODUCED.
#
# The mixed set was run on campaign01, all sixty of its cases, and then measured.
# Nothing in this
# problem can concentrate: no species is fed above 1.0 mol/L and the reaction
# only consumes, so any value above the feed came from a boundary. Over the pore
# voxels of the last frame:
#
#     Pe 0.1    20 of 20 runs clean, worst 1.018 x feed
#     Pe 0.5    16 of 20 clean, worst 1.165
#     Pe 2.0     0 of 20 clean, median 1.46, worst 67.3
#
# The worst case, gid 1 at Pe 2, put B at 67 x its feed ON THE x = 0 PLANE,
# decaying inward over four voxels. That plane is B's Neumann, and it is also
# where the FLOW ENTERS. A zero-gradient condition at an inflow is ill-posed:
# whatever the plane holds is advected into the domain, and the plane then
# re-reads its value from the layer it just fed. At Pe 2 that loop has gain
# above one and it runs away.
#
# At a genuine outflow it is milder and still wrong. In gid 0 at Pe 2, species A
# is Dirichlet at its inlet and Neumann at the outlet, and 52 per cent of its
# pore voxels finished above 1.05 x the feed with a mean of 1.051 -- a passive
# tracer, at Damkohler 0.1, exceeding the only concentration it was ever given.
#
# Hence the default. The dirichlet set holds every plane, so nothing can be
# manufactured at one, and it is the only set with closed-form answers to check
# against. Its zeros are not artificial: "no A in the B reservoir" is true by
# construction.
BOUNDARY_SETS = {
    "mixed": {
        "A_LEFT_T": "Dirichlet", "A_LEFT_V": "A0",
        "A_RIGHT_T": "Neumann",  "A_RIGHT_V": "0.0",
        "B_LEFT_T": "Neumann",   "B_LEFT_V": "0.0",
        "B_RIGHT_T": "Dirichlet", "B_RIGHT_V": "B0",
        "C_LEFT_T": "Neumann",   "C_LEFT_V": "0.0",
        "C_RIGHT_T": "Neumann",  "C_RIGHT_V": "0.0",
    },
    "dirichlet": {
        "A_LEFT_T": "Dirichlet", "A_LEFT_V": "A0",
        "A_RIGHT_T": "Dirichlet", "A_RIGHT_V": "0.0",
        "B_LEFT_T": "Dirichlet", "B_LEFT_V": "0.0",
        "B_RIGHT_T": "Dirichlet", "B_RIGHT_V": "B0",
        "C_LEFT_T": "Dirichlet", "C_LEFT_V": "0.0",
        "C_RIGHT_T": "Dirichlet", "C_RIGHT_V": "0.0",
    },
}
_boundaries = ["dirichlet"]         # set once by build(), read by case_tokens()

L_M = NX * DX_UM * 1e-6          # characteristic length, metres
CHAR_LEN = NX * DX_UM            # in the same unit as <dx>, so micrometres

# The nominal timestep, used only as a fallback and as a sanity bound on the
# measured one: refNu = (tau - 1/2)/3 at Pe = 0.
DT_NOMINAL = ((TAU - 0.5) / 3.0) * (DX_UM * 1e-6) ** 2 / D


def k_from_da(da):
    """Da = k_r L^2 / D with k_r = r(A0,A0)/A0 = k A0, so k = Da D / (A0 L^2).

    The pseudo-first-order collapse is not a convenience: k L^2 / D is
    dimensionless only for a FIRST-order constant, and a second-order k has units
    that do not cancel. A0 is the right reference because nothing starts in the
    domain, so the feed is the only composition the reaction ever sees at full
    strength. defineAbioticKinetics.hh and the patched startup banner use the
    same collapse, so the number printed by the run is this number back again.
    """
    return da * D / (A0 * L_M * L_M)


_run_factor = [RUN_FACTOR]          # set once by build(), read by everything else
_n_snap = [N_SNAP]                  # likewise, so a smoke campaign can ask for fewer


def t_end_seconds(pe):
    return _run_factor[0] * (L_M * L_M / D) / (pe + 1.0)


def dt_estimate(porosity_block):
    """Fallback only, used when no binary was available to measure the real one.

    refNu = PoreMeanU * L / Pe, and the flow was already scaled so that the mean
    OUTLET x-velocity gives the target Peclet. PoreMeanU is the mean |u| over the
    whole box, which exceeds that outlet average because |u| >= |u_x| and because
    a tortuous pore space carries transverse flow. Fifteen per cent is what that
    is worth on packs of this kind; it is a model of the solver rather than the
    solver, which is exactly why `build --complab` does not use it.
    """
    del porosity_block            # the correction is not a function of porosity
    return DT_NOMINAL * 1.15


def steps_for(pe, dt):
    """(ade_max_iT, vtk_interval) giving exactly N_SNAP snapshots.

    The loop writes when iT % V == 0 and runs iT = 0 .. ade_max_iT-1, so
    ade_max_iT = V*(N_SNAP-1) + 1 puts the last write on the last step.
    """
    n = _n_snap[0]
    total = t_end_seconds(pe) / dt
    vtk = max(1, int(round(total / (n - 1))))
    return vtk * (n - 1) + 1, vtk


# ------------------------------------------------------------------ render ---
def render(tokens):
    src = open(TEMPLATE, encoding="utf8").read()
    for k, v in tokens.items():
        src = src.replace("@%s@" % k, str(v))
    left = sorted(set(re.findall(r"@([A-Z0-9_]+)@", src)))
    if left:
        raise SystemExit("template token(s) never filled: %s" % ", ".join(left))
    return src


def case_tokens(gid, run_id, pe, da, dt, porosity, split):
    ade_max, vtk = steps_for(pe, dt)
    bc = dict(BOUNDARY_SETS[_boundaries[0]])
    for k, v in list(bc.items()):
        if v == "A0":
            bc[k] = "%g" % A0
        elif v == "B0":
            bc[k] = "%g" % B0
    tok = {
        "RUN_ID": "run_%04d" % run_id, "GID": gid,
        "POROSITY": "%.4f" % porosity, "SPLIT": split,
        "PE": "%g" % pe, "DA": "%g" % da, "K_ABIO": "%.6g" % k_from_da(da),
        "NX": NX, "NY": NY, "NZ": NZ, "NVOXEL": NX * NY * NZ,
        "DX": "%g" % DX_UM, "CHAR_LEN": "%g" % CHAR_LEN, "L_M": "%.4g" % L_M,
        "TAU": "%g" % TAU, "D": "%g" % D,
        "A0": "%g" % A0, "B0": "%g" % B0,
        "DELTA_P": "%g" % DELTA_P_SEED, "NS_MAX": NS_MAX,
        "ADE_MAX": ade_max, "VTK_INTERVAL": vtk, "N_SNAP": _n_snap[0],
        "BOUNDARIES": _boundaries[0],
    }
    tok.update(bc)
    return tok


ENV_SH = """#!/bin/sh
# Read by defineAbioticKinetics.hh through getenv(), once, on its first call.
# It prints a [KIN] line naming all three back; if that line is missing or shows
# different numbers, the binary was built against the default header and is NOT
# running this campaign's chemistry.
#
# k     = Da * D / (A0 * L^2), L/(mol s).  Da = %(da)g here.
# dt    = the transport timestep this case will actually take, seconds. It has to
#         match the [ADE] dt the run prints, or the per-step cap is calibrated
#         against a step the solver is not taking.
# maxfr = no single step may consume more than this fraction of the scarcer
#         reactant. It keeps the step from overshooting so far that the v1.3.2
#         positivity clamp has to fire; a clamped step is no longer the reaction
#         the case was configured for.
export PRT_KABIO=%(k)r
export PRT_DT=%(dt)r
export PRT_MAXFRAC=%(mf)r
"""

RUN_ONE = """#!/usr/bin/env bash
# One case. Never exits non-zero: it writes status.json either way, so a Slurm
# array task that dies cannot take the collector's triage down with it.
set -u
CASE="${1:?usage: run_one.sh <case directory> [complab binary]}"
BIN="${2:-%(bin)s}"
cd "$CASE" || exit 0
mkdir -p output
. ./env.sh
START=$(date +%%s)
"$BIN" CompLaB.xml > output/run.log 2>&1
RC=$?
END=$(date +%%s)
STATE=ok
REASON=""
if [ $RC -ne 0 ]; then
  STATE=failed; REASON="exit code $RC"
elif ! grep -q '\\[KIN\\] abiotic A+B->C' output/run.log; then
  STATE=failed; REASON="no [KIN] line: the binary is not reading PRT_KABIO"
fi
cat > status.json <<EOF
{"state": "$STATE", "reason": "$REASON", "exit_code": $RC, "wall_s": $((END-START))}
EOF
exit 0
"""

SBATCH = """#!/bin/bash
# ============================================================================
#  ABC sweep  --  Slurm ARRAY submit script, one task per case
#  Written by make_campaign.py. Edit the generator, not this copy.
# ============================================================================
#
#  Submit from the campaign directory, the one holding runs.csv:
#
#      sbatch submit.sbatch
#
#  ONE PROCESS PER CASE, not one big job. Array tasks are separate processes, so
#  a crash in one cannot take the rest down, and run_one.sh never propagates a
#  failure: it writes status.json either way and exits 0. That is what lets the
#  collector triage the campaign afterwards instead of losing it.
#
#      sbatch submit.sbatch                 all %(ncase)d cases
#      sbatch --array=3,7,12 submit.sbatch  three of them again
#
#  make_campaign.py retry prints the --array list for whatever did not finish.
# ============================================================================
#SBATCH --job-name=abc_sweep
#SBATCH --partition=%(partition)s
#SBATCH --array=0-%(last)d%%%(throttle)d      # %(ncase)d cases, at most %(throttle)d at a time
#SBATCH --nodes=1
#SBATCH --ntasks=1                            # SERIAL: see the note below
#SBATCH --cpus-per-task=1
#SBATCH --mem=%(mem)s
#SBATCH --time=%(walltime)s
#SBATCH --output=logs/abc_%%A_%%a.out         # %%A = array job id, %%a = task id
#SBATCH --error=logs/abc_%%A_%%a.err
%(mail)s
# ---------------------------------------------------------------------------
# WHY SERIAL
#
#   Each case is 32x32x32 = 32768 voxels. Split across even 8 ranks that is 4096
#   voxels each, and the ranks spend longer talking than computing. The array
#   gives the parallelism instead: %(throttle)d cases at once, each one process.
#
#   Serial also avoids the interconnect entirely, which matters on Sapelo2: UCX
#   1.12.1 in the foss/2022a stack aborts during connection setup on some nodes
#   with  ib_iface.c:742  Assertion `gid->global.interface_id != 0' failed.
# ---------------------------------------------------------------------------

cd "$SLURM_SUBMIT_DIR" || { echo "cannot cd to $SLURM_SUBMIT_DIR"; exit 1; }
mkdir -p logs

# ---------------------------------------------------------------------------
# MODULES. Load exactly what the binary was BUILT with. Without this the job
# starts in the default environment and dies on a missing libmpi.so before the
# solver has read anything.
# ---------------------------------------------------------------------------
module purge
%(modules)s

CASE=$(awk -F, -v n="$SLURM_ARRAY_TASK_ID" 'NR>1 && $1==n {print $2}' runs.csv)
if [ -z "$CASE" ]; then
    echo "no case for array index $SLURM_ARRAY_TASK_ID in runs.csv"
    exit 0
fi

echo "============================================================"
echo " array job  : $SLURM_ARRAY_JOB_ID  task $SLURM_ARRAY_TASK_ID"
echo " case       : $CASE"
echo " host       : $(hostname)"
echo " directory  : $(pwd)"
echo " started    : $(date)"
echo "============================================================"
module list 2>&1
echo "============================================================"

BIN="%(bin)s"
if [ ! -x "$BIN" ]; then
    echo "ERROR: $BIN is not there or is not executable."
    echo
    echo "This path was baked in when the campaign was built. A relative './complab'"
    echo "means the campaign was built WITHOUT --complab, so the generator had no"
    echo "binary to point at and every task in this array will fail the same way."
    echo "Rebuild with --complab /full/path/to/complab, which also measures the"
    echo "timestep, or put the binary beside this script."
    exit 1
fi

bash ./run_one.sh "runs/$CASE" "$BIN"

echo "============================================================"
echo " finished   : $(date)"
echo " status     : $(cat "runs/$CASE/status.json" 2>/dev/null)"
echo "============================================================"
"""


def write_exec(path, text):
    # newline="\n" ON PURPOSE. These three are shell scripts that will be read by
    # bash, and a campaign built on Windows and copied to a Linux cluster would
    # otherwise arrive with CRLF line endings and die at
    # "bad interpreter: /usr/bin/env bash^M" before running anything.
    with open(path, "w", newline="\n") as f:
        f.write(text)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


# --------------------------------------------------------------- calibrate ---
DT_RE = re.compile(r"\[ADE\]\s+dt=([0-9.eE+-]+)\s*s/iter")


def measure_dt(binary, geom_dat, pe, scratch, porosity=0.5, keep=False):
    """Run the solver for one transport step and read the dt it printed."""
    os.makedirs(os.path.join(scratch, "input"), exist_ok=True)
    os.makedirs(os.path.join(scratch, "output"), exist_ok=True)
    shutil.copyfile(geom_dat, os.path.join(scratch, "input", "geometry.dat"))
    # Any Damkohler will do: this runs ONE step purely to read the [ADE] dt line,
    # and dt comes from the flow field, not from the chemistry. DA_LIST[0] rather
    # than the caller's list because measure_dt is called from calibrate() before
    # any list is in scope here.
    tok = case_tokens(0, 0, pe, DA_LIST[0], DT_NOMINAL, porosity, "calib")
    tok["ADE_MAX"] = 1
    tok["VTK_INTERVAL"] = 1000000          # one write at iT = 0 and no more
    with open(os.path.join(scratch, "CompLaB.xml"), "w", encoding="utf8") as f:
        f.write(render(tok))
    env = dict(os.environ, PRT_KABIO="0", PRT_DT=repr(DT_NOMINAL),
               PRT_MAXFRAC=repr(MAX_RATE_FRACTION))
    log = os.path.join(scratch, "output", "run.log")
    with open(log, "w") as f:
        rc = subprocess.call([os.path.abspath(binary), "CompLaB.xml"],
                             cwd=scratch, env=env, stdout=f, stderr=subprocess.STDOUT)
    text = open(log, encoding="utf8", errors="replace").read()
    m = DT_RE.search(text)
    if not keep:
        shutil.rmtree(scratch, ignore_errors=True)
    if m is None:
        tail = "\n".join(text.strip().splitlines()[-12:])
        raise RuntimeError("no '[ADE] dt=' line (exit %d). Last of the log:\n%s" % (rc, tail))
    return float(m.group(1))


def calibrate(geoms, binary, workdir, pelist=PE_LIST):
    """{"gid,pe": dt} measured, one short run per pair."""
    cal = {}
    print("calibrating the timestep: %d geometries x %d Peclet = %d short runs"
          % (len(geoms), len(pelist), len(geoms) * len(pelist)))
    for g in geoms:
        for pe in pelist:
            key = "%d,%g" % (g["gid"], pe)
            scratch = os.path.join(workdir, "calib_%d_%g" % (g["gid"], pe))
            dt = measure_dt(binary, g["dat"], pe, scratch, g["porosity_block"])
            cal[key] = dt
            ratio = dt / DT_NOMINAL
            flag = "" if 0.5 <= ratio <= 5.0 else "   <-- far from nominal, check the log"
            print("  gid %3d  Pe %-5g  dt = %.6g s  (%.2f x nominal)%s"
                  % (g["gid"], pe, dt, ratio, flag))
    return cal


# ------------------------------------------------------------------- build ---
def read_geometries(gdir, include_channel=False, gids=None):
    idx = os.path.join(gdir, "index.csv")
    if not os.path.isfile(idx):
        raise SystemExit("no index.csv in %s -- run tools/make_geometries.py first" % gdir)
    out = []
    with open(idx, newline="") as f:
        for r in csv.DictReader(f):
            if r["kind"] == "straight_duct" and not include_channel:
                continue
            if r["kind"] != "straight_duct" and include_channel:
                continue
            if gids is not None and int(r["gid"]) not in gids:
                continue
            out.append({"gid": int(r["gid"]), "name": r["name"],
                        "porosity": float(r["porosity"]),
                        "porosity_block": float(r["porosity_block"]),
                        "tortuosity": float(r["tortuosity"]),
                        "split": r["split"],
                        "dat": os.path.join(gdir, r["name"], "geometry.dat")})
    if not out:
        raise SystemExit("no geometries selected from %s%s"
                         % (idx, "" if gids is None else " for --gids %s" % sorted(gids)))
    return out


def build(args):
    _run_factor[0] = float(args.run_factor)
    _n_snap[0] = int(args.snapshots)
    _boundaries[0] = args.boundaries
    if _n_snap[0] < 2:
        raise SystemExit("--snapshots must be at least 2")
    geoms = read_geometries(args.geometries, include_channel=args.channel,
                            gids=set(args.gids) if args.gids else None)
    pelist = tuple(args.pe) if args.pe else PE_LIST
    dalist = tuple(args.da) if args.da else DA_LIST

    out = args.out
    os.makedirs(os.path.join(out, "runs"), exist_ok=True)
    os.makedirs(os.path.join(out, "logs"), exist_ok=True)

    cal, calibrated = {}, False
    calpath = os.path.join(out, "calibration.json")
    if args.calibration and os.path.isfile(args.calibration):
        cal = json.load(open(args.calibration))
        calibrated = True
        print("using the timesteps in %s" % args.calibration)
    elif args.complab:
        cal = calibrate(geoms, args.complab, os.path.join(out, "_calib"), pelist)
        calibrated = True
        shutil.rmtree(os.path.join(out, "_calib"), ignore_errors=True)
    if calibrated:
        json.dump(cal, open(calpath, "w"), indent=2, sort_keys=True)
    else:
        print()
        print("  NOT CALIBRATED. No --complab binary and no --calibration file, so the")
        print("  run lengths below come from an estimate of the solver's timestep rather")
        print("  than from the solver. The campaign will run; the physical duration of")
        print("  each case will be approximate and will differ a little between")
        print("  geometries, which weakens the meaning of t_norm across the dataset.")
        print("  Re-run with --complab ./complab, or run the calibrate subcommand later.")
        print()

    rows, run_id = [], 0
    for g in geoms:
        for pe in pelist:
            dt = cal.get("%d,%g" % (g["gid"], pe), dt_estimate(g["porosity_block"]))
            ade_max, vtk = steps_for(pe, dt)
            for da in dalist:
                name = "run_%04d" % run_id
                cdir = os.path.join(out, "runs", name)
                os.makedirs(os.path.join(cdir, "input"), exist_ok=True)
                os.makedirs(os.path.join(cdir, "output"), exist_ok=True)

                tok = case_tokens(g["gid"], run_id, pe, da, dt, g["porosity"], g["split"])
                with open(os.path.join(cdir, "CompLaB.xml"), "w", encoding="utf8") as f:
                    f.write(render(tok))
                write_exec(os.path.join(cdir, "env.sh"),
                           ENV_SH % {"k": k_from_da(da), "dt": dt,
                                     "mf": MAX_RATE_FRACTION, "da": da})
                shutil.copyfile(g["dat"], os.path.join(cdir, "input", "geometry.dat"))

                params = {
                    # --- what collect_complab_output.py requires, exactly ---
                    "gid": g["gid"], "run_id": run_id,
                    "pe": pe, "da_bio": 0.0, "da_abio": da,
                    "ks_ac_norm": 0.0, "ks_a_norm": 0.0, "y_norm": 0.0,
                    "nx": NX, "ny": NY, "nz": NZ,
                    "species": SPECIES, "microbes": [],
                    # --- everything else this sweep wants recorded ---
                    "case": name, "split": g["split"],
                    "porosity": g["porosity"], "porosity_block": g["porosity_block"],
                    "tortuosity": g["tortuosity"],
                    "k_abio": k_from_da(da), "A0": A0, "B0": B0, "D": D,
                    "dx_um": DX_UM, "tau": TAU,
                    "L_m": L_M, "dt_s": dt, "dt_calibrated": calibrated,
                    "ade_max_iT": ade_max, "vtk_interval": vtk,
                    "snapshots": _n_snap[0],
                    "t_end_s": ade_max * dt, "t_end_target_s": t_end_seconds(pe),
                    "reaction": "A + B -> C,  R = k[A][B]",
                    "boundaries": _boundaries[0],
                }
                json.dump(params, open(os.path.join(cdir, "params.json"), "w"), indent=2)

                rows.append({"index": run_id, "case": name, "gid": g["gid"],
                             "split": g["split"], "porosity": g["porosity"],
                             "pe": pe, "da": da, "k_abio": k_from_da(da),
                             "dt_s": dt, "ade_max_iT": ade_max,
                             "vtk_interval": vtk, "snapshots": _n_snap[0],
                             "t_end_s": ade_max * dt})
                run_id += 1

    with open(os.path.join(out, "runs.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    binpath = os.path.abspath(args.complab) if args.complab else "./complab"
    write_exec(os.path.join(out, "run_one.sh"), RUN_ONE % {"bin": binpath})
    longest = max(r["t_end_s"] for r in rows)
    hours = max(1, int(math.ceil(longest / 3600.0 * 4)))      # generous
    mods = "\n".join("module load %s" % m for m in args.modules)
    mail = ("#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=%s\n"
            "#\n#  One mail for the ARRAY as a whole, not one per task: Slurm only\n"
            "#  mails per task when ARRAY_TASKS is in --mail-type.\n" % args.email
            if args.email else
            "#  No --mail-user was given, so this job mails nobody. Rebuild with\n"
            "#  --email you@uga.edu if you want to be told when it finishes.\n")
    write_exec(os.path.join(out, "submit.sbatch"),
               SBATCH % {"last": len(rows) - 1, "throttle": args.throttle,
                         "ncase": len(rows), "bin": binpath,
                         "partition": args.partition, "mem": args.mem,
                         "modules": mods, "mail": mail,
                         "walltime": args.time or ("%02d:00:00" % min(hours, 48))})

    json.dump({
        "geometries": args.geometries, "n_geometries": len(geoms),
        "pe": list(pelist), "da": list(dalist), "runs": len(rows),
        "nx": NX, "ny": NY, "nz": NZ, "dx_um": DX_UM, "tau": TAU, "D": D,
        "A0": A0, "B0": B0, "L_m": L_M, "snapshots": _n_snap[0],
        "run_factor": _run_factor[0], "max_rate_fraction": MAX_RATE_FRACTION,
        "delta_p_seed": DELTA_P_SEED, "dt_calibrated": calibrated,
        "species": SPECIES, "reaction": "A + B -> C,  R = k[A][B]",
        "boundaries": _boundaries[0],
        "samples": len(rows) * (_n_snap[0] - 1),
        "snapshots_written": len(rows) * _n_snap[0],
    }, open(os.path.join(out, "campaign.json"), "w"), indent=2)

    print("%d cases in %s/runs   (%d geometries x %d Pe x %d Da), boundaries: %s"
          % (len(rows), out, len(geoms), len(pelist), len(dalist), _boundaries[0]))
    print("%d snapshots each, of which the one at iteration 0 is identically zero\n"
          "and is dropped by the collector: %d usable samples"
          % (_n_snap[0], len(rows) * (_n_snap[0] - 1)))
    print()
    print("  %-6s %-22s %-16s %-12s %s"
          % ("Pe", "k  L/(mol s)", "steps", "mean dt", "physical end"))
    for pe in pelist:
        rs = [r for r in rows if r["pe"] == pe]
        print("  %-6g %-22s %-16s %-12s %.4g s"
              % (pe, "%.4g to %.4g" % (min(r["k_abio"] for r in rs),
                                       max(r["k_abio"] for r in rs)),
                 "%d to %d" % (min(r["ade_max_iT"] for r in rs),
                               max(r["ade_max_iT"] for r in rs)),
                 "%.4g s" % (sum(r["dt_s"] for r in rs) / len(rs)),
                 t_end_seconds(pe)))
    print()
    print("next:  python3 tools/check_campaign.py check --campaign %s --complab %s"
          % (out, binpath))
    return 0


# ------------------------------------------------------- calibrate / status ---
def calibrate_cmd(args):
    camp = args.campaign
    cj = json.load(open(os.path.join(camp, "campaign.json")))
    geoms = read_geometries(cj["geometries"], include_channel=False)
    cal = calibrate(geoms, args.complab, os.path.join(camp, "_calib"),
                    tuple(cj["pe"]))
    shutil.rmtree(os.path.join(camp, "_calib"), ignore_errors=True)
    path = os.path.join(camp, "calibration.json")
    json.dump(cal, open(path, "w"), indent=2, sort_keys=True)
    print()
    print("wrote %s" % path)
    print("now rebuild so the cases use it:")
    print("  python3 tools/make_campaign.py build --geometries %s --out %s "
          "--calibration %s" % (cj["geometries"], camp, path))
    return 0


def status(args):
    camp = args.campaign
    rows = list(csv.DictReader(open(os.path.join(camp, "runs.csv"), newline="")))
    done, failed, pending = [], [], []
    for r in rows:
        sp = os.path.join(camp, "runs", r["case"], "status.json")
        if not os.path.isfile(sp):
            pending.append(r)
            continue
        try:
            st = json.load(open(sp))
        except Exception:
            failed.append((r, "unreadable status.json"))
            continue
        if st.get("state") == "ok":
            done.append((r, st))
        else:
            failed.append((r, st.get("reason") or "marked failed"))

    print("%d cases: %d done, %d failed, %d not started"
          % (len(rows), len(done), len(failed), len(pending)))
    if done:
        w = sorted(float(st.get("wall_s", 0)) for _, st in done)
        print("  wall time  median %.0f s   max %.0f s   total %.2f core-hours"
              % (w[len(w) // 2], w[-1], sum(w) / 3600.0))
    for r, why in failed[:20]:
        print("  FAILED %s  gid=%s Pe=%s Da=%s  %s"
              % (r["case"], r["gid"], r["pe"], r["da"], why))
    if len(failed) > 20:
        print("  ... and %d more" % (len(failed) - 20))
    return 1 if failed else 0


def retry(args):
    camp = args.campaign
    rows = list(csv.DictReader(open(os.path.join(camp, "runs.csv"), newline="")))
    again = []
    for r in rows:
        sp = os.path.join(camp, "runs", r["case"], "status.json")
        ok = False
        if os.path.isfile(sp):
            try:
                ok = json.load(open(sp)).get("state") == "ok"
            except Exception:
                ok = False
        if not ok:
            again.append(r["index"])
    if not again:
        print("nothing to retry: every case reports ok")
        return 0
    print("%d case(s) to run again:" % len(again))
    print("  sbatch --array=%s submit.sbatch" % ",".join(again))
    print("or, one at a time:")
    for i in again[:5]:
        c = [r for r in rows if r["index"] == i][0]["case"]
        print("  ./run_one.sh runs/%s" % c)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="write the whole campaign")
    b.add_argument("--geometries", default="geometries")
    b.add_argument("--out", default="campaign")
    b.add_argument("--complab", default=None,
                   help="path to the built solver; given, the timestep is MEASURED")
    b.add_argument("--calibration", default=None, help="reuse a calibration.json")
    b.add_argument("--channel", action="store_true",
                   help="build the straight duct instead of the sphere packs")
    b.add_argument("--pe", type=float, nargs="+", default=None)
    b.add_argument("--da", type=float, nargs="+", default=None)
    b.add_argument("--boundaries", choices=sorted(BOUNDARY_SETS), default="dirichlet",
                   help="'dirichlet' (the default) holds every face: two reservoirs "
                        "with the rock between them. 'mixed' leaves each species' exit "
                        "face open, which reads as the conventional choice and was "
                        "MEASURED TO FAIL on this problem: see the note in the source. "
                        "Use it only knowing that.")
    b.add_argument("--gids", type=int, nargs="+", default=None,
                   help="only these geometry ids. Use it for a smoke campaign: "
                        "--gids 0 --pe 0.02 --da 1.0 is one geometry and one case.")
    b.add_argument("--snapshots", type=int, default=N_SNAP,
                   help="frames per run, counting the one at iteration 0. The "
                        "sweep uses %d; a smoke run wants 3 or 4." % N_SNAP)
    b.add_argument("--run-factor", type=float, default=RUN_FACTOR,
                   help="how far past one domain crossing each case runs. 1.25 is "
                        "the sweep default. The analytic duct case needs steady "
                        "state, so build it with something like 6.")
    b.add_argument("--throttle", type=int, default=12,
                   help="max concurrent Slurm array tasks. 12 puts the 270-case "
                        "sweep at about 2.5 hours of wall clock for roughly 31 "
                        "core-hours of work. Raise it if the queue is generous, "
                        "lower it to be a quieter neighbour.")
    b.add_argument("--partition", default="batch", help="Slurm partition")
    b.add_argument("--mem", default="16gb", help="memory per array task")
    b.add_argument("--time", default=None,
                   help="wall clock per array task, HH:MM:SS. Default is computed "
                        "from the longest run with a generous margin.")
    b.add_argument("--email", default=None,
                   help="address for END and FAIL mail. One message for the array "
                        "as a whole, not one per case.")
    b.add_argument("--modules", nargs="+", default=["foss/2022a"],
                   help="modules the batch job loads. These must be what the "
                        "binary was BUILT with, or it dies on a missing library.")
    b.set_defaults(func=build)

    c = sub.add_parser("calibrate", help="measure the timestep of an existing campaign")
    c.add_argument("--campaign", default="campaign")
    c.add_argument("--complab", required=True)
    c.set_defaults(func=calibrate_cmd)

    s = sub.add_parser("status", help="what ran, what failed")
    s.add_argument("--campaign", default="campaign")
    s.set_defaults(func=status)

    r = sub.add_parser("retry", help="the --array list for whatever did not finish")
    r.add_argument("--campaign", default="campaign")
    r.set_defaults(func=retry)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
