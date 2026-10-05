#!/usr/bin/env python3
"""Check an ABC campaign before you spend a night on it, and triage it after.

    python3 tools/check_campaign.py check --campaign campaign --complab ./complab
    python3 tools/check_campaign.py clamp --campaign campaign

`check` is the gate. It reads every case, and where a binary is given it runs one
of them for a single step under a sentinel rate constant. Fatal findings stop the
submission; warnings do not.

THE CHECK THAT MATTERS MOST
---------------------------
    Does the binary actually read PRT_KABIO?

env.sh sets it in every case, but a solver built against
config/kinetics/defineAbioticKinetics.default.hh has its rate compiled in and
ignores the environment completely. Every case then runs IDENTICAL chemistry
while its params.json records a different Damkohler. Every run succeeds. The
dataset builds. It teaches the network that Damkohler does nothing, and nothing
anywhere says so.

So this runs the solver once with

    PRT_KABIO = 1234.5678

and looks for that number in the [KIN] line of its own log. A build that does not
echo it is the wrong build, and that is fatal.

The same short run is used for two more end-to-end checks, which cost nothing
once it has been paid for:

    the Damkohler banner   Da_r should come back as the Da the case asked for.
                           It will not unless patch/apply_damkohler_feed_reference.py
                           has been applied, because stock v1.3.2 builds its
                           reference composition from <initial_concentration> and
                           every one of those is zero here.
    the achieved Peclet    [NS] Pe achieved should match <Peclet>. A large gap
                           means the flow solver did not converge or the
                           permeability estimate was poor.

EVERYTHING ELSE IT CHECKS
-------------------------
    geometry.dat length against <nx><ny><nz>, which is the failure nothing else
        catches: readGeometry() stops when it has what it wants and ignores the
        rest, so a wrong-sized file gives a truncated pore space that still runs
    percolation along x, and sealed y and z faces
    no code-0 (no dynamics) voxel face-adjacent to pore
    XML well-formedness
    <abiotic_rate_scale> pinned at 1.0 while PRT_KABIO varies
    PRT_DT against the timestep params.json recorded
    ade_max_iT = vtk_interval * (snapshots - 1) + 1, so the last write lands on
        the last step and every run holds the same number of snapshots
    <delta_P> above zero, or the solver silently sets Pe = 0 and the case becomes
        pure diffusion
    cell Peclet, Pe * dx / L, below 2
    the reaction against the per-step cap: k A0 dt / A0 versus PRT_MAXFRAC
    disk footprint
    at least one geometry in each of train, validation and test

`clamp` is for afterwards. It reads the [CLAMP] line out of every run.log and
ranks the runs by how much reaction had to be refused. Above about one per cent,
the run is not the reaction it was configured for. This matters because the
collector silently drops any run more than 5 per cent negative into failures.csv:
without this step the first sign of trouble is a dataset quietly smaller than the
campaign that produced it.
"""
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import xml.dom.minidom as minidom

import numpy as np

# THE SENTINEL RATE CONSTANT.
#
# Distinctive enough to be unmistakable in the [KIN] line, and small enough NOT
# to trip the rate law's per-step cap. The first version of this used 1234.5678,
# which is four orders of magnitude above anything the sweep runs: the cap fired,
# the banner correctly reported the CAPPED rate, and this tool then complained
# that the banner disagreed with an expectation that had ignored the cap. The
# check below now accounts for the cap as well, so either value works; this one
# keeps the arithmetic in the range the campaign actually uses.
SENTINEL_K = 0.12345678
SOLID, WALL, PORE = 0, 1, 2


class Report(object):
    def __init__(self):
        self.fatal, self.warn, self.note = [], [], []

    def F(self, s): self.fatal.append(s)

    def W(self, s): self.warn.append(s)

    def N(self, s): self.note.append(s)

    def show(self):
        for s in self.note:
            print("      %s" % s)
        for s in self.warn:
            print("  WARN  %s" % s)
        for s in self.fatal:
            print("  FATAL %s" % s)
        print()
        if self.fatal:
            print("%d fatal, %d warning. DO NOT SUBMIT until the fatal ones are fixed."
                  % (len(self.fatal), len(self.warn)))
            return 1
        if self.warn:
            print("no fatal findings, %d warning(s). Read them, then submit."
                  % len(self.warn))
            return 0
        print("no findings. The campaign is ready to submit.")
        return 0


def xml_text(dom, path):
    """Value of a nested tag, by tag-name path, or None."""
    node = dom
    for name in path:
        kids = [n for n in node.childNodes
                if n.nodeType == n.ELEMENT_NODE and n.tagName == name]
        if not kids:
            return None
        node = kids[0]
    return "".join(c.data for c in node.childNodes
                   if c.nodeType == c.TEXT_NODE).strip()


def load_dat(path, nx, ny, nz):
    v = np.loadtxt(path, dtype=np.int32)
    if v.size != nx * ny * nz:
        return None, v.size
    return v.reshape((nx, ny, nz)), v.size


def percolates(a):
    pore = (a == PORE)
    nx, ny, nz = a.shape
    seen = np.zeros_like(pore)
    stack = [(0, y, z) for y in range(ny) for z in range(nz) if pore[0, y, z]]
    for p in stack:
        seen[p] = True
    while stack:
        x, y, z = stack.pop()
        for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                           (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            u, v, w = x + dx, y + dy, z + dz
            if 0 <= u < nx and 0 <= v < ny and 0 <= w < nz \
                    and pore[u, v, w] and not seen[u, v, w]:
                seen[u, v, w] = True
                stack.append((u, v, w))
    return bool(seen[nx - 1].any())


def bare_solid_touching_pore(a):
    pore = (a == PORE)
    pad = np.pad(pore, 1, mode="constant", constant_values=False)
    touch = np.zeros_like(pore)
    for ax in range(3):
        for s in (-1, 1):
            sl = [slice(1, -1)] * 3
            sl[ax] = slice(1 + s, pad.shape[ax] - 1 + s)
            touch |= pad[tuple(sl)]
    return int(((a == SOLID) & touch).sum())


def env_values(path):
    out = {}
    for line in open(path, encoding="utf8"):
        m = re.match(r"\s*export\s+(\w+)=([^\s#]+)", line)
        if m:
            out[m.group(1)] = m.group(2)
    return out


# ------------------------------------------------------------------- check ---
def check(args):
    R = Report()
    camp = args.campaign
    try:
        cj = json.load(open(os.path.join(camp, "campaign.json")))
        rows = list(csv.DictReader(open(os.path.join(camp, "runs.csv"), newline="")))
    except Exception as e:
        print("cannot read the campaign: %s" % e)
        return 1

    print("campaign %s: %d cases, %d geometries, Pe %s, Da %s"
          % (camp, len(rows), cj["n_geometries"], cj["pe"], cj["da"]))
    if not cj.get("dt_calibrated"):
        R.W("the timestep was estimated, not measured. Run "
            "`make_campaign.py calibrate --complab ...` and rebuild, or accept that "
            "each case's physical duration is approximate and unequal between "
            "geometries, which weakens t_norm across the dataset.")

    # --- per case, structural --------------------------------------------
    geom_cache, splits, total_bytes = {}, {}, 0
    for r in rows:
        cdir = os.path.join(camp, "runs", r["case"])
        incomplete = False
        for f in ("CompLaB.xml", "env.sh", "params.json", "input/geometry.dat"):
            if not os.path.isfile(os.path.join(cdir, f)):
                R.F("%s is missing %s" % (r["case"], f))
                incomplete = True
        if incomplete:
            continue

        try:
            dom = minidom.parse(os.path.join(cdir, "CompLaB.xml")).documentElement
        except Exception as e:
            R.F("%s: CompLaB.xml is not well-formed XML (%s)" % (r["case"], e))
            continue

        p = json.load(open(os.path.join(cdir, "params.json")))
        env = env_values(os.path.join(cdir, "env.sh"))
        splits.setdefault(p["split"], set()).add(p["gid"])

        nx = int(xml_text(dom, ["LB_numerics", "domain", "nx"]))
        ny = int(xml_text(dom, ["LB_numerics", "domain", "ny"]))
        nz = int(xml_text(dom, ["LB_numerics", "domain", "nz"]))
        if (nx, ny, nz) != (p["nx"], p["ny"], p["nz"]):
            R.F("%s: XML says %dx%dx%d, params.json says %dx%dx%d"
                % (r["case"], nx, ny, nz, p["nx"], p["ny"], p["nz"]))

        dat = os.path.join(cdir, "input", "geometry.dat")
        # Keyed on the CONTENT, not the size. Nine cases share each geometry, so
        # caching is worth it; keying on size alone made an edited file read back
        # as whichever unedited file happened to be the same length, which is
        # precisely the kind of silent substitution this tool exists to catch.
        key = (hashlib.md5(open(dat, "rb").read()).hexdigest(), nx, ny, nz)
        if key not in geom_cache:
            a, n = load_dat(dat, nx, ny, nz)
            geom_cache[key] = (a, n)
        a, n = geom_cache[key]
        if a is None:
            R.F("%s: geometry.dat holds %d integers, the XML declares %d. "
                "readGeometry() will NOT refuse this; it stops when it has what it "
                "wants and the run gets a truncated pore space."
                % (r["case"], n, nx * ny * nz))
            continue
        if not percolates(a):
            R.F("%s: the pore space does not connect x=0 to x=%d. A sealed pore "
                "space has no pressure solution and the flow solver will spend its "
                "whole budget failing to find one." % (r["case"], nx - 1))
        if (a[:, 0, :] == PORE).any() or (a[:, -1, :] == PORE).any() \
                or (a[:, :, 0] == PORE).any() or (a[:, :, -1] == PORE).any():
            R.F("%s: open pore on a y or z face. Those four faces get no boundary "
                "condition: the lattices stream off the block there and read back an "
                "envelope nothing updates." % r["case"])
        nb = bare_solid_touching_pore(a)
        if nb:
            R.F("%s: %d no-dynamics (code 0) voxels are face-adjacent to pore. "
                "Those are holes the lattice streams into and never out of; a grain "
                "needs a bounce-back shell." % (r["case"], nb))

        # --- can both reactants actually get in? --------------------------
        #
        # The boundary TYPE is per species and per face, so a mixed set is legal
        # but is NOT what the sweep uses: the sweep holds every face. What is not
        # legal is a species whose feed
        # face is open: A or B would then never enter, the rate k[A][B] would be
        # zero everywhere for the whole run, and the case would finish
        # successfully with an empty field. Nothing else in the chain notices
        # until verify_dataset.py finds that no C was made, days later.
        feeds = {}
        for iS, nm in enumerate(p["species"]):
            sub = ["chemistry", "substrate%d" % iS]
            held = []
            for side in ("left", "right"):
                t = (xml_text(dom, sub + [side + "_boundary_type"]) or "").lower()
                v = float(xml_text(dom, sub + [side + "_boundary_condition"]) or 0)
                if t == "dirichlet" and v > 0:
                    held.append(side)
            feeds[nm] = held
        for nm in ("A", "B"):
            if nm in feeds and not feeds[nm]:
                R.F("%s: %s has no Dirichlet face with a positive value, so it "
                    "never enters the domain. The rate is k[A][B], so the whole "
                    "run would produce nothing." % (r["case"], nm))
        if feeds.get("A") and feeds.get("B") and feeds["A"] == feeds["B"]:
            R.W("%s: A and B are both fed from the %s face. With one diffusivity "
                "they then satisfy the same equation with the same boundary "
                "conditions and are the same field, which collapses the problem to "
                "one species reacting with itself."
                % (r["case"], feeds["A"][0]))

        # --- the two scaling paths ---------------------------------------
        scale = xml_text(dom, ["simulation_mode", "abiotic_rate_scale"])
        if scale is None or abs(float(scale) - 1.0) > 1e-12:
            R.F("%s: <abiotic_rate_scale> is %s, not 1.0, while PRT_KABIO also "
                "varies. The two multiply, so the Damkohler in params.json is not "
                "the Damkohler that would run." % (r["case"], scale))

        # --- timestep agreement ------------------------------------------
        if "PRT_DT" not in env:
            R.F("%s: env.sh does not set PRT_DT" % r["case"])
        elif abs(float(env["PRT_DT"]) - p["dt_s"]) > 1e-12 * max(1.0, p["dt_s"]):
            R.F("%s: PRT_DT is %s but params.json recorded dt = %g"
                % (r["case"], env["PRT_DT"], p["dt_s"]))
        if "PRT_KABIO" not in env:
            R.F("%s: env.sh does not set PRT_KABIO" % r["case"])
        elif abs(float(env["PRT_KABIO"]) - p["k_abio"]) > 1e-9 * p["k_abio"]:
            R.F("%s: PRT_KABIO is %s but params.json recorded k = %g"
                % (r["case"], env["PRT_KABIO"], p["k_abio"]))

        # --- snapshot arithmetic -----------------------------------------
        ade = int(xml_text(dom, ["LB_numerics", "iteration", "ade_max_iT"]))
        vtk = int(xml_text(dom, ["IO", "save_VTK_interval"]))
        want = vtk * (p["snapshots"] - 1) + 1
        if ade != want:
            R.F("%s: ade_max_iT is %d; %d snapshots at interval %d needs %d, or the "
                "last write does not land on the last step and the runs hold "
                "different numbers of frames"
                % (r["case"], ade, p["snapshots"], vtk, want))

        # --- the seed that turns the flow off ----------------------------
        dp = float(xml_text(dom, ["LB_numerics", "delta_P"]))
        if not dp > 0:
            R.F("%s: <delta_P> is %g. complab_functions.hh sets Pe = 0 when the seed "
                "is below threshold, so the flow solver never runs and the case is "
                "pure diffusion, silently." % (r["case"], dp))

        # --- resolution ---------------------------------------------------
        pe_cell = p["pe"] * p["dx_um"] / (p["L_m"] * 1e6)
        if pe_cell > 2.0:
            R.W("%s: cell Peclet %.2f is above 2. The advection term will ring."
                % (r["case"], pe_cell))

        # --- will the cap fire? -------------------------------------------
        # R = k A0 B0 at the feed; the fraction of A0 consumed in one step is
        # k B0 dt, and the cap trips at PRT_MAXFRAC.
        frac = p["k_abio"] * p["B0"] * p["dt_s"]
        mf = float(env.get("PRT_MAXFRAC", 0.25))
        if frac > mf:
            R.F("%s: one step would consume %.3f of the feed concentration, past the "
                "PRT_MAXFRAC cap of %.3f. The rate law would throttle itself and the "
                "run would not be the Damkohler it claims. Shorten the step."
                % (r["case"], frac, mf))
        elif frac > 0.1 * mf:
            R.W("%s: one step consumes %.3f of the feed, within a factor of ten of "
                "the %.3f cap." % (r["case"], frac, mf))

        # WHAT A RUN ACTUALLY LEAVES ON DISK. The first version of this estimate
        # was four times too small, which is precisely the wrong direction for a
        # number somebody sizes a scratch quota from. Counted properly:
        #
        #   * Palabos writes Density with writeData<T> and T is double, so eight
        #     bytes a voxel, not four, and VtkImageOutput3D base64-encodes it,
        #     which adds a third again.
        #   * there are snapshots+1 frames, not snapshots: the loop writes at
        #     iT = 0, V, ... 20V, and then PHASE 7 writes one more at
        #     iT = ade_max_iT. Every run gets the same extra frame, so the
        #     collector is unaffected.
        #   * the rate fields double the concentration volume. They come on
        #     automatically whenever abiotic kinetics is enabled; there is no
        #     XML switch for them.
        #   * the flow file holds velocityNorm as float plus a 3-component
        #     velocity, so 16 bytes a voxel.
        #   * PHASE 7 writes .chk checkpoints unconditionally, whatever
        #     <save_CHK_interval> says. Those are on the GHOST grid, nx+2 wide,
        #     and D3Q7 for each solute plus D3Q19 for the flow.
        nvox = nx * ny * nz
        nghost = (nx + 2) * ny * nz
        frames = p["snapshots"] + 1
        b64 = 1.37
        nsp = len(p["species"])
        total_bytes += frames * nsp * nvox * 8 * b64          # concentrations
        total_bytes += frames * nsp * nvox * 8 * b64          # rate fields
        total_bytes += frames * nvox * 16 * b64               # flow
        total_bytes += nvox * 8 * b64                         # the mask, once
        total_bytes += (nsp + 1) * nghost * 7 * 8             # solute and mask .chk
        total_bytes += nghost * 19 * 8                        # the flow .chk

    # --- campaign-wide ----------------------------------------------------
    # A campaign built with --channel is the analytic reference case: one duct,
    # run on its own and compared against a closed form. It is not training data
    # and has no split, so requiring one would be a finding against a campaign
    # that is doing exactly what it should.
    analytic_only = set(splits) <= {"analytic"}
    if analytic_only:
        R.N("analytic reference campaign: %d geometry, no train/val/test split"
            % sum(len(v) for v in splits.values()))
        R.N("after it runs: tools/analytic_channel.py --run %s/runs/run_0000"
            % camp)
    else:
        for want in ("train", "val", "test"):
            if not splits.get(want):
                R.F("no geometry in the '%s' split. The split is by geometry, and "
                    "a split with none in it cannot be scored." % want)
            else:
                R.N("%-5s %d geometr%s" % (want, len(splits[want]),
                                           "y" if len(splits[want]) == 1 else "ies"))
    overlap = [g for s in splits for g in splits[s]
               if sum(g in splits[t] for t in splits) > 1]
    if overlap:
        R.F("geometries %s appear in more than one split. A network tested on a pore "
            "space it trained on is being scored on memory, not transfer."
            % sorted(set(overlap)))
    R.N("expected output: about %.1f GB across the campaign (%.0f MB per run, "
        "counting .vti, the rate fields and the closing .chk files)"
        % (total_bytes / 1e9, total_bytes / 1e6 / max(len(rows), 1)))
    if total_bytes > 50e9:
        R.W("that is a lot of disk. Check your quota with `quota` or `df -h` on "
            "the filesystem you are writing to before submitting.")

    # --- the run that proves the binary ------------------------------------
    if args.complab:
        rc = probe(args.complab, camp, rows[0], R, keep=args.keep_probe)
        if rc is None:
            R.F("the probe run did not produce a usable log, so nothing below it "
                "could be checked.")
    else:
        R.W("no --complab given, so the one check that actually matters was skipped: "
            "whether the binary reads PRT_KABIO at all. A build against the default "
            "header runs identical chemistry in every case and says nothing.")

    return R.show()


def probe(binary, camp, row, R, keep=False):
    """One case, one transport step, under a sentinel rate constant."""
    src = os.path.join(camp, "runs", row["case"])
    dst = os.path.join(camp, "_probe")
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(src, dst)
    xmlp = os.path.join(dst, "CompLaB.xml")
    text = open(xmlp, encoding="utf8").read()
    text = re.sub(r"<ade_max_iT>\d+</ade_max_iT>", "<ade_max_iT>1</ade_max_iT>", text)
    text = re.sub(r"<save_VTK_interval>\d+</save_VTK_interval>",
                  "<save_VTK_interval>1000000</save_VTK_interval>", text)
    open(xmlp, "w", encoding="utf8").write(text)

    p = json.load(open(os.path.join(dst, "params.json")))
    env = dict(os.environ, PRT_KABIO=repr(SENTINEL_K),
               PRT_DT=repr(p["dt_s"]), PRT_MAXFRAC="0.25")
    logp = os.path.join(dst, "probe.log")
    print("  probing with %s (one step, PRT_KABIO=%g)" % (row["case"], SENTINEL_K))
    with open(logp, "w") as f:
        rc = subprocess.call([os.path.abspath(binary), "CompLaB.xml"],
                             cwd=dst, env=env, stdout=f, stderr=subprocess.STDOUT)
    log = open(logp, encoding="utf8", errors="replace").read()

    # 1. does the binary read the environment?
    m = re.search(r"\[KIN\] abiotic A\+B->C\s+k_abio=([0-9.eE+-]+)", log)
    if m is None:
        R.F("the probe log has no '[KIN] abiotic A+B->C' line. This binary was NOT "
            "built against kinetics/defineAbioticKinetics.hh from this package, so "
            "it does not read PRT_KABIO and every case in the campaign would run the "
            "same chemistry while claiming different Damkohler numbers.")
    elif abs(float(m.group(1)) - SENTINEL_K) > 1e-3 * SENTINEL_K:
        R.F("the probe set PRT_KABIO=%g and the binary reported k_abio=%s. It is not "
            "reading the environment." % (SENTINEL_K, m.group(1)))
    else:
        R.N("the binary reads PRT_KABIO: probe echoed %s" % m.group(1))

    # 2. is the Damkohler patch in?
    if "Cref = feed" in log:
        das = re.findall(r"Da_r = ([0-9.eE+-]+)", log)
        if das:
            # What the banner SHOULD say at the sentinel rate.
            #
            #   R    = k A0 B0                       the rate law at the feed
            #   but  = scarcest * MAXFRAC / dt       if the per-step cap bites
            #   k_r  = R / A0                        the pseudo-first-order collapse
            #   Da   = k_r L^2 / D
            #
            # The cap has to be in here. Leave it out and a sentinel large enough
            # to be unmistakable is also large enough to be capped, and the tool
            # reports a disagreement that is entirely its own.
            A0 = float(p.get("A0", 1.0))
            B0 = float(p.get("B0", A0))
            dt = float(p["dt_s"])
            Dif = float(p["D"])
            L = float(p["L_m"])
            maxfrac = 0.25
            rate = SENTINEL_K * A0 * B0
            cap = min(A0, B0) * maxfrac / dt if dt > 0 else float("inf")
            capped = rate > cap
            if capped:
                rate = cap
            da_expect = (rate / A0) * L * L / Dif
            da_probe = float(das[0])
            if abs(da_probe - da_expect) > 0.02 * max(da_expect, 1e-30):
                R.W("the banner reported Da_r = %.4g where the sentinel rate should "
                    "give %.4g%s. Check <characteristic_length> and the "
                    "diffusivities."
                    % (da_probe, da_expect,
                       " (the per-step cap is active at this rate)" if capped else ""))
            else:
                R.N("the Damkohler banner agrees with the campaign arithmetic "
                    "(Da_r = %.4g at the sentinel rate%s)"
                    % (da_probe, ", capped" if capped else ""))
    elif "Damkohler: not reported" in log or "Cref = C0" in log:
        R.W("the Damkohler banner is building its reference from "
            "<initial_concentration>, which is zero in every case here, so it will "
            "report nothing useful for the whole campaign. Apply "
            "patch/apply_damkohler_feed_reference.py and rebuild.")

    # 3. did the flow solver hit the Peclet it was asked for?
    m = re.search(r"\[NS\] Pe achieved=([0-9.eE+-]+)\s*\(target=([0-9.eE+-]+)\)", log)
    if m:
        got, want = float(m.group(1)), float(m.group(2))
        if want > 0 and abs(got - want) > 0.05 * want:
            R.W("Pe achieved %.4g against a target of %.4g. The flow solve may not "
                "have converged; raise <ns_max_iT1>." % (got, want))
        else:
            R.N("flow solver hit its Peclet target (%.4g vs %.4g)" % (got, want))
    else:
        R.W("the probe log has no '[NS] Pe achieved' line, so the flow solve cannot "
            "be confirmed.")

    if rc != 0:
        R.F("the probe run exited %d. The log is at %s" % (rc, logp))
    if not keep:
        shutil.rmtree(dst, ignore_errors=True)
    else:
        R.N("probe kept at %s" % dst)
    return rc


# ------------------------------------------------------------------- clamp ---
CLAMP_RE = re.compile(
    r"\[CLAMP\] reaction increments applied: (\d+), of which (\d+) \(([\d.]+)%\).*?"
    r"Total refused: ([0-9.eE+-]+), which is ([\d.]+)%", re.S)


def clamp(args):
    camp = args.campaign
    rows = list(csv.DictReader(open(os.path.join(camp, "runs.csv"), newline="")))
    found, missing = [], []
    for r in rows:
        log = os.path.join(camp, "runs", r["case"], "output", "run.log")
        if not os.path.isfile(log):
            missing.append(r["case"])
            continue
        m = CLAMP_RE.search(open(log, encoding="utf8", errors="replace").read())
        if not m:
            missing.append(r["case"])
            continue
        found.append((float(m.group(5)), r, m))

    if not found:
        print("no [CLAMP] line in any run.log. Either nothing has run yet, or the "
              "binary predates v1.3.2 and has no clamp, in which case a high-Damkohler "
              "run can go negative and the collector will drop it silently.")
        return 1

    # by the refused fraction only: the rows after it are a dict and a match
    # object, and a tie on the fraction would otherwise fall through to comparing
    # those. Every clean campaign is one long tie at zero.
    found.sort(key=lambda x: x[0], reverse=True)
    print("%d runs with a [CLAMP] report. The worst %d:"
          % (len(found), min(len(found), 15)))
    print("  %-10s %-6s %-6s %10s  %s" % ("case", "Pe", "Da", "refused %", "verdict"))
    bad = 0
    for frac, r, m in found[:15]:
        if frac > 1.0:
            verdict = "NOT the reaction it was configured for; run it again shorter"
            bad += 1
        elif frac > 0.1:
            verdict = "noticeable, worth a look"
        else:
            verdict = "fine"
        print("  %-10s %-6s %-6s %10.4f  %s"
              % (r["case"], r["pe"], r["da"], frac, verdict))
    bad = sum(1 for f, _, _ in found if f > 1.0)
    if len(found) > 15:
        print("  ... and %d more, none worse than %.4f%%"
              % (len(found) - 15, found[15][0]))
    if missing:
        print()
        print("%d run(s) with no [CLAMP] line: %s%s"
              % (len(missing), ", ".join(missing[:8]),
                 " ..." if len(missing) > 8 else ""))
    print()
    if bad:
        print("%d run(s) refused more than one per cent of the consumption asked for. "
              "Those are not the Damkohler they claim." % bad)
        return 1
    print("every run applied at least 99 per cent of the reaction it asked for.")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check", help="the gate, before submitting")
    c.add_argument("--campaign", default="campaign")
    c.add_argument("--complab", default=None)
    c.add_argument("--keep-probe", action="store_true")
    c.set_defaults(func=check)
    k = sub.add_parser("clamp", help="triage the finished runs")
    k.add_argument("--campaign", default="campaign")
    k.set_defaults(func=clamp)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
