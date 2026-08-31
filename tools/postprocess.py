#!/usr/bin/env python3
"""
Turn a CompLB3D output folder into numbers.

WHY THIS EXISTS

CompLB3D writes .vti snapshots and nothing else. No CSV, no time series, no
totals, no porosity, no breakthrough curve. Every quantitative statement about
a run -- did mass balance, did biomass converge, did the pore space close --
therefore had to be extracted by a script the user wrote themselves, which is
exactly the gap this closes.

Point it at an output folder. It reads every snapshot, works out what the
fields are from the filenames and the XML, and writes:

    summary.csv       one row per timestep: total, mean, min, max of every
                      field, plus porosity and pore volume
    profiles.csv      the x-profile of every field at the last step, averaged
                      over y and z
    massbalance.txt   conservation checks, with a verdict
    plots/*.png       the same, drawn, if matplotlib is installed

USAGE

    python3 postprocess.py path/to/run
    python3 postprocess.py path/to/run --conserve "A+C" --conserve "B+C"
    python3 postprocess.py path/to/run --no-plots

The folder is expected to contain CompLaB.xml and output/. Substrate and
microbe names come from the XML, so the columns are called `acetate` rather
than `subsLattice0`.

REQUIREMENTS

    numpy               required
    matplotlib          optional, for the plots
    vtk                 not required; there is a built-in reader

Part of the CompLaB program.  GNU Affero General Public License v3 or later.
"""

import argparse
import os
import re
import sys
import xml.etree.ElementTree as ET

try:
    import numpy as np
except ImportError:
    sys.exit("ERROR: numpy is required.  python3 -m pip install numpy")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Piping to `head` closes stdout early; without this that shows up as a
# BrokenPipeError traceback rather than as nothing at all.
try:
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)
except (ImportError, ValueError):
    pass

from vtireader import read_vti, VtiError            # noqa: E402


# ---------------------------------------------------------------------------
def read_case(folder):
    """Field names and material numbers, from CompLaB.xml if it is there."""
    info = {"subs": [], "bio": [], "pore": [2], "solid": 0, "bounce": 1,
            "seeds": {}, "dx": None}
    xml = os.path.join(folder, "CompLaB.xml")
    if not os.path.exists(xml):
        return info
    root = ET.parse(xml).getroot()

    def txt(p, d=None):
        e = root.find(p)
        return e.text.strip() if e is not None and e.text else d

    ns = int(txt("chemistry/number_of_substrates", "0") or 0)
    for i in range(ns):
        info["subs"].append(txt("chemistry/substrate%d/name_of_substrates" % i,
                                "substrate%d" % i))
    nm = int(txt("microbiology/number_of_microbes", "0") or 0)
    for i in range(nm):
        info["bio"].append(txt("microbiology/microbe%d/name_of_microbes" % i,
                               "microbe%d" % i))

    mats = root.find("LB_numerics/domain/material_numbers")
    if mats is not None:
        for el in mats:
            vals = [int(v) for v in (el.text or "").split()]
            if el.tag == "pore":
                info["pore"] = vals
            elif el.tag == "solid":
                info["solid"] = vals[0]
            elif el.tag == "bounce_back":
                info["bounce"] = vals[0]
            elif el.tag.startswith("microbe"):
                info["seeds"][el.tag] = vals
    d = txt("LB_numerics/domain/dx")
    u = txt("LB_numerics/domain/unit", "um")
    if d:
        info["dx"] = float(d) * {"m": 1.0, "mm": 1e-3, "um": 1e-6}.get(u, 1e-6)
    return info


def scan(outdir):
    """Group the snapshots by field and iteration.

    Filenames look like  subsLattice0_0000500.vti  --  a stem, an optional
    field index, and a zero-padded iteration."""
    pat = re.compile(r"^(?P<stem>[A-Za-z_]+?)(?P<idx>\d*)_(?P<it>\d+)\.vti$")
    found = {}
    for fn in sorted(os.listdir(outdir)):
        m = pat.match(fn)
        if not m:
            continue
        key = (m.group("stem"), int(m.group("idx") or 0))
        found.setdefault(key, {})[int(m.group("it"))] = os.path.join(outdir, fn)
    return found


def label(stem, idx, info):
    if stem.startswith("subs") and idx < len(info["subs"]):
        return info["subs"][idx]
    if stem.startswith("bio") and idx < len(info["bio"]):
        return info["bio"][idx]
    if stem.startswith("mask"):
        return "mask"
    if stem.startswith("ns"):
        return "flow"
    return "%s%d" % (stem, idx)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        prog="postprocess.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Turn a CompLB3D output folder into CSV, checks and plots.")
    ap.add_argument("folder", help="run folder, containing CompLaB.xml and output/")
    ap.add_argument("--output", default=None, help="output subfolder (default: <folder>/output)")
    ap.add_argument("--conserve", action="append", default=[], metavar="EXPR",
                    help="a sum that must stay constant, e.g. \"A+C\". Repeat "
                         "for several. Names are the <name_of_substrates>.")
    ap.add_argument("--tol", type=float, default=1e-6,
                    help="relative tolerance for the conservation check")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    folder = os.path.abspath(args.folder)
    outdir = args.output or os.path.join(folder, "output")
    if not os.path.isdir(outdir):
        sys.exit("ERROR: no output folder at %s.\n"
                 "       Point me at the RUN folder, the one holding "
                 "CompLaB.xml and output/." % outdir)

    info = read_case(folder)
    groups = scan(outdir)
    if not groups:
        sys.exit("ERROR: no .vti files in %s.\n"
                 "       Either the run wrote nothing (check "
                 "<save_VTK_interval>, and that <track_performance> is false, "
                 "because it suppresses all output), or it has not run yet."
                 % outdir)

    iters = sorted(set(it for g in groups.values() for it in g))
    print("Run      : %s" % folder)
    print("Snapshots: %d iteration(s), %d field(s)" % (len(iters), len(groups)))

    # ---- the mask, for porosity and for masking the other fields ----------
    maskkey = next((k for k in groups if k[0].startswith("mask")), None)
    pore_set = set(info["pore"]) | set(
        v for vals in info["seeds"].values() for v in vals)

    rows = []
    last_fields = {}
    dims = None

    for it in iters:
        row = {"iteration": it}
        mask = None
        if maskkey and it in groups[maskkey]:
            try:
                dims, _, _, arr = read_vti(groups[maskkey][it])
                mask = np.rint(list(arr.values())[0]).astype(int)
            except VtiError as e:
                print("  warning: could not read the mask at step %d: %s" % (it, e))

        if mask is not None:
            openv = np.isin(mask, list(pore_set))
            row["porosity"] = float(openv.mean())
            row["pore_voxels"] = int(openv.sum())
        else:
            openv = None

        for (stem, idx), byit in sorted(groups.items()):
            if stem.startswith("mask") or it not in byit:
                continue
            name = label(stem, idx, info)
            try:
                dims, _, _, arr = read_vti(byit[it])
            except VtiError as e:
                print("  warning: %s: %s" % (os.path.basename(byit[it]), e))
                continue
            for aname, a in arr.items():
                a = np.asarray(a, dtype=float)
                if a.ndim == 4:                     # a vector field, e.g. velocity
                    a = np.linalg.norm(a, axis=-1)
                col = name if len(arr) == 1 else "%s_%s" % (name, aname)
                sel = a[openv] if (openv is not None and openv.shape == a.shape) else a.ravel()
                if sel.size == 0:
                    continue
                row["%s_total" % col] = float(sel.sum())
                row["%s_mean" % col] = float(sel.mean())
                row["%s_min" % col] = float(sel.min())
                row["%s_max" % col] = float(sel.max())
                if it == iters[-1]:
                    last_fields[col] = (a, openv)
        rows.append(row)
        if not args.quiet:
            print("  step %-9d %s" % (it, "porosity %.4f" % row["porosity"]
                                      if "porosity" in row else ""))

    # ---- summary.csv ------------------------------------------------------
    cols = ["iteration"] + sorted(c for c in
                                  {k for r in rows for k in r} - {"iteration"})
    spath = os.path.join(folder, "summary.csv")
    with open(spath, "w") as f:
        f.write("# CompLB3D run summary, one row per snapshot\n")
        f.write("# totals and means are over OPEN voxels only where a mask "
                "was available\n")
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join("%.17g" % r[c] if c in r else "" for c in cols) + "\n")
    print("\nWrote %s   (%d rows, %d columns)" % (spath, len(rows), len(cols)))

    # ---- profiles.csv -----------------------------------------------------
    if last_fields:
        ppath = os.path.join(folder, "profiles.csv")
        names = sorted(last_fields)
        nx = list(last_fields.values())[0][0].shape[2]
        with open(ppath, "w") as f:
            f.write("# x-profiles at the last snapshot, averaged over y and z\n")
            f.write("# open voxels only; solid voxels are excluded from the mean\n")
            f.write("x," + ",".join(names) + "\n")
            for x in range(nx):
                vals = []
                for n in names:
                    a, openv = last_fields[n]
                    if openv is not None and openv.shape == a.shape:
                        sl, ol = a[:, :, x], openv[:, :, x]
                        vals.append(float(sl[ol].mean()) if ol.any() else float("nan"))
                    else:
                        vals.append(float(a[:, :, x].mean()))
                f.write("%d," % x + ",".join("%.17g" % v for v in vals) + "\n")
        print("Wrote %s   (%d rows)" % (ppath, nx))

    # ---- mass balance -----------------------------------------------------
    mpath = os.path.join(folder, "massbalance.txt")
    lines = ["CompLB3D mass balance check", "=" * 60, ""]
    verdicts = []

    def series(expr):
        """Total of a sum like 'A+C' at every snapshot."""
        terms = [t.strip() for t in expr.split("+")]
        out = []
        for r in rows:
            s, ok = 0.0, True
            for t in terms:
                k = "%s_total" % t
                if k not in r:
                    ok = False
                    break
                s += r[k]
            out.append(s if ok else None)
        return terms, out

    for expr in args.conserve:
        terms, ser = series(expr)
        if any(v is None for v in ser):
            lines.append("%-24s SKIPPED: no column for one of %s"
                         % (expr, ", ".join(terms)))
            lines.append("%-24s available: %s" % ("", ", ".join(
                sorted(c[:-6] for c in cols if c.endswith("_total")))))
            continue
        first, last = ser[0], ser[-1]
        scale = max(abs(first), abs(last), 1e-300)
        drift = abs(last - first) / scale
        good = drift <= args.tol
        verdicts.append(good)
        lines.append("%-24s first %.10g  last %.10g  relative drift %.3e  %s"
                     % (expr, first, last, drift, "PASS" if good else "FAIL"))

    # always-on sanity checks
    lines.append("")
    lines.append("Sanity checks")
    for c in cols:
        if c.endswith("_min"):
            worst = min(r[c] for r in rows if c in r)
            if worst < -1e-9:
                lines.append("  %-28s went NEGATIVE: %.3e  <- not physical"
                             % (c[:-4], worst))
                verdicts.append(False)
    if "porosity" in cols:
        ps = [r["porosity"] for r in rows if "porosity" in r]
        lines.append("  porosity %.4f -> %.4f  (%s)"
                     % (ps[0], ps[-1],
                        "unchanged" if abs(ps[-1] - ps[0]) < 1e-12 else
                        "fell" if ps[-1] < ps[0] else "rose"))
    if not any(l.startswith("  ") and "NEGATIVE" in l for l in lines):
        lines.append("  no field went negative")

    lines.append("")
    lines.append("VERDICT: " + ("PASS" if all(verdicts) else
                                "FAIL" if verdicts else "no checks requested"))
    with open(mpath, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("Wrote %s" % mpath)
    print()
    print("\n".join(lines[2:]))

    # ---- plots ------------------------------------------------------------
    if not args.no_plots:
        try:
            make_plots(folder, rows, cols, last_fields)
        except ImportError:
            print("\n(matplotlib not installed, so no plots. "
                  "python3 -m pip install matplotlib)")

    return 0 if all(verdicts) else 1


def make_plots(folder, rows, cols, last_fields):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pdir = os.path.join(folder, "plots")
    os.makedirs(pdir, exist_ok=True)
    it = [r["iteration"] for r in rows]

    totals = [c for c in cols if c.endswith("_total")]
    if totals:
        fig, ax = plt.subplots(figsize=(7, 4.2))
        for c in totals:
            ax.plot(it, [r.get(c, float("nan")) for r in rows], label=c[:-6])
        ax.set_xlabel("iteration")
        ax.set_ylabel("total over open voxels")
        ax.legend(fontsize=8, ncol=2)
        ax.set_title("Totals")
        fig.tight_layout()
        fig.savefig(os.path.join(pdir, "totals.png"), dpi=130)
        plt.close(fig)

    if "porosity" in cols:
        fig, ax = plt.subplots(figsize=(7, 3.6))
        ax.plot(it, [r.get("porosity", float("nan")) for r in rows], color="#9C3B2E")
        ax.set_xlabel("iteration")
        ax.set_ylabel("porosity")
        ax.set_title("Porosity")
        fig.tight_layout()
        fig.savefig(os.path.join(pdir, "porosity.png"), dpi=130)
        plt.close(fig)

    if last_fields:
        fig, ax = plt.subplots(figsize=(7, 4.2))
        for n, (a, openv) in sorted(last_fields.items()):
            nx = a.shape[2]
            prof = []
            for x in range(nx):
                if openv is not None and openv.shape == a.shape:
                    sl, ol = a[:, :, x], openv[:, :, x]
                    prof.append(sl[ol].mean() if ol.any() else float("nan"))
                else:
                    prof.append(a[:, :, x].mean())
            ax.plot(range(nx), prof, label=n)
        ax.set_xlabel("x (voxels)")
        ax.set_ylabel("value, averaged over y and z")
        ax.legend(fontsize=8, ncol=2)
        ax.set_title("Profiles at the last snapshot")
        fig.tight_layout()
        fig.savefig(os.path.join(pdir, "profiles.png"), dpi=130)
        plt.close(fig)

    print("Wrote %s/*.png" % pdir)


if __name__ == "__main__":
    sys.exit(main())
