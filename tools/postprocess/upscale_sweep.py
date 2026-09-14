#!/usr/bin/env python3
"""Sweep aggregate size and bulk composition, and turn the runs into an upscaling law.

WHAT THIS IS FOR

A pore-scale run resolves one aggregate and reports one effectiveness factor:

    eta = <r> inside the aggregate / r(C_bulk)

That is one point. A continuum model needs the curve -- eta as a function of the
things that set it -- because it will be asked for a rate at sizes and
compositions this one run did not visit.

Two groups set it here. The first is the Thiele modulus, which is the classical
one and says how far the reaction outruns diffusion across the aggregate:

    phi = R sqrt(k/D),   k = r(C_bulk)/C_bulk

The second is the thermodynamic gate at the bulk composition, F_T(C_bulk). This
one is not in any textbook treatment, and it is why the sweep has two axes rather
than one. With Monod kinetics alone the core of an aggregate runs slowly. With
the gate on it STOPS, at the depth where its own products have raised the free
energy past what the organism can use. That is a moving internal boundary, and a
first-order Thiele analysis does not anticipate it: measured eta falls below the
classical curve, and by more as the bulk sits closer to the threshold.

WHAT IT DOES

For every (radius, bulk composition) pair it writes a geometry and a CompLaB.xml,
runs the solver, and reads the effectiveness factor out of output/upscaling.csv.
It then fits the ratio of measured to classical eta against F_T(C_bulk), which is
the correction a continuum model needs on top of the textbook curve, and writes
the whole table out so the fit can be checked rather than trusted.

    python3 offline/upscale.py --radii 4,6,8 --bulk 0.01,0.02 --iters 6000

Standard library only. Run it from a staged case directory, after ./pipeline.sh
has built the solver once:

    ./scripts/setup_case.sh 19_thermodynamic_gate run/sweep
    cd run/sweep && ./pipeline.sh && python3 offline/upscale.py
"""
import argparse
import csv
import math
import os
import re
import subprocess
import sys

HS_INDEX = 2      # substrate2 is sulfide: the product that shuts the gate


def build_geometry(nx, ny, nz, radius, pore=2, wall=1, bio=3):
    """A sphere of the given radius, centred, in a box closed on four sides.

    x = 0 and x = NX-1 are left open: those two carry the inlet and the outlet.
    The other four faces are wall, because the solver gives them no boundary
    condition of its own -- see src/complab3d_outerfaces.hh.
    """
    cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
    cells = []
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                if y in (0, ny - 1) or z in (0, nz - 1):
                    cells.append(wall)
                    continue
                d2 = (x - cx) ** 2 + (y - cy) ** 2
                cells.append(bio if d2 <= radius * radius else pore)
    return cells


def set_tag(xml, tag, value, occurrence=0):
    """Replace the nth <tag>...</tag> in the document."""
    hits = list(re.finditer(r"<%s>[^<]*</%s>" % (tag, tag), xml))
    if len(hits) <= occurrence:
        return xml
    m = hits[occurrence]
    return xml[:m.start()] + "<%s>%s</%s>" % (tag, value, tag) + xml[m.end():]


def substrate_block(xml, index):
    m = re.search(r"<substrate%d>.*?</substrate%d>" % (index, index), xml, re.S)
    return m


def make_case(base_xml, nx, ny, nz, geom_name, hs_bulk, iters, out_dir, interval):
    """One variant of the case: a geometry, a HELD bulk, a run length.

    EVERY SPECIES IS HELD, NOT JUST THE SUBSTRATES. That is what "the aggregate sits
    in a bulk of composition C" means, and it is the condition every
    effectiveness-factor definition assumes. The shipped case does not do this -- it
    feeds methane from one end and lets everything leave at the other -- so the water
    around the aggregate drifts, the gate at the bulk drifts with it, and eta never
    stops moving: it fell from 1.000 to 0.918 over the first thousand steps and was
    still falling at four diffusion times.

    Holding the products does not impose the answer. The gate closes because sulfide
    accumulates INSIDE the aggregate relative to the bulk, and that gradient is free
    to develop however it likes. What is fixed is only the outside, which is the
    single number a continuum grid block would have.

    Fixing the bulk sulfide is also what makes the second axis of this sweep a
    control rather than an observation: F_T(C_bulk) is set by it, so the sweep can
    ask what happens at a given distance from the energy threshold instead of
    waiting to see where a drifting run ends up.
    """
    x = base_xml
    x = set_tag(x, "nx", nx)
    x = set_tag(x, "ny", ny)
    x = set_tag(x, "nz", nz)
    x = set_tag(x, "filename", geom_name)
    x = set_tag(x, "ade_max_iT", iters)
    x = set_tag(x, "output_path", out_dir)
    x = set_tag(x, "interval", interval)
    x = set_tag(x, "save_VTK_interval", max(iters, 1))    # one volume, at the end
    # Hold the catalyst. eta is a property of the concentration profile, which relaxes
    # on R^2/D; biomass moves on its own much slower timescale and NOT uniformly -- the
    # rim is better fed than the core -- so without this eta keeps drifting after the
    # profile has settled and no run ever reports steady. Every run in this sweep was
    # still moving at five diffusion times until this was switched on.
    x = set_tag(x, "freeze_biomass", "true")

    idx = 0
    while True:
        blk = substrate_block(x, idx)
        if not blk:
            break
        b = blk.group(0)
        if idx == HS_INDEX and hs_bulk is not None:
            val = repr(hs_bulk)
        else:
            m = re.search(r"<initial_concentration>([^<]*)</initial_concentration>", b)
            val = m.group(1).strip() if m else None
        if val is not None:
            nb = set_tag(b, "initial_concentration", val)
            nb = set_tag(nb, "left_boundary_type", "Dirichlet")
            nb = set_tag(nb, "right_boundary_type", "Dirichlet")
            nb = set_tag(nb, "left_boundary_condition", val)
            nb = set_tag(nb, "right_boundary_condition", val)
            x = x[:blk.start()] + nb + x[blk.end():]
        idx += 1
    return x


def diffusion_iterations(xml, radius_voxels, times):
    """How many steps are `times` diffusion times across an aggregate of this radius.

    The concentration profile inside an aggregate relaxes on tau = R^2/D, and eta is
    a property of that profile, so a run shorter than a few tau is reporting a
    transient. Everything needed is in the case file:

        dt = nu_ref dx^2 / D_ref,  nu_ref = cs2 (tau_ref - 1/2),  cs2 = 1/4 for D3Q7

    so the answer is in iterations rather than seconds, which is what the solver wants.
    """
    def num(tag, default=None, block=None):
        m = re.search(r"<%s>([^<]*)</%s>" % (tag, tag), block if block else xml)
        return float(m.group(1)) if m else default
    tau_ref = num("tau", 0.8)
    dx = num("dx", 10.0) * 1e-6                     # micrometres in the file, metres here
    b0 = substrate_block(xml, 0)
    dpore = dfilm = 1.5e-9
    if b0:
        dpore = num("in_pore", 1.5e-9, b0.group(0))
        dfilm = num("in_biofilm", dpore, b0.group(0))
    nu_ref = 0.25 * (tau_ref - 0.5)
    dt = nu_ref * dx * dx / dpore
    R = radius_voxels * dx
    tdiff = R * R / dfilm
    return max(500, int(math.ceil(times * tdiff / dt)))


def read_last_point(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        rows = [r for r in csv.DictReader(l for l in f if not l.startswith("#"))]
    return {k: float(v) for k, v in rows[-1].items()} if rows else None


def least_squares_line(xs, ys):
    """y = a + b x, by hand, so this file needs nothing installed."""
    n = len(xs)
    if n < 2:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    a = my - b * mx
    ss_tot = sum((y - my) ** 2 for y in ys)
    ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return a, b, r2


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--radii", default="4,6,8", help="aggregate radii, voxels")
    ap.add_argument("--bulk", default="1e-5,3e-4,1e-3",
                    help="bulk SULFIDE, mol/L: the product that sets F_T at the bulk")
    ap.add_argument("--iters", type=int, default=0,
                    help="iterations per run; 0 = set it from the diffusion time")
    ap.add_argument("--diffusion-times", type=float, default=4.0,
                    help="how many R^2/D to run when --iters is 0")
    ap.add_argument("--nx", type=int, default=32)
    ap.add_argument("--ny", type=int, default=32)
    ap.add_argument("--nz", type=int, default=10)
    ap.add_argument("--solver", default="./complab")
    ap.add_argument("--out", default="upscale_sweep.csv")
    a = ap.parse_args()

    radii = [float(v) for v in a.radii.split(",")]
    bulks = [float(v) for v in a.bulk.split(",")]

    if not os.path.exists(a.solver):
        print("no solver at %s -- run ./pipeline.sh first, which builds it" % a.solver,
              file=sys.stderr)
        return 1
    base = open("CompLaB.xml").read()
    if "<upscaling>" not in base:
        print("CompLaB.xml has no <upscaling> block, so there is nothing to collect.",
              file=sys.stderr)
        return 1

    os.makedirs("input", exist_ok=True)
    rows = []
    for R in radii:
        geom = "sweep_R%g.dat" % R
        cells = build_geometry(a.nx, a.ny, a.nz, R)
        with open(os.path.join("input", geom), "w") as f:
            for c in cells:
                f.write("%d\n" % c)
        nbio = cells.count(3)
        for C in bulks:
            tag = "R%g_C%g" % (R, C)
            out_dir = "sweep_%s" % tag
            os.makedirs(out_dir, exist_ok=True)
            iters = a.iters or diffusion_iterations(base, R, a.diffusion_times)
            xml = "CompLaB_%s.xml" % tag
            open(xml, "w").write(
                make_case(base, a.nx, a.ny, a.nz, geom, C, iters, out_dir + "/",
                          max(1, iters // 8)))
            print("== R = %g voxels (%d aggregate voxels), bulk HS = %g mol/L, %d iterations"
                  % (R, nbio, C, iters))
            r = subprocess.run([a.solver, xml], capture_output=True, text=True)
            if r.returncode != 0:
                print("   solver failed:\n" + r.stdout[-1500:], file=sys.stderr)
                continue
            steady = "NOT STEADY" not in r.stdout
            p = read_last_point(os.path.join(out_dir, "upscaling.csv"))
            if not p:
                print("   no upscaling record written", file=sys.stderr)
                continue
            p["radius_voxels"] = R
            p["bulk_setpoint"] = C
            p["steady"] = 1 if steady else 0
            rows.append(p)
            print("   phi %.4g   F_T(bulk) %.4g   eta %.4g   eta_classical %.4g   %s"
                  % (p["thiele"], p["F_T_bulk"], p["eta"], p["eta_classical"],
                     "steady" if steady else "STILL MOVING"))

    if not rows:
        print("no runs produced a record", file=sys.stderr)
        return 1

    keys = sorted(rows[0].keys())
    with open(a.out, "w") as f:
        f.write("# CompLB3D pore-to-aggregate sweep\n")
        f.write("# eta is measured; eta_classical is the sphere with first-order kinetics.\n")
        f.write("# Their ratio is what the gate does to the classical result.\n")
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print("\n" + "=" * 78)
    print("THE UPSCALING LAW")
    print("=" * 78)
    used = [r for r in rows if r["steady"] and r["eta_classical"] > 0]
    if len(used) < 2:
        print("  Fewer than two runs reached steady state, so there is nothing to fit.")
        print("  Raise --iters and run it again. %d of %d runs were still moving."
              % (len(rows) - len(used), len(rows)))
        print("\n  wrote %s" % a.out)
        return 0

    xs = [r["F_T_bulk"] for r in used]
    ys = [r["eta"] / r["eta_classical"] for r in used]
    fit = least_squares_line(xs, ys)
    print("  A continuum model that uses the classical Thiele curve is wrong by the")
    print("  factor below, and the factor depends on how close the BULK sits to the")
    print("  energy threshold. That dependence is what the gate adds and what no")
    print("  first-order treatment contains.\n")
    print("    %-10s %-12s %-10s %-12s %-10s" % ("R (vox)", "bulk HS", "phi", "F_T(bulk)",
                                                 "eta/eta_cl"))
    for r in used:
        print("    %-10g %-12g %-10.4g %-12.4g %-10.4g"
              % (r["radius_voxels"], r["bulk_setpoint"], r["thiele"],
                 r["F_T_bulk"], r["eta"] / r["eta_classical"]))
    if fit:
        a0, b0, r2 = fit
        print("\n    eta / eta_classical  =  %.4g  +  %.4g * F_T(bulk)     (R^2 = %.4f)"
              % (a0, b0, r2))
        print("\n  Read that as: at F_T(bulk) = 1, where the bulk is far above the")
        print("  threshold, the correction is %.3f. As the bulk approaches the" % (a0 + b0))
        print("  threshold the aggregate loses its core entirely and the correction")
        print("  falls towards %.3f." % a0)
    print("\n  wrote %s" % a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
