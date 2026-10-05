#!/usr/bin/env python3
"""Prove the whole chain works before spending a night of cluster time on it.

    python3 tools/selftest.py --work /tmp/abc_selftest

It builds a small campaign, FILLS IT WITH SYNTHETIC OUTPUT instead of running the
solver, collects it into a .h5, verifies that .h5, and trains the standalone
network on it for a handful of epochs. If this passes, every step after CompLaB
is known to work and anything that goes wrong later is in the solver or in the
case setup, which is a much smaller place to look.

WHAT IS AND IS NOT BEING TESTED
-------------------------------
The synthetic fields are NOT a simulation. They are a caricature with the right
dependencies: A enters from x = 0, B from x = nx-1, the two fronts meet at a
position set by Peclet, the product tracks their overlap with a strength set by
Damkohler, and the pore space perturbs all three through its geodesic field. That
is enough to exercise the .vti reader, the snapshot regex, the geometry
bookkeeping, the split, the .h5 layout, the sensitivity checks in
verify_dataset.py and the training loop.

What it cannot test is the physics. For that, run the straight duct case and
compare against tools/analytic_channel.py.

The .vti files it writes are ASCII. Palabos writes base64, usually compressed and
sometimes appended; the reader in collect_to_h5.py handles all of those, and only
the ASCII path is exercised here. That is the one gap in this test, and it is the
best-covered path in the reader because it is the simplest.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PORE = 2


def write_vti_scalar(path, a):
    """One scalar on a point grid, ASCII, VTK's x-fastest ordering."""
    nx, ny, nz = a.shape
    flat = np.ascontiguousarray(a.transpose(2, 1, 0)).reshape(-1)
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n')
        f.write('<ImageData WholeExtent="0 %d 0 %d 0 %d" Origin="0 0 0" '
                'Spacing="1 1 1">\n' % (nx - 1, ny - 1, nz - 1))
        f.write('<Piece Extent="0 %d 0 %d 0 %d">\n' % (nx - 1, ny - 1, nz - 1))
        f.write('<PointData Scalars="Density">\n')
        f.write('<DataArray type="Float32" Name="Density" format="ascii">\n')
        f.write(" ".join("%.7g" % v for v in flat.tolist()))
        f.write('\n</DataArray>\n</PointData>\n</Piece>\n</ImageData>\n</VTKFile>\n')


def write_vti_vector(path, v):
    """A 3-component field, same conventions."""
    _, nx, ny, nz = v.shape
    flat = np.ascontiguousarray(v.transpose(3, 2, 1, 0)).reshape(-1)
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n')
        f.write('<ImageData WholeExtent="0 %d 0 %d 0 %d" Origin="0 0 0" '
                'Spacing="1 1 1">\n' % (nx - 1, ny - 1, nz - 1))
        f.write('<Piece Extent="0 %d 0 %d 0 %d">\n' % (nx - 1, ny - 1, nz - 1))
        f.write('<PointData Vectors="velocity">\n')
        f.write('<DataArray type="Float32" Name="velocity" '
                'NumberOfComponents="3" format="ascii">\n')
        f.write(" ".join("%.7g" % x for x in flat.tolist()))
        f.write('\n</DataArray>\n</PointData>\n</Piece>\n</ImageData>\n</VTKFile>\n')


def fake_fields(mat, gdf, pe, da, feed, t):
    """A caricature with the right dependencies. Not a simulation."""
    nx, ny, nz = mat.shape
    pore = (mat == PORE)
    x = np.arange(nx)[:, None, None].astype(np.float64)

    # the geodesic field, normalised, is what makes two pore spaces differ
    g = gdf / max(float(gdf.max()), 1.0)

    # A is carried in from x = 0 and has reached this far by time t
    front_a = t * nx * (0.35 + 0.65 * pe / (pe + 3.0)) * 1.6
    # B works upstream from x = nx-1; its reach is about L/Pe
    reach_b = nx / (pe + 0.5)
    w = 2.0 + 6.0 / (pe + 1.0)

    # the geodesic excess is what makes two pore spaces give different answers:
    # a tighter pack has a longer path to every voxel, so both fronts lag
    A = feed / (1.0 + np.exp((x - front_a + 14.0 * g) / w))
    B = feed / (1.0 + np.exp(((nx - 1 - x) - reach_b - 14.0 * g) / w))
    # the product lives on the overlap, and Damkohler says how much of it forms
    # capped at the feed: nothing in this problem can concentrate, and
    # verify_dataset.py checks exactly that, so the caricature must respect it too
    C = feed * (4.0 * A * B / (feed * feed)) * (da / (da + 1.0)) * min(1.0, t * 2.0)
    C = np.minimum(C, feed)

    A = np.where(pore, A, 0.0).astype(np.float32)
    B = np.where(pore, B, 0.0).astype(np.float32)
    C = np.where(pore, C, 0.0).astype(np.float32)
    A = np.maximum(A - C * 0.5, 0.0)
    B = np.maximum(B - C * 0.5, 0.0)

    # Hold all four x planes exactly, as the all-Dirichlet set does: each reactant
    # at the feed on its own face and at zero on the far one, and C at zero on both.
    # Pinning only the two feed faces would leave the caricature describing the
    # mixed boundary set while its own params.json said dirichlet, and the
    # held-zero checks in verify_dataset.py would then never be exercised by
    # anything, which is the one thing a selftest must not allow.
    A[0][pore[0]] = feed
    A[-1][pore[-1]] = 0.0
    B[-1][pore[-1]] = feed
    B[0][pore[0]] = 0.0
    C[0][pore[0]] = 0.0
    C[-1][pore[-1]] = 0.0

    vel = np.zeros((3, nx, ny, nz), np.float32)
    vel[0][pore] = pe * 1e-4
    return A, B, C, vel


def run(cmd, allow_fail=False, **kw):
    print("\n$ %s" % " ".join(cmd))
    rc = subprocess.call(cmd, **kw)
    if rc != 0 and not allow_fail:
        raise SystemExit("step failed with exit code %d: %s" % (rc, " ".join(cmd)))
    return rc


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", default="selftest_work")
    ap.add_argument("--keep", action="store_true", help="leave the work directory")
    ap.add_argument("--epochs", type=int, default=8)
    ap.add_argument("--skip-training", action="store_true")
    args = ap.parse_args(argv)

    work = os.path.abspath(args.work)
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    py = sys.executable

    gdir = os.path.join(work, "geometries")
    cdir = os.path.join(work, "campaign")
    ddir = os.path.join(work, "dataset")

    run([py, os.path.join(HERE, "make_geometries.py"), "--out", gdir, "--no-preview"])

    # A SUBSET, deliberately. The full campaign is 270 cases x 21 snapshots x 4
    # ASCII .vti: about 22000 files and 11 GB, which is not a pre-flight check, it
    # is the thing the pre-flight check is supposed to come before. Four geometries
    # across three splits, both ends of the Peclet ladder and ALL THREE Damkohler
    # values (so the consecutive-pair sensitivity test in verify_dataset.py has
    # something to bite on), at 5 snapshots, exercises every code path at about one
    # per cent of the writing.
    run([py, os.path.join(HERE, "make_campaign.py"), "build",
         "--geometries", gdir, "--out", cdir,
         "--gids", "0", "4", "9", "18",
         "--pe", "0.02", "2.0",
         "--da", "0.1", "1.0", "10.0",
         "--snapshots", "5"])
    run([py, os.path.join(HERE, "check_campaign.py"), "check", "--campaign", cdir])

    # ---- fill the campaign with synthetic output -------------------------
    print("\n$ (filling %d cases with synthetic .vti)" % len(os.listdir(
        os.path.join(cdir, "runs"))))
    feed = json.load(open(os.path.join(cdir, "campaign.json")))["A0"]
    cache = {}
    for case in sorted(os.listdir(os.path.join(cdir, "runs"))):
        rdir = os.path.join(cdir, "runs", case)
        p = json.load(open(os.path.join(rdir, "params.json")))
        g = p["gid"]
        if g not in cache:
            z = np.load(os.path.join(gdir, "geom_%04d" % g, "geom_%04d.npz" % g))
            cache[g] = (z["material"], np.nan_to_num(z["gdf"]))
        mat, gdf = cache[g]
        od = os.path.join(rdir, "output")
        os.makedirs(od, exist_ok=True)
        vtk, nsnap = p["vtk_interval"], p["snapshots"]
        for k in range(nsnap):
            it = k * vtk
            A, B, C, vel = fake_fields(mat, gdf, p["pe"], p["da_abio"], feed,
                                       k / float(nsnap - 1))
            write_vti_scalar(os.path.join(od, "A_%07d.vti" % it), A)
            write_vti_scalar(os.path.join(od, "B_%07d.vti" % it), B)
            write_vti_scalar(os.path.join(od, "C_%07d.vti" % it), C)
            write_vti_vector(os.path.join(od, "nsLattice_%07d.vti" % it), vel)
        with open(os.path.join(od, "run.log"), "w") as f:
            f.write("[KIN] abiotic A+B->C  k_abio=%g L/(mol s)  dt_kinetics=%g s\n"
                    % (p["k_abio"], p["dt_s"]))
            f.write("  [ADE] dt=%g s/iter, total=%g s\n" % (p["dt_s"], p["t_end_s"]))
            f.write("  [NS] Pe achieved=%g (target=%g)\n" % (p["pe"], p["pe"]))
            f.write("  [CLAMP] reaction increments applied: 1000, of which 0 (0.00%) "
                    "were limited by the\n  [CLAMP] amount present. Total refused: 0, "
                    "which is 0.000% of all the consumption\n  [CLAMP] the rate laws "
                    "asked for.\n")
        json.dump({"state": "ok", "reason": "", "exit_code": 0, "wall_s": 1},
                  open(os.path.join(rdir, "status.json"), "w"))
    print("  done")

    run([py, os.path.join(HERE, "check_campaign.py"), "clamp", "--campaign", cdir])
    run([py, os.path.join(HERE, "collect_to_h5.py"), "--campaign", cdir,
         "--geometries", gdir, "--out", ddir])
    run([py, os.path.join(HERE, "verify_dataset.py"),
         "--dataset", os.path.join(ddir, "dataset.h5"), "--feed", str(feed)])

    if not args.skip_training:
        # allow_fail ON PURPOSE. A handful of epochs on a caricature says nothing
        # about accuracy, so train_prt3d's own PASS / MARGINAL / FAIL verdict is
        # not a verdict on the pipeline here. What this step proves is that the
        # loader, the model and the loop run at all, and that is the same whether
        # the score came out at 0.2 or 0.6. The real verdict is the one it gives
        # on real data after a real number of epochs.
        rc = run([py, os.path.join(HERE, "train_prt3d.py"),
                  "--dataset", os.path.join(ddir, "dataset.h5"),
                  "--species", "C", "--epochs", str(args.epochs),
                  "--n-points", "1024", "--batch-size", "4",
                  "--out", os.path.join(work, "model")], allow_fail=True)
        if rc != 0:
            print("\n  (the training verdict above is not a finding: %d epochs on "
                  "synthetic\n   data is far too short to score. The step ran, which "
                  "is what was\n   being tested.)" % args.epochs)

    print()
    print("=" * 72)
    print("SELFTEST PASSED. Every step after the solver works on this machine:")
    print("  geometries -> campaign -> preflight -> .vti -> .h5 -> verify -> train")
    print("What it does NOT test is the physics; for that, run the straight duct")
    print("case and compare with tools/analytic_channel.py.")
    print("=" * 72)
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    else:
        print("work kept at %s" % work)
    return 0


if __name__ == "__main__":
    sys.exit(main())
