#!/usr/bin/env python3
"""Build the pore spaces the ABC sweep runs on.

    python3 tools/make_geometries.py --out geometries

Writes 31 geometries, each in its own directory, in the exact layout
PRT-DeepONet-3D's own tools expect:

    geometries/
      geom_0000/geom_0000.npz     material, material_jung, gdf, edt, meta
      geom_0000/geometry.dat      what CompLaB3D reads
      geom_0000/geometry.dims     nx ny nz, for a human
      geom_0000/geom_0000.png     three orthogonal slices
      ...
      geom_0900/                  the straight duct, the analytic reference case
      index.csv  index.json       gid, porosity, split, one row each

    30 sphere packs      one grain family, ten porosities, three seeds each
    one straight duct    run on its own, compared against a closed form

That layout is not decoration. 3D/tools/collect_complab_output.py reads
    <geometries>/geom_%04d/geom_%04d.npz
and pulls `material`, `gdf` and `edt` straight out of it, so a campaign built
from these files can be collected by the tool already in the GUI, unchanged, as
well as by tools/collect_to_h5.py in this package.

WHAT A .dat FILE IS, EXACTLY
----------------------------
Plain integers, one per line, no header. Nothing in the file records its own
shape, so a file of the wrong length is not refused: readGeometry() stops when
it has what it wants and ignores the rest, leaving a run with a truncated pore
space that still looks plausible. tools/check_campaign.py checks the length
against the XML for exactly that reason.

The ordering is x slowest, z fastest:

    for x: for y: for z: write value

which is NumPy C order over an array shaped (nx, ny, nz), and which is Palabos
IndexOrdering::forward as readGeometry() reads it. It is what the shipped
CompLaB examples write (examples/03_abiotic_kinetics/preprocess.py) and what
PRT's own build_geometry_3d.py writes. It is NOT what src/complab3d_geometry.hh
writes -- that in-solver generator packs x fastest, a pre-existing disagreement
with its own reader. Do not use <generate> in the XML for this sweep; use these
files.

THE FILE IS nx SLICES AND THE XML SAYS nx
-----------------------------------------
complab_functions.hh does nx += 2 after reading <nx>, and readGeometry() fills
x = 1 .. nx-2 from the file and duplicates the first and last slice into the two
ghost columns. So <nx>32</nx> means a file 32 slices wide and a transport domain
32 voxels long. Write 32, declare 32.

MATERIAL CODES, AND WHY A GRAIN HAS TWO OF THEM
-----------------------------------------------
    2   open pore, water
    1   inert wall, bounce-back
    0   declared solid phase, no dynamics

A grain is code 0 inside and code 1 in the one-voxel shell that touches pore.
That is not cosmetic. Code 0 voxels are given NoDynamics: they are not simulated
at all, and a code 0 voxel face-adjacent to pore is a hole the lattice streams
into and never streams back out of. The bounce-back shell is the wall. The
generator asserts that no code 0 voxel touches pore before it writes anything.

It is also the convention every PRT-DeepONet-3D geometry already uses
(build_geometry_3d.py), and `material` is an input the network sees, so a dataset
built with a different convention would not be comparable with the ones already
collected.

THE FOUR FACES NOBODY GIVES A BOUNDARY CONDITION
------------------------------------------------
x = 0 and x = nx-1 get the inlet and outlet named per substrate in the XML. The
four sides in y and z get nothing at all, and that is not a boundary condition:
the advection-diffusion lattices stream off the block wherever those faces are
open pore and read back an envelope nothing updates.

So the generator draws the confining wall into the geometry itself: y = 0,
y = ny-1, z = 0 and z = nz-1 are bounce-back in every file. That is what the
solver's own note recommends over <outer_faces>sealed</outer_faces>, because the
wall is then visible in the file, visible in the preview, and part of the
geometry the network is shown rather than something done to it afterwards.

The pore space is a duct of cross-section (ny-2) x (nz-2) packed with grains.
POROSITY IS REPORTED OVER THAT DUCT INTERIOR, because the interior is the packing
and the walls are the container. The whole-block figure is in index.csv too.

ISOLATED POCKETS ARE REMOVED
----------------------------
Only the pore cluster that touches both x = 0 and x = nx-1 is kept; anything else
becomes wall. A pocket with no connection to either face holds water that no flow
and no boundary condition ever reaches. In this sweep every initial concentration
is zero, so such a pocket would sit at exactly zero for the whole run and teach
the network that some pores are simply dead. Dead ENDS are kept: they are part of
the spanning cluster, they fill by diffusion, and they are a real part of how a
pore space behaves.

A seed whose spanning cluster is empty does not percolate and is rejected; the
generator says so and tries the next seed.
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

try:
    from scipy import ndimage
except ImportError:                       # SciPy is used for exactly one thing
    ndimage = None

SOLID, WALL, PORE = 0, 1, 2          # CompLaB <material_numbers>

# --- the sweep's fixed geometry family ---------------------------------------
NX = NY = NZ = 32
GRAIN_RADIUS = 3.0
# Ten porosity levels, three independent packings each, so 30 geometries.
# campaign01 used five levels and two packings, which gave the surrogate only two
# held-out porosities to be judged on and no way to tell a porosity effect from a
# single packing's quirks. Three packings per level separates those.
POROSITIES = (0.32, 0.35, 0.38, 0.41, 0.44, 0.47, 0.50, 0.53, 0.56, 0.59)
SEEDS_PER_POROSITY = 3
CHANNEL_GID = 900

# Which porosities go to which split. Held out BY POROSITY, so a validation or test
# geometry is not merely a different packing of a structure the network has already
# seen: no porosity level appears in two splits.
#
# Both validation levels and both test levels are INTERIOR to the training range
# [0.32, 0.59]:
#     0.35 sits between training 0.32 and 0.38
#     0.50 sits between training 0.47 and 0.53
#     0.41 sits between training 0.38 and 0.44
#     0.56 sits between training 0.53 and 0.59
# so this asks for interpolation, which is the honest question for an operator
# network. Extrapolating below 0.32 or above 0.59 is a different experiment and
# deserves its own geometries rather than a corner of this one.
SPLIT_OF_POROSITY = {0.32: "train", 0.35: "val",   0.38: "train",
                     0.41: "test",  0.44: "train", 0.47: "train",
                     0.50: "val",   0.53: "train", 0.56: "test",
                     0.59: "train"}


# ------------------------------------------------------------------ helpers ---
def interior(a):
    """The duct interior: everything the confining wall does not occupy."""
    return a[:, 1:-1, 1:-1]


def seal_faces(a, code=WALL):
    a[:, 0, :] = code
    a[:, -1, :] = code
    a[:, :, 0] = code
    a[:, :, -1] = code
    return a


def touches_pore(a):
    """True where a voxel is face-adjacent to at least one pore voxel."""
    pore = (a == PORE)
    pad = np.pad(pore, 1, mode="constant", constant_values=False)
    touch = np.zeros_like(pore)
    for ax in range(3):
        for s in (-1, 1):
            sl = [slice(1, -1)] * 3
            sl[ax] = slice(1 + s, pad.shape[ax] - 1 + s)
            touch |= pad[tuple(sl)]
    return touch


def apply_interface_shell(a):
    """Grain interior -> SOLID, the shell that touches pore -> WALL.

    Called last, after the pore space is final, because the shell is defined by
    what is pore and removing a pocket changes that.
    """
    out = np.where(a == PORE, PORE, SOLID).astype(np.int32)
    out[(out != PORE) & touches_pore(out)] = WALL
    return out


def no_bare_solid_touching_pore(a):
    return not ((a == SOLID) & touches_pore(a)).any()


# ---------------------------------------------------------------- geometry ---
def sphere_pack(n_grains, seed, radius=GRAIN_RADIUS):
    """n_grains overlapping spheres, centres uniform, wrapped in x only.

    Wrapped in x because the inlet and outlet planes are where the chemistry
    enters, and a grain-free margin there would put a clear channel in front of
    every case. Not wrapped in y or z: those are walls, and a grain may sit
    against one.
    """
    rng = np.random.default_rng(seed)
    a = np.full((NX, NY, NZ), PORE, dtype=np.int32)

    cx = rng.uniform(0, NX, n_grains)
    cy = rng.uniform(-radius, NY + radius, n_grains)
    cz = rng.uniform(-radius, NZ + radius, n_grains)

    xs = np.arange(NX)[:, None, None]
    ys = np.arange(NY)[None, :, None]
    zs = np.arange(NZ)[None, None, :]
    r2 = radius * radius

    for i in range(n_grains):
        dx = np.abs(xs - cx[i])
        dx = np.minimum(dx, NX - dx)                 # periodic in x
        d2 = dx * dx + (ys - cy[i]) ** 2 + (zs - cz[i]) ** 2
        a[d2 <= r2] = WALL

    return seal_faces(a)


def spanning_cluster(a):
    """Pore voxels reachable from x = 0, and whether they reach x = nx-1.

    Six-connected flood fill, iterative so a 32^3 grid cannot blow the stack.
    """
    pore = (a == PORE)
    seen = np.zeros_like(pore)
    stack = [(0, y, z) for y in range(NY) for z in range(NZ) if pore[0, y, z]]
    for p in stack:
        seen[p] = True
    while stack:
        x, y, z = stack.pop()
        for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                           (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            u, v, w = x + dx, y + dy, z + dz
            if 0 <= u < NX and 0 <= v < NY and 0 <= w < NZ \
                    and pore[u, v, w] and not seen[u, v, w]:
                seen[u, v, w] = True
                stack.append((u, v, w))
    return seen, bool(seen[NX - 1].any())


def clean(a):
    """Keep only the spanning cluster, then set the interface shell."""
    seen, perc = spanning_cluster(a)
    if not perc:
        return a, False
    out = a.copy()
    out[(a == PORE) & (~seen)] = WALL
    return apply_interface_shell(out), True


def porosity(a):
    inner = interior(a)
    return float((inner == PORE).sum()) / float(inner.size)


def block_porosity(a):
    return float((a == PORE).sum()) / float(a.size)


def build_at_porosity(target, seed, tol=0.004, verbose=True):
    """Bisect, then scan locally, on grain count until the CLEANED pack hits the
    target porosity.

    Measuring the cleaned porosity rather than the raw one matters: removing
    isolated pockets lowers porosity, and at the tight end it lowers it by more
    than the spacing between two targets.
    """
    lo, hi = 0, 16
    while True:
        a, perc = clean(sphere_pack(hi, seed))
        if not perc or porosity(a) < target:
            break
        lo = hi
        hi *= 2
        if hi > 8192:
            raise RuntimeError("cannot reach porosity %.3f at seed %d" % (target, seed))

    best = None
    for _ in range(24):
        mid = (lo + hi) // 2
        if mid == lo:
            break
        a, perc = clean(sphere_pack(mid, seed))
        phi = porosity(a) if perc else -1.0
        if perc and (best is None or abs(phi - target) < abs(best[1] - target)):
            best = (mid, phi, a)
        if not perc or phi < target:
            hi = mid
        else:
            lo = mid

    # Bisection lands within one grain of the crossing, but one grain is about
    # 0.4% of porosity here and the targets are 3% apart, so the local scan is
    # what actually gets inside the tolerance.
    if best is not None:
        for n in range(max(1, best[0] - 3), best[0] + 4):
            if n == best[0]:
                continue
            a, perc = clean(sphere_pack(n, seed))
            if not perc:
                continue
            phi = porosity(a)
            if abs(phi - target) < abs(best[1] - target):
                best = (n, phi, a)

    if best is None:
        return None
    n, phi, a = best
    if verbose:
        flag = "" if abs(phi - target) <= tol else "   <-- off target by %+.4f" % (phi - target)
        print("    %4d grains  porosity %.4f (target %.2f)%s" % (n, phi, target, flag))
    return {"grid": a, "grains": n, "porosity": phi}


def straight_duct():
    """No grains at all: a 30 x 30 open duct inside the confining wall.

    The analytic reference. At Pe = 0 the counter-diffusion of A and B has a
    closed form; with flow the duct has a known laminar profile.
    tools/analytic_channel.py compares against both.
    """
    return apply_interface_shell(seal_faces(np.full((NX, NY, NZ), PORE, dtype=np.int32)))


# -------------------------------------------------- the two distance fields ---
def geodesic_from_inlet(a):
    """Tortuous distance from the inlet plane THROUGH the pore space, in voxels.

    Six-neighbour relaxation, not fast marching, and deliberately so. Relaxation
    steps only along axes, so its paths run about 18 per cent longer than the
    eikonal solution -- but every PRT-DeepONet dataset computes it this way, and
    the trunk learns THAT scale. A geometry whose distance field is a fifth
    smaller than everything the model trained on gets a quietly worse prediction
    with no error and no warning. Consistency beats smoothness.

    Same algorithm as build_practice_dataset.geodesic(), reproduced here so this
    package stands on its own.
    """
    inf = np.float32(1e9)
    pore = (a == PORE)
    d = np.where(pore, inf, np.nan).astype(np.float32)
    d[0][pore[0]] = 0.0
    for _ in range(4 * a.shape[0]):
        prev = d.copy()
        for ax in range(3):
            for sh in (1, -1):
                nb = np.roll(d, sh, axis=ax)
                sl = [slice(None)] * 3
                sl[ax] = 0 if sh > 0 else -1
                nb[tuple(sl)] = inf
                d = np.where(pore, np.fmin(d, np.nan_to_num(nb, nan=inf) + 1.0), np.nan)
        if np.nanmax(np.abs(np.nan_to_num(d - prev, nan=0.0))) < 1e-6:
            break
    return np.nan_to_num(d, nan=0.0).astype(np.float32)


def _edt1d(f):
    """The lower envelope of the parabolas f(q) + (x-q)^2, along the last axis.

    Felzenszwalb and Huttenlocher's exact distance transform, which is what makes
    the separable three-pass version below exact rather than an approximation.
    Vectorised over every other axis.
    """
    n = f.shape[-1]
    flat = f.reshape(-1, n).copy()
    out = np.empty_like(flat)
    v = np.zeros(n, np.intp)
    z = np.empty(n + 1, np.float64)
    for r in range(flat.shape[0]):
        row = flat[r]
        k = 0
        v[0] = 0
        z[0], z[1] = -np.inf, np.inf
        for q in range(1, n):
            while True:
                s = ((row[q] + q * q) - (row[v[k]] + v[k] * v[k])) / (2.0 * q - 2.0 * v[k])
                if s <= z[k] and k > 0:
                    k -= 1
                    continue
                break
            k += 1
            v[k] = q
            z[k] = s
            z[k + 1] = np.inf
        k = 0
        for q in range(n):
            while z[k + 1] < q:
                k += 1
            out[r, q] = (q - v[k]) ** 2 + row[v[k]]
    return out.reshape(f.shape)


def euclidean_to_solid(a):
    """Straight-line distance to the nearest non-pore voxel.

    SciPy when it is there; an exact separable transform when it is not. This is
    the only thing in the package SciPy was needed for, and requiring a whole
    dependency for one call is how a script comes to be unrunnable on the machine
    it is most wanted on.
    """
    pore = (a == PORE)
    if ndimage is not None:
        return ndimage.distance_transform_edt(pore).astype(np.float32)
    big = np.float64(a.size) ** 2
    d = np.where(pore, big, 0.0).astype(np.float64)
    d = _edt1d(d)                                   # along z
    d = _edt1d(d.transpose(0, 2, 1)).transpose(0, 2, 1)     # along y
    d = _edt1d(d.transpose(2, 1, 0)).transpose(2, 1, 0)     # along x
    return np.sqrt(d).astype(np.float32)


def to_jung_codes(a):
    """CompLaB (0 solid, 1 interface, 2 pore) -> Jung (0 solid, 1 pore, 2 interface)."""
    j = np.zeros_like(a)
    j[a == PORE] = 1
    j[a == WALL] = 2
    return j


def tortuosity(gdf, a):
    """Mean geodesic length of the outlet plane, over the straight-line length."""
    out = gdf[NX - 1][a[NX - 1] == PORE]
    out = out[out > 0]
    if out.size == 0:
        return float("nan")
    return float(out.mean() / (NX - 1))


# --------------------------------------------------------------------- I/O ---
def write_dat(path, a):
    """x slowest, z fastest. C order over (nx, ny, nz) is exactly that."""
    assert a.shape == (NX, NY, NZ)
    flat = np.ascontiguousarray(a, dtype=np.int32).reshape(-1)
    # newline="\n" so a file written on Windows is byte-identical to one written
    # on Linux. Palabos would parse either, but a geometry that changes size when
    # it crosses machines is a thing nobody should have to discover.
    with open(path, "w", newline="\n") as f:
        f.write("\n".join(str(int(v)) for v in flat.tolist()))
        f.write("\n")
    return flat.size


def read_dat(path):
    """The inverse. Used by the round-trip check below and by check_campaign.py."""
    return np.loadtxt(path, dtype=np.int32).reshape((NX, NY, NZ))


def preview(path, a):
    """Three orthogonal slices, so a wrong axis order is visible rather than latent."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap, BoundaryNorm
    except Exception:
        return False
    cmap = ListedColormap(["#222222", "#888888", "#ffffff"])   # solid, wall, pore
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    fig, ax = plt.subplots(1, 3, figsize=(9.6, 3.4))
    for k, (img, title) in enumerate((
            (a[:, :, NZ // 2].T, "z mid-slice   x right, y up"),
            (a[:, NY // 2, :].T, "y mid-slice   x right, z up"),
            (a[NX // 2, :, :].T, "x mid-slice   y right, z up"))):
        ax[k].imshow(img, origin="lower", cmap=cmap, norm=norm, interpolation="nearest")
        ax[k].set_title(title, fontsize=8)
        ax[k].set_xticks([]); ax[k].set_yticks([])
    fig.suptitle("white pore (2)    grey bounce-back (1)    black no-dynamics solid (0)",
                 fontsize=8, y=0.03)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return True


def emit(outdir, gid, grid, row, want_preview=True):
    """One geometry, in the layout collect_complab_output.py reads."""
    d = os.path.join(outdir, "geom_%04d" % gid)
    os.makedirs(d, exist_ok=True)

    assert no_bare_solid_touching_pore(grid), \
        "geom_%04d has a no-dynamics voxel touching pore" % gid

    gdf = geodesic_from_inlet(grid)
    edt = euclidean_to_solid(grid)
    row = dict(row)
    row["tortuosity"] = round(tortuosity(gdf, grid), 4)

    np.savez_compressed(os.path.join(d, "geom_%04d.npz" % gid),
                        material=grid.astype(np.uint8),
                        material_jung=to_jung_codes(grid).astype(np.uint8),
                        gdf=gdf, edt=edt, meta=json.dumps(row))
    n = write_dat(os.path.join(d, "geometry.dat"), grid)
    with open(os.path.join(d, "geometry.dims"), "w") as f:
        f.write("%d %d %d\n" % grid.shape)
    if want_preview:
        preview(os.path.join(d, "geom_%04d.png" % gid), grid)
    row["values"] = n
    return row


# -------------------------------------------------------------------- main ---
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="geometries", help="output directory")
    ap.add_argument("--no-preview", action="store_true", help="skip the PNG previews")
    args = ap.parse_args(argv)

    os.makedirs(args.out, exist_ok=True)
    rows = []
    gid = 0

    print("sphere packs, grain radius %.0f voxels, %d x %d x %d"
          % (GRAIN_RADIUS, NX, NY, NZ))
    # SEEDS ARE BLOCKED PER LEVEL, 100 apart, and never shared.
    #
    # The obvious scheme, seed = 1000 + round(100 * porosity), gives bases only 3
    # apart now that the levels are 0.03 apart, while each level consumes a seed per
    # ATTEMPT rather than per success. One non-percolating seed at a tight porosity
    # would push that level into its neighbour's first seed, and sphere_pack draws
    # its centres from that seed, so the two packs would be nested: one pack's grains
    # a prefix of the other's. Across a split boundary (0.35 validation next to 0.38
    # train, 0.56 test next to 0.59 train) that is leakage no gid-based split check
    # could ever see. A block of 100 per level cannot collide.
    seed_used = {}
    for level_index, phi_target in enumerate(POROSITIES):
        print("  porosity %.2f  ->  %s" % (phi_target, SPLIT_OF_POROSITY[phi_target]))
        made = 0
        seed = 1000 + 100 * level_index
        while made < SEEDS_PER_POROSITY:
            if seed >= 1000 + 100 * level_index + 100:
                raise SystemExit(
                    "porosity %.2f exhausted its 100 seeds without %d percolating "
                    "packs. Widen the block or loosen the tolerance."
                    % (phi_target, SEEDS_PER_POROSITY))
            if seed in seed_used:
                raise SystemExit("seed %d reused (porosity %.2f and %.2f). Two packs "
                                 "drawn from one seed are nested, not independent."
                                 % (seed, seed_used[seed], phi_target))
            got = build_at_porosity(phi_target, seed)
            if got is None:
                print("    seed %d does not percolate at this porosity, trying the next"
                      % seed)
                seed += 1
                continue
            seed_used[seed] = phi_target
            rows.append(emit(args.out, gid, got["grid"], {
                "gid": gid, "name": "geom_%04d" % gid, "kind": "spheres",
                "nx": NX, "ny": NY, "nz": NZ,
                "grain_radius": GRAIN_RADIUS, "grains": got["grains"], "seed": seed,
                "porosity_target": phi_target,
                "porosity": round(got["porosity"], 6),
                "porosity_block": round(block_porosity(got["grid"]), 6),
                "split": SPLIT_OF_POROSITY[phi_target]},
                want_preview=not args.no_preview))
            gid += 1
            made += 1
            seed += 1

    duct = straight_duct()
    rows.append(emit(args.out, CHANNEL_GID, duct, {
        "gid": CHANNEL_GID, "name": "geom_%04d" % CHANNEL_GID, "kind": "straight_duct",
        "nx": NX, "ny": NY, "nz": NZ,
        "grain_radius": 0.0, "grains": 0, "seed": 0,
        "porosity_target": 1.0,
        "porosity": round(porosity(duct), 6),
        "porosity_block": round(block_porosity(duct), 6),
        "split": "analytic"}, want_preview=not args.no_preview))
    print("straight duct  gid %d, porosity %.4f over the interior"
          % (CHANNEL_GID, porosity(duct)))

    fields = ["gid", "name", "kind", "nx", "ny", "nz", "values", "grain_radius",
              "grains", "seed", "porosity_target", "porosity", "porosity_block",
              "tortuosity", "split"]
    with open(os.path.join(args.out, "index.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(args.out, "index.json"), "w") as f:
        json.dump({"nx": NX, "ny": NY, "nz": NZ,
                   "codes": {"solid_no_dynamics": SOLID, "bounce_back": WALL, "pore": PORE},
                   "ordering": "x slowest, z fastest (C order over (nx,ny,nz))",
                   "geometries": rows}, f, indent=2)

    # --- the round trip, run every time rather than trusted ------------------
    bad = 0
    for r in rows:
        p = os.path.join(args.out, r["name"], "geometry.dat")
        back = read_dat(p)
        npz = np.load(os.path.join(args.out, r["name"], r["name"] + ".npz"))
        if back.shape != (NX, NY, NZ):
            print("  ROUND TRIP FAILED (shape) %s" % r["name"]); bad += 1; continue
        if (back != npz["material"]).any():
            print("  ROUND TRIP FAILED (.dat vs .npz) %s" % r["name"]); bad += 1; continue
        if abs(porosity(back) - r["porosity"]) > 1e-6:      # index.csv holds 6 places
            print("  ROUND TRIP FAILED (porosity) %s" % r["name"]); bad += 1; continue
        if (back[:, 0, :] == PORE).any() or (back[:, -1, :] == PORE).any() \
                or (back[:, :, 0] == PORE).any() or (back[:, :, -1] == PORE).any():
            print("  ROUND TRIP FAILED (open y/z face) %s" % r["name"]); bad += 1; continue
        if not no_bare_solid_touching_pore(back):
            print("  ROUND TRIP FAILED (bare solid) %s" % r["name"]); bad += 1; continue
        # slice k of the file must be plane x = k
        first = np.loadtxt(p, dtype=np.int32)[:NY * NZ].reshape((NY, NZ))
        if (first != back[0]).any():
            print("  ROUND TRIP FAILED (ordering) %s" % r["name"]); bad += 1

    print()
    print("wrote %d geometries to %s/  (%d integers each)"
          % (len(rows), args.out, NX * NY * NZ))
    print("index.csv and index.json list gid, porosity, tortuosity and split")
    if bad:
        print("%d geometries FAILED the round trip" % bad)
        return 1
    print("all %d round-tripped: shape, .npz agreement, porosity, sealed y/z faces,"
          % len(rows))
    print("no bare solid against pore, and x-slice ordering")
    return 0


if __name__ == "__main__":
    sys.exit(main())
