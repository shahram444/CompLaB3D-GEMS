#!/usr/bin/env python3
"""
Get a pore space into CompLB3D, from wherever it is now.

WHY THIS EXISTS

CompLB3D reads its geometry as a bare list of integers with no header, and
until now the only way to produce one was to write the loop yourself and get
the ordering right. That is a silent failure: the wrong loop order gives a file
of exactly the right length that describes a transposed pore space, and nothing
downstream complains.

This tool covers the whole path in: build a geometry, or import one from a CT
image stack or a numpy array; check it before you spend a day of compute on it;
seed microbes into it; and export it to VTK so you can look at it.

    geometry.py create   ...   make one from a recipe
    geometry.py import   ...   from an image stack, .npy or raw binary
    geometry.py inspect  ...   porosity, percolation, dead pore space
    geometry.py seed     ...   add microbe patches to an existing geometry
    geometry.py export   ...   .dat to .vti, to look at in ParaView

MATERIAL NUMBERS, the convention used throughout:

    0  solid        a grain: no dynamics, and what dissolution eats
    1  bounce_back  a wall: reflects the fluid
    2  pore         water
    3, 4, ...       microbe seeds

EXAMPLES

    # a 100^3 random porous medium at 35% porosity, with walls on y
    geometry.py create random --nx 100 --ny 100 --nz 100 --porosity 0.35 \\
        --walls y -o input/geometry.dat

    # a sphere pack
    geometry.py create spheres --nx 80 --ny 80 --nz 80 --radius 8 --n 40 \\
        -o input/geometry.dat

    # from a segmented CT stack: one image per z-slice, dark = solid
    geometry.py import ct_slices/ --threshold 128 -o input/geometry.dat

    # always do this before running anything big
    geometry.py inspect input/geometry.dat --nx 100 --ny 100 --nz 100

Part of the CompLaB program.  GNU Affero General Public License v3 or later.
"""

import argparse
import os
import sys

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


SOLID, WALL, PORE = 0, 1, 2


# ===========================================================================
# reading and writing the .dat
# ===========================================================================
def write_dat(path, g):
    """g is (nz, ny, nx). The file is x fastest, then y, then z.

    That order is the whole contract with the solver, and it is not recorded
    anywhere in the file, so it is written in exactly one place: here."""
    if os.path.dirname(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    nz, ny, nx = g.shape
    with open(path, "w") as f:
        for z in range(nz):
            for y in range(ny):
                f.write("\n".join(str(int(v)) for v in g[z, y, :]))
                f.write("\n")
    return nx * ny * nz


def read_dat(path, nx, ny, nz):
    vals = np.loadtxt(path, dtype=int)
    if vals.size != nx * ny * nz:
        raise SystemExit(
            "ERROR: %s has %d values but nx*ny*nz = %d.\n"
            "       Either the dimensions are wrong or the file is truncated."
            % (path, vals.size, nx * ny * nz))
    return vals.reshape((nz, ny, nx))


def add_walls(g, axis):
    """Bounce-back walls on the two faces of one axis.

    Never on x by default: x is where the pressure boundary drives the flow,
    and walling it off gives a domain with no inlet."""
    if axis in ("y", "xy", "yz", "all"):
        g[:, 0, :] = WALL
        g[:, -1, :] = WALL
    if axis in ("z", "xz", "yz", "all"):
        g[0, :, :] = WALL
        g[-1, :, :] = WALL
    if axis in ("x", "xy", "xz", "all"):
        g[:, :, 0] = WALL
        g[:, :, -1] = WALL
    return g


# ===========================================================================
# generators
# ===========================================================================
def gen_channel(nx, ny, nz, **kw):
    return np.full((nz, ny, nx), PORE, dtype=int)


def gen_spheres(nx, ny, nz, radius=6, n=20, seed=0, **kw):
    """Randomly placed overlapping spheres of solid.

    Overlapping on purpose: a non-overlapping pack is a much harder packing
    problem and gives a narrower porosity range."""
    rng = np.random.default_rng(seed)
    g = np.full((nz, ny, nx), PORE, dtype=int)
    zz, yy, xx = np.mgrid[0:nz, 0:ny, 0:nx]
    for _ in range(n):
        cx, cy, cz = rng.integers(0, nx), rng.integers(0, ny), rng.integers(0, nz)
        g[(xx - cx) ** 2 + (yy - cy) ** 2 + (zz - cz) ** 2 <= radius ** 2] = SOLID
    return g


def gen_cylinders(nx, ny, nz, radius=4, n=8, seed=0, **kw):
    """Cylinders through the full z depth: a quasi-2D grain pack.

    Useful because every z-plane is identical, so the answer must not depend on
    z. That makes it a free correctness check on any run."""
    rng = np.random.default_rng(seed)
    g = np.full((nz, ny, nx), PORE, dtype=int)
    yy, xx = np.mgrid[0:ny, 0:nx]
    for _ in range(n):
        cx, cy = rng.integers(0, nx), rng.integers(0, ny)
        disc = (xx - cx) ** 2 + (yy - cy) ** 2 <= radius ** 2
        g[:, disc] = SOLID
    return g


def _smooth(a, passes=2):
    """A cheap box blur, so a random field has structure at more than one voxel.

    scipy.ndimage.gaussian_filter would be better, but this needs no scipy and
    the difference does not matter for generating a test medium."""
    for _ in range(passes):
        s = np.zeros_like(a)
        cnt = 0
        for dz in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    s += np.roll(np.roll(np.roll(a, dz, 0), dy, 1), dx, 2)
                    cnt += 1
        a = s / cnt
    return a


def gen_random(nx, ny, nz, porosity=0.4, seed=0, roughness=2, **kw):
    """Threshold a smoothed random field at whatever level gives the target
    porosity. Exact by construction, because the threshold is a quantile."""
    rng = np.random.default_rng(seed)
    f = _smooth(rng.random((nz, ny, nx)), passes=int(roughness))
    thr = np.quantile(f, 1.0 - porosity)
    g = np.where(f >= thr, PORE, SOLID)
    return g.astype(int)


def gen_fracture(nx, ny, nz, aperture=6, roughness=2.0, seed=0, **kw):
    """A rough-walled fracture through solid rock, running along x."""
    rng = np.random.default_rng(seed)
    g = np.full((nz, ny, nx), SOLID, dtype=int)
    mid = ny / 2.0
    wob = np.zeros((nz, nx))
    for k in (1, 2, 3):
        ph = rng.random() * 2 * np.pi
        wob += (roughness / k) * np.sin(2 * np.pi * k * np.arange(nx) / nx + ph)[None, :]
    yy = np.arange(ny)[None, :, None]
    centre = (mid + wob)[:, None, :]
    g[np.abs(yy - centre) <= aperture / 2.0] = PORE
    return g


def gen_layered(nx, ny, nz, layers=4, porosity=0.5, seed=0, **kw):
    """Alternating high and low porosity bands across the flow direction.

    The simplest geometry that produces preferential flow, which is what most
    people actually want to demonstrate."""
    rng = np.random.default_rng(seed)
    g = np.full((nz, ny, nx), PORE, dtype=int)
    edges = np.linspace(0, ny, layers + 1).astype(int)
    for i in range(layers):
        p = porosity * (1.4 if i % 2 == 0 else 0.5)
        sl = rng.random((nz, edges[i + 1] - edges[i], nx)) > min(p, 1.0)
        g[:, edges[i]:edges[i + 1], :][sl] = SOLID
    return g


GENERATORS = {
    "channel": gen_channel,
    "spheres": gen_spheres,
    "cylinders": gen_cylinders,
    "random": gen_random,
    "fracture": gen_fracture,
    "layered": gen_layered,
}


# ===========================================================================
# inspection
# ===========================================================================
def _label(binary):
    """Connected components of a boolean array, 6-connected.

    Uses scipy when it is there and a plain BFS when it is not, so this runs on
    a bare login node."""
    try:
        from scipy import ndimage
        return ndimage.label(binary)[0]
    except ImportError:
        pass
    lab = np.zeros(binary.shape, dtype=np.int32)
    cur = 0
    nz, ny, nx = binary.shape
    idx = np.argwhere(binary)
    seen = set()
    pos = {tuple(p): True for p in map(tuple, idx)}
    for start in map(tuple, idx):
        if start in seen:
            continue
        cur += 1
        stack = [start]
        seen.add(start)
        while stack:
            z, y, x = stack.pop()
            lab[z, y, x] = cur
            for dz, dy, dx in ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                               (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                n = (z + dz, y + dy, x + dx)
                if (0 <= n[0] < nz and 0 <= n[1] < ny and 0 <= n[2] < nx
                        and n in pos and n not in seen):
                    seen.add(n)
                    stack.append(n)
    return lab


def inspect(g, quiet=False):
    nz, ny, nx = g.shape
    total = g.size
    codes, counts = np.unique(g, return_counts=True)
    openv = g >= PORE
    porosity = openv.mean()

    print("Geometry: %d x %d x %d  =  %d voxels" % (nx, ny, nz, total))
    print("\nMaterial numbers present")
    names = {SOLID: "solid", WALL: "bounce_back", PORE: "pore"}
    for c, n in zip(codes, counts):
        nm = names.get(int(c), "microbe seed" if c > PORE else "?")
        print("  %-3d  %-14s %8d  %6.2f%%" % (c, nm, n, 100.0 * n / total))
    print("\nPorosity (everything >= 2 counts as open): %.4f" % porosity)

    problems = []

    # percolation along x, the flow direction
    lab = _label(openv)
    inlet = set(np.unique(lab[:, :, 0])) - {0}
    outlet = set(np.unique(lab[:, :, -1])) - {0}
    spanning = inlet & outlet
    if spanning:
        vol = sum(int((lab == s).sum()) for s in spanning)
        print("\nPercolation along x: YES, %d spanning cluster(s), %.1f%% of "
              "the open space" % (len(spanning), 100.0 * vol / max(openv.sum(), 1)))
    else:
        print("\nPercolation along x: NO")
        problems.append(
            "No open path from the inlet face to the outlet face. The flow "
            "solver will not converge to anything meaningful, and a pressure "
            "drop across a sealed domain is not a well-posed problem.")

    # dead pore space
    ncl = int(lab.max())
    if ncl > 1:
        sizes = np.bincount(lab.ravel())[1:]
        dead = sum(int(s) for i, s in enumerate(sizes, start=1) if i not in spanning)
        if dead:
            frac = 100.0 * dead / max(openv.sum(), 1)
            print("Isolated pore space: %d voxel(s), %.1f%% of the open space, "
                  "in %d cluster(s)" % (dead, frac, ncl - len(spanning)))
            if frac > 5:
                problems.append(
                    "%.1f%% of the pore space is isolated from the inlet. "
                    "Solutes there can only diffuse in through nothing, so it "
                    "stays at its initial concentration forever and drags "
                    "every domain-averaged number you compute." % frac)
    else:
        print("Isolated pore space: none")

    # z-independence, worth knowing because every shipped example relies on it
    if nz > 1:
        quasi2d = all((g[k] == g[0]).all() for k in range(1, nz))
        print("Quasi-2D (every z-plane identical): %s" % ("yes" if quasi2d else "no"))
        if quasi2d:
            print("  so the answer must not depend on z; that is a free "
                  "correctness check on any run")

    if problems:
        print("\nPROBLEMS")
        for p in problems:
            print("  - %s" % p)
    else:
        print("\nNo problems found.")
    return len(problems)


# ===========================================================================
def main():
    ap = argparse.ArgumentParser(
        prog="geometry.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Create, import, inspect and export CompLB3D geometries.",
        epilog="See the header of this file for worked examples.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def dims(p):
        p.add_argument("--nx", type=int, required=True)
        p.add_argument("--ny", type=int, required=True)
        p.add_argument("--nz", type=int, required=True)

    c = sub.add_parser("create", help="build a geometry from a recipe")
    c.add_argument("kind", choices=sorted(GENERATORS))
    dims(c)
    c.add_argument("--porosity", type=float, default=0.4)
    c.add_argument("--radius", type=float, default=6)
    c.add_argument("--n", type=int, default=20, help="how many grains")
    c.add_argument("--aperture", type=float, default=6)
    c.add_argument("--roughness", type=float, default=2)
    c.add_argument("--layers", type=int, default=4)
    c.add_argument("--seed", type=int, default=0)
    c.add_argument("--walls", default="y", choices=["none", "x", "y", "z", "xy", "xz", "yz", "all"])
    c.add_argument("-o", "--output", default="geometry.dat")
    c.add_argument("--vti", action="store_true", help="also write a .vti to look at")

    i = sub.add_parser("import", help="from an image stack, .npy or raw binary")
    i.add_argument("source", help="a directory of slice images, a .npy, or a raw file")
    i.add_argument("--threshold", type=float, default=None,
                   help="values at or above this are PORE (default: Otsu)")
    i.add_argument("--invert", action="store_true", help="bright means solid instead")
    i.add_argument("--raw-dims", nargs=3, type=int, metavar=("NX", "NY", "NZ"),
                   help="required for a raw binary source")
    i.add_argument("--raw-dtype", default="uint8")
    i.add_argument("--downsample", type=int, default=1, help="take every Nth voxel")
    i.add_argument("--walls", default="none", choices=["none", "x", "y", "z", "xy", "xz", "yz", "all"])
    i.add_argument("-o", "--output", default="geometry.dat")
    i.add_argument("--vti", action="store_true")

    s = sub.add_parser("inspect", help="porosity, percolation, dead pore space")
    s.add_argument("dat")
    dims(s)

    d = sub.add_parser("seed", help="add microbe patches to an existing geometry")
    d.add_argument("dat")
    dims(d)
    d.add_argument("--code", type=int, required=True, help="material number, 3 or above")
    d.add_argument("--box", nargs=6, type=int, required=True,
                   metavar=("X0", "X1", "Y0", "Y1", "Z0", "Z1"))
    d.add_argument("-o", "--output", required=True)

    e = sub.add_parser("export", help=".dat to .vti, for ParaView")
    e.add_argument("dat")
    dims(e)
    e.add_argument("-o", "--output", default=None)

    args = ap.parse_args()

    # -------------------------------------------------------------- create
    if args.cmd == "create":
        g = GENERATORS[args.kind](
            args.nx, args.ny, args.nz,
            porosity=args.porosity, radius=args.radius, n=args.n,
            aperture=args.aperture, roughness=args.roughness,
            layers=args.layers, seed=args.seed)
        if args.walls != "none":
            g = add_walls(g, args.walls)
        n = write_dat(args.output, g)
        print("Wrote %s   (%d values)" % (args.output, n))
        print()
        bad = inspect(g)
        if args.vti:
            from vtireader import write_vti
            v = os.path.splitext(args.output)[0] + ".vti"
            write_vti(v, g, name="material")
            print("\nWrote %s   open it in ParaView and Threshold on 'material'" % v)
        return 1 if bad else 0

    # -------------------------------------------------------------- import
    if args.cmd == "import":
        g = load_source(args)
        n = write_dat(args.output, g)
        print("Wrote %s   (%d values)" % (args.output, n))
        print("\nAdd these to CompLaB.xml:")
        print("    <nx>%d</nx> <ny>%d</ny> <nz>%d</nz>"
              % (g.shape[2], g.shape[1], g.shape[0]))
        print()
        bad = inspect(g)
        if args.vti:
            from vtireader import write_vti
            v = os.path.splitext(args.output)[0] + ".vti"
            write_vti(v, g, name="material")
            print("\nWrote %s" % v)
        return 1 if bad else 0

    # ------------------------------------------------------------- inspect
    if args.cmd == "inspect":
        return 1 if inspect(read_dat(args.dat, args.nx, args.ny, args.nz)) else 0

    # ---------------------------------------------------------------- seed
    if args.cmd == "seed":
        if args.code <= PORE:
            sys.exit("ERROR: a microbe seed code must be 3 or above; 0, 1 and 2 "
                     "are solid, bounce_back and pore.")
        g = read_dat(args.dat, args.nx, args.ny, args.nz)
        x0, x1, y0, y1, z0, z1 = args.box
        sel = np.zeros(g.shape, bool)
        sel[max(z0, 0):z1 + 1, max(y0, 0):y1 + 1, max(x0, 0):x1 + 1] = True
        # Only pore voxels: seeding onto a wall puts biomass where the solver
        # never looks, and it silently vanishes on the first step.
        sel &= (g == PORE)
        g[sel] = args.code
        write_dat(args.output, g)
        print("Seeded %d voxel(s) with material number %d" % (int(sel.sum()), args.code))
        if sel.sum() == 0:
            print("WARNING: nothing was seeded. The box may be outside the "
                  "domain, or every voxel in it is solid or wall.")
        print("Wrote %s" % args.output)
        return 0

    # -------------------------------------------------------------- export
    if args.cmd == "export":
        from vtireader import write_vti
        g = read_dat(args.dat, args.nx, args.ny, args.nz)
        out = args.output or (os.path.splitext(args.dat)[0] + ".vti")
        write_vti(out, g, name="material")
        print("Wrote %s" % out)
        return 0


def _otsu(a):
    """Otsu's threshold: the level that best separates the histogram in two.

    Written out rather than imported so this works without scikit-image."""
    hist, edges = np.histogram(a.ravel(), bins=256)
    hist = hist.astype(float)
    w = np.cumsum(hist)
    m = np.cumsum(hist * ((edges[:-1] + edges[1:]) / 2))
    tot = m[-1]
    denom = w * (w[-1] - w)
    denom[denom == 0] = 1.0
    between = (tot * w - m * w[-1]) ** 2 / denom
    return float((edges[:-1] + edges[1:])[np.argmax(between)] / 2)


def load_source(args):
    src = args.source

    if os.path.isdir(src):
        try:
            from PIL import Image
        except ImportError:
            sys.exit("ERROR: reading an image stack needs Pillow.\n"
                     "       python3 -m pip install pillow")
        files = sorted(f for f in os.listdir(src)
                       if f.lower().endswith((".png", ".tif", ".tiff",
                                              ".jpg", ".jpeg", ".bmp")))
        if not files:
            sys.exit("ERROR: no image files in %s" % src)
        slices = [np.asarray(Image.open(os.path.join(src, f)).convert("L"))
                  for f in files]
        vol = np.stack(slices, axis=0).astype(float)
        print("Read %d slice(s) of %d x %d from %s"
              % (len(files), vol.shape[2], vol.shape[1], src))

    elif src.endswith(".npy"):
        vol = np.load(src).astype(float)
        if vol.ndim != 3:
            sys.exit("ERROR: %s has shape %s; a 3D array is required"
                     % (src, vol.shape))
        print("Read %s, shape (nz, ny, nx) = %s" % (src, vol.shape))

    else:
        if not args.raw_dims:
            sys.exit("ERROR: a raw binary source needs --raw-dims NX NY NZ.")
        nx, ny, nz = args.raw_dims
        vol = np.fromfile(src, dtype=np.dtype(args.raw_dtype)).astype(float)
        if vol.size != nx * ny * nz:
            sys.exit("ERROR: %s holds %d values but nx*ny*nz = %d"
                     % (src, vol.size, nx * ny * nz))
        vol = vol.reshape((nz, ny, nx))
        print("Read %s as %s, shape (nz, ny, nx) = %s"
              % (src, args.raw_dtype, vol.shape))

    if args.downsample > 1:
        d = args.downsample
        vol = vol[::d, ::d, ::d]
        print("Downsampled by %d to %s" % (d, vol.shape))

    thr = args.threshold if args.threshold is not None else _otsu(vol)
    if args.threshold is None:
        print("Threshold chosen automatically (Otsu): %.4g" % thr)
    pore = vol >= thr
    if args.invert:
        pore = ~pore
    g = np.where(pore, PORE, SOLID).astype(int)
    print("Porosity after thresholding: %.4f" % pore.mean())
    if args.walls != "none":
        g = add_walls(g, args.walls)
    return g


if __name__ == "__main__":
    sys.exit(main())
