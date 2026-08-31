#!/usr/bin/env python3
"""
15_precip_and_dissolution  --  PRE-PROCESSING

Builds the pore space this case runs on: a slot lined with one phase while a second forms from solution.

  reads   nothing
  writes  input/geometry.dat   24 x 24 x 6 = 3456 voxels, one integer per line

Self-contained: no imports beyond the standard library, nothing from tools/.
Run it before the solver, or let pipeline.sh run it for you:

    python3 preprocess.py

THE GEOMETRY IS A STAGGERED SLOT. Grain blocks alternate between the two
walls every four voxels along x, so the flow has to weave and the field is
not trivially symmetric. Small enough to run in seconds, awkward enough to
be a real test of a boundary condition.

The material numbers are the ones declared in CompLaB.xml:
    0  a declared solid phase (dissolves, or is added to)
    1  inert wall, never changes
    2  open pore, water
    3, 4, ...  open pore holding an initial biomass patch
"""
import os
import sys
from collections import deque

NX, NY, NZ = 24, 24, 6
OUT = os.path.join("input", "geometry.dat")

PORE, WALL = 2, 1
MINERAL = 0                      # the declared phase, lining the grain surfaces


def build():
    """One integer per voxel, indexed [x][y][z]."""
    g = [[[PORE] * NZ for _ in range(NY)] for _ in range(NX)]

    # --- the staggered grain blocks -------------------------------------
    for x in range(NX):
        if x % 4 == 0:                 # block against the y = 0 wall
            for y in range(0, 4):
                for z in range(NZ):
                    g[x][y][z] = WALL
        if x % 4 == 3:                 # block against the y = NY-1 wall
            for y in range(NY - 4, NY):
                for z in range(NZ):
                    g[x][y][z] = WALL

    # --- the reactive phase: the grain layer that FACES the pore ---------
    # One voxel thick, on the inner face of each block only. Mineral buried
    # inside a grain can never be reached by the water, so declaring it would
    # cost solve time and change nothing. The rest of each block stays inert
    # wall, which is what <bounce_back> in CompLaB.xml binds to.
    for x in range(NX):
        if x % 4 == 0:                 # this block's inner face is at y = 3
            for z in range(NZ):
                g[x][3][z] = MINERAL
        if x % 4 == 3:                 # ... and at y = NY-4 for the other
            for z in range(NZ):
                g[x][NY - 4][z] = MINERAL

    return g


def percolates(g):
    """Can water get from the x = 0 face to the x = NX-1 face? Breadth first."""
    open_ = lambda x, y, z: g[x][y][z] != WALL and g[x][y][z] != MINERAL
    seen = set()
    q = deque((0, y, z) for y in range(NY) for z in range(NZ) if open_(0, y, z))
    seen.update(q)
    while q:
        x, y, z = q.popleft()
        if x == NX - 1:
            return True
        for dx, dy, dz in ((1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)):
            a, b, c = x + dx, y + dy, z + dz
            if (0 <= a < NX and 0 <= b < NY and 0 <= c < NZ
                    and (a, b, c) not in seen and open_(a, b, c)):
                seen.add((a, b, c))
                q.append((a, b, c))
    return False


def main():
    g = build()

    os.makedirs("input", exist_ok=True)
    with open(OUT, "w") as f:
        for x in range(NX):                 # x-slices, y then z inside each
            for y in range(NY):
                for z in range(NZ):
                    f.write("%d\n" % g[x][y][z])

    counts = {}
    for x in range(NX):
        for y in range(NY):
            for z in range(NZ):
                counts[g[x][y][z]] = counts.get(g[x][y][z], 0) + 1
    total = NX * NY * NZ

    label = {0: "declared solid phase", 1: "inert wall", 2: "open pore"}
    print("wrote %s  (%d x %d x %d = %d voxels)" % (OUT, NX, NY, NZ, total))
    for m in sorted(counts):
        name = label.get(m, "biomass patch, microbe%d" % (m - 3))
        print("  material %d  %6d  %5.1f%%   %s" % (m, counts[m], 100.0 * counts[m] / total, name))

    fluid = sum(n for m, n in counts.items() if m != 1 and m != 0)
    print("  porosity   %.4f" % (float(fluid) / total))

    if percolates(g):
        print("  percolates from inlet to outlet: yes")
    else:
        print("  percolates from inlet to outlet: NO", file=sys.stderr)
        print("  A pressure drop across a sealed domain is not a well-posed", file=sys.stderr)
        print("  problem and the solver will refuse it. Fix the geometry.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
