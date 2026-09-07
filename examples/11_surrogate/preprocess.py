#!/usr/bin/env python3
"""
11_surrogate  --  PRE-PROCESSING

Builds the pore space this case runs on: the same geometry as 09 and 10, so the surrogate can be held against the linear program.

  reads   nothing
  writes  input/geometry.dat   24 x 26 x 8 = 4992 voxels, one integer per line
                               ([v1.3] 24 x 24 x 6 is the shape BEFORE the wall
                               padding explained under pad_closed_faces() below,
                               which grows y and z by two each. The file on disk
                               holds 4992 values and CompLaB.xml declares
                               ny 26, nz 8.)

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
BIOFILM_3 = 3


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

    # --- the initial biomass patch for material 3 ------------------------
    for x in range(NX):
        for y in range(NY):
            if ((x % 4) == 0 and 4 <= y <= 6):
                for z in range(NZ):
                    if g[x][y][z] == PORE:
                        g[x][y][z] = BIOFILM_3

    return g


def percolates(g):
    """Can water get from the x = 0 face to the x = NX-1 face? Breadth first."""
    open_ = lambda x, y, z: g[x][y][z] != WALL
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


def pad_closed_faces(g):
    """Wrap the domain in the wall it was always assumed to have.

    The solver conditions two faces of the box: x = 0 and x = NX-1 get the inlet
    and the outlet named per substrate in CompLaB.xml. The four sides -- y = 0,
    y = NY-1, z = 0, z = NZ-1 -- get nothing at all. Not a wall, not a symmetry
    plane, not periodicity. A lattice Boltzmann field with no condition on a face
    streams off the edge of the block there and reads back an envelope nothing
    updates, so what crosses is undefined.

    Every generator in this repository used to draw its grains and forget the box,
    and this one used to be among them. It now ADDS a wall layer rather than
    converting pore into wall: a face that is already closed is left alone, and a
    face that is open gains one new layer of inert wall outside everything the
    case declared. The pore space, the grains and the patches are untouched, so
    porosity, the mineral inventory and every documented count are exactly what
    they were. What changes is only NY and NZ, which grow by two where a pair of
    faces needed closing, and those are printed at the end for CompLaB.xml.

    Padding rather than converting matters most on a thin slab. Walling the two
    z faces of a six-deep domain would spend a third of the pore space on the
    boundary condition; adding two layers spends none of it.
    """
    nx, ny, nz = len(g), len(g[0]), len(g[0][0])
    closed = (WALL, 0)
    yopen = any(g[x][y][z] not in closed
                for x in range(nx) for y in (0, ny - 1) for z in range(nz))
    zopen = any(g[x][y][z] not in closed
                for x in range(nx) for y in range(ny) for z in (0, nz - 1))
    py, pz = (1 if yopen else 0), (1 if zopen else 0)
    if not (py or pz):
        return g, nx, ny, nz

    NY2, NZ2 = ny + 2 * py, nz + 2 * pz
    out = [[[WALL] * NZ2 for _ in range(NY2)] for _ in range(nx)]
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                out[x][y + py][z + pz] = g[x][y][z]
    return out, nx, NY2, NZ2

def main():
    #   NY and NZ below are the padded dimensions, which is what the geometry
    #   file now holds and what CompLaB.xml has to declare. The module constants
    #   NX/NY/NZ stay as the case's own dimensions, so build() is unchanged.
    g, NX, NY, NZ = pad_closed_faces(build())

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

    fluid = sum(n for m, n in counts.items() if m != 1)
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
