#!/usr/bin/env python3
"""
20_upscaling  --  PRE-PROCESSING

Builds the pore space: one microbial aggregate, a circular cylinder spanning the
full depth in z, sitting in open water.
That is the whole geometry, because that is the whole question — how much of a
resolved aggregate is actually working, and what single rate a continuum model
should use in its place.

  reads   nothing
  writes  input/geometry.dat   one integer per voxel, z fastest

Self-contained: standard library only, nothing from tools/. Run it before the
solver, or let pipeline.sh run it for you:

    python3 preprocess.py

THE DOMAIN IS CLOSED ON FOUR SIDES. x = 0 and x = NX-1 carry the boundary
conditions named per substrate in CompLaB.xml. The solver gives the other four
faces nothing at all — not a wall, not a symmetry plane, not periodicity — so a
face left as open water has no boundary condition, and the lattice streams off
the block there and reads back an envelope nothing updates. The two y faces and
the two z faces are therefore inert wall, and NY and NZ are two larger than the
water they hold. Adding the layers rather than converting water into wall keeps
the aggregate and the open water at exactly the volumes this case declares.

WHY A ROUND AGGREGATE, AND WHY ONE OF THEM. The effectiveness factor has an
exact classical result for a sphere with first-order kinetics, printed by the
solver next to the measured value. Starting from a geometry close to the one the
textbook can also handle is what makes the disagreement meaningful.

[v1.3] What is drawn is a circular CYLINDER, not a sphere: the membership test
below is on x and y only, inside the z loop, so every z layer between the walls
is filled to the same radius. A sphere is not available at this size in any case,
the drawn diameter being 12 voxels with only 8 z layers of water. Part of any
departure between the measured and the classical curve is therefore a shape
effect, and not only the chemistry.

The material numbers are the ones declared in CompLaB.xml:
    1  inert wall, never changes
    2  open pore, water — the BULK
    3  the aggregate: pore holding biomass
"""
import os
import sys

NX, NY, NZ = 32, 34, 10          # NY and NZ include the added wall layers
RADIUS = 6.0                     # voxels; 6 x 10 um = 60 um radius, 120 um across

PORE, WALL, BIO = 2, 1, 3
OUT = os.path.join("input", "geometry.dat")


def build():
    """One integer per voxel, indexed [x][y][z]."""
    cx, cy = (NX - 1) / 2.0, (NY - 1) / 2.0
    g = [[[PORE] * NZ for _ in range(NY)] for _ in range(NX)]
    for x in range(NX):
        for y in range(NY):
            for z in range(NZ):
                if y in (0, NY - 1) or z in (0, NZ - 1):
                    g[x][y][z] = WALL
                    continue
                if (x - cx) ** 2 + (y - cy) ** 2 <= RADIUS * RADIUS:
                    g[x][y][z] = BIO
    return g


def main():
    g = build()

    os.makedirs("input", exist_ok=True)
    with open(OUT, "w") as f:
        for x in range(NX):
            for y in range(NY):
                for z in range(NZ):
                    f.write("%d\n" % g[x][y][z])

    counts = {}
    for x in range(NX):
        for y in range(NY):
            for z in range(NZ):
                counts[g[x][y][z]] = counts.get(g[x][y][z], 0) + 1
    total = NX * NY * NZ

    label = {1: "inert wall", 2: "open water (the bulk)", 3: "the aggregate"}
    print("wrote %s  (%d x %d x %d = %d voxels)" % (OUT, NX, NY, NZ, total))
    for m in sorted(counts):
        print("  material %d  %6d  %5.1f%%   %s"
              % (m, counts[m], 100.0 * counts[m] / total, label.get(m, "?")))

    nbio = counts.get(BIO, 0)
    # The radius the solver will use for the Thiele modulus when <radius> is 0:
    # the sphere of the same volume as the aggregate actually drawn. Printing it
    # here means the number in the report can be checked against the geometry
    # rather than taken on trust.
    r_eff = (3.0 * nbio / (4.0 * 3.141592653589793)) ** (1.0 / 3.0)
    print("  aggregate      %d voxels, equivalent sphere radius %.2f voxels = %.0f um"
          % (nbio, r_eff, r_eff * 10))
    print("  drawn radius   %.1f voxels = %.0f um" % (RADIUS, RADIUS * 10))
    print("     (the two differ because the sphere is drawn on a cubic lattice and")
    print("      spans the full depth in z; the solver uses the equivalent radius,")
    print("      which is the one the volume actually supports)")

    if counts.get(PORE, 0) == 0:
        print("  no open water left: the aggregate fills the domain", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
