#!/usr/bin/env python3
"""
19_thermodynamic_gate  --  PRE-PROCESSING

Builds the pore space this case runs on: a single spherical microbial
aggregate sitting in open water.

  reads   nothing
  writes  input/geometry.dat   32 x 32 x 10 = 10240 voxels, one integer per line

Self-contained: nothing beyond the standard library, nothing from tools/.
Run it before the solver, or let pipeline.sh run it for you:

    python3 preprocess.py

WHY A SPHERE AND NOT A SLOT.  The gate this case exists to demonstrate closes
where a reaction's own products accumulate faster than diffusion removes them.
That needs a region with a long diffusive path to open water, which is exactly
what the interior of a millimetre-scale aggregate is, and exactly what the
pore-throat geometries in examples 01 to 18 are not.  ANME-SRB consortia are
observed as roughly spherical aggregates 100 to 300 um across (Orphan et al.
2002); the one below is 120 um at dx = 10 um.

The material numbers are the ones declared in CompLaB.xml:
    0  a declared solid phase
    1  inert wall, never changes
    2  open pore, water
    3  open pore holding the initial biomass patch
"""
import os
import sys

NX, NY, NZ = 32, 32, 10

PORE = 2
WALL = 1
BIO = 3

# The aggregate: a sphere of radius R voxels centred in x and y, spanning the
# full depth in z so the third direction adds no escape route the picture does
# not show.
CX, CY = (NX - 1) / 2.0, (NY - 1) / 2.0
R = 6.0                     # voxels; 6 x 10 um = 60 um radius, 120 um across


def build():
    """Return the mask as a flat list in the solver's z-fastest order.

    NZ below is two larger than the aggregate slab needs, because the two z faces
    are wall. The solver conditions only x = 0 and x = NX-1, so a z face left as
    open water has no boundary condition at all: the lattice streams off the block
    there and reads back an envelope nothing updates. Nothing in this geometry
    varies with z, so the closure costs nothing except the two layers -- which are
    ADDED here rather than taken out of the water, so the aggregate and the open
    water are exactly the volumes the case declares."""
    cells = []
    for x in range(NX):
        for y in range(NY):
            for z in range(NZ):
                # The four faces the solver gives no boundary condition to are
                # inert wall, so the domain is a closed slab of water with the
                # aggregate in the middle of it. x = 0 and x = NX-1 stay open:
                # those two get the inlet and the outlet named per substrate.
                if y == 0 or y == NY - 1 or z == 0 or z == NZ - 1:
                    cells.append(WALL)
                    continue
                d2 = (x - CX) ** 2 + (y - CY) ** 2
                cells.append(BIO if d2 <= R * R else PORE)
    return cells


def report(cells):
    n = len(cells)
    bio = cells.count(BIO)
    wall = cells.count(WALL)
    pore = cells.count(PORE)
    print("  geometry      %d x %d x %d = %d voxels" % (NX, NY, NZ, n))
    print("  aggregate     %d voxels (%.1f%%), radius %.0f um" % (bio, 100.0 * bio / n, R * 10))
    print("  open water    %d voxels (%.1f%%)" % (pore, 100.0 * pore / n))
    print("  inert wall    %d voxels" % wall)
    # The deepest voxel inside the aggregate is the one the gate should shut
    # first.  Saying how far it is from open water sets the expectation.
    print("  deepest point %.0f um from open water" % (R * 10))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "input", "geometry.dat")
    os.makedirs(os.path.dirname(out), exist_ok=True)

    cells = build()
    if len(cells) != NX * NY * NZ:
        print("preprocess.py: built %d cells, expected %d" % (len(cells), NX * NY * NZ))
        return 1

    with open(out, "w") as f:
        for c in cells:
            f.write("%d\n" % c)

    print("== 19_thermodynamic_gate pre-processing")
    report(cells)
    print("  wrote         input/geometry.dat")
    return 0


if __name__ == "__main__":
    sys.exit(main())
