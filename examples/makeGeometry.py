#!/usr/bin/env python3
"""
Build the tiny geometries the example cases use.

CompLB3D reads a geometry as ONE INTEGER PER LINE, x fastest, then y, then z:

    for z: for y: for x:  write mask(x, y, z)

Every example here is 24 x 24 x 6 with z periodic, which makes the domain
quasi-2D: every z-plane sees identical conditions. That is deliberate for a test
suite. It keeps the runs to seconds, and it means a result that depends on z is
a bug rather than physics.

Material numbers, the same in every case:

    0  solid        no dynamics, and the phase dissolution eats
    1  bounce_back  the channel walls
    2  pore         water
    3  microbe0 seed
    4  microbe1 seed

Run with no arguments; it writes into ./input_shared/.

Part of the CompLaB program.  GNU Affero General Public License v3 or later.
"""

import os

NX, NY, NZ = 24, 24, 6
SOLID, WALL, PORE, MIC0, MIC1 = 0, 1, 2, 3, 4


def blank():
    """Open channel: bounce-back walls on the two y faces, pore everywhere else.

    The walls are on y, not x, because flow is driven along x by <delta_P> and
    the x faces have to stay open for the pressure boundary to mean anything."""
    g = [[[PORE for _ in range(NX)] for _ in range(NY)] for _ in range(NZ)]
    for z in range(NZ):
        for x in range(NX):
            g[z][0][x] = WALL
            g[z][NY - 1][x] = WALL
    return g


def add_grains(g, centres, radius, code=SOLID):
    """Drop cylindrical grains through the full z depth.

    Cylinders rather than spheres, again to keep the domain quasi-2D: a sphere
    would make the middle z-planes different from the outer ones, and then the
    "every plane must agree" check stops working."""
    for z in range(NZ):
        for y in range(NY):
            for x in range(NX):
                for (cx, cy) in centres:
                    if (x - cx) ** 2 + (y - cy) ** 2 <= radius * radius:
                        g[z][y][x] = code
    return g


def add_patch(g, x0, x1, y0, y1, code):
    """Seed a microbe over a rectangular patch, pore voxels only.

    Seeding on top of a wall would put biomass where the solver never visits,
    and the initial biomass would silently vanish."""
    for z in range(NZ):
        for y in range(max(y0, 0), min(y1 + 1, NY)):
            for x in range(max(x0, 0), min(x1 + 1, NX)):
                if g[z][y][x] == PORE:
                    g[z][y][x] = code
    return g


def write(g, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for z in range(NZ):
            for y in range(NY):
                for x in range(NX):
                    f.write("%d\n" % g[z][y][x])
    counts = {}
    for z in range(NZ):
        for y in range(NY):
            for x in range(NX):
                counts[g[z][y][x]] = counts.get(g[z][y][x], 0) + 1
    print("%-40s %s" % (os.path.basename(path),
                        "  ".join("mat%d:%d" % (k, counts[k]) for k in sorted(counts))))


def main():
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "input_shared")

    # 1. open channel: flow, diffusion, abiotic chemistry, equilibrium
    write(blank(), os.path.join(out, "channel.dat"))

    # 2. channel with four grains: anything where the flow field should be
    #    interesting, and the substrate for precipitation and dissolution
    grains = add_grains(blank(), [(7, 8), (7, 16), (16, 8), (16, 16)], 3)
    write(grains, os.path.join(out, "grains.dat"))

    # 3. one microbe seeded as a biofilm patch on the lower wall
    one = add_patch(blank(), 4, 9, 1, 3, MIC0)
    write(one, os.path.join(out, "one_microbe.dat"))

    # 4. two microbes, one on each wall, so they compete through the solute
    #    field rather than by touching
    two = add_patch(blank(), 4, 9, 1, 3, MIC0)
    two = add_patch(two, 14, 19, NY - 4, NY - 2, MIC1)
    write(two, os.path.join(out, "two_microbes.dat"))

    print("\n%d voxels per file (%d x %d x %d)" % (NX * NY * NZ, NX, NY, NZ))


if __name__ == "__main__":
    main()
