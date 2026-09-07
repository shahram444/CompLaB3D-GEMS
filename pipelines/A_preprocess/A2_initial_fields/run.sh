#!/usr/bin/env bash
# A2 - starting fields, for cases that need more than a uniform value.
#   in : geometry.dat
#   out: the same geometry with a microbe seed region written into it
#
# [v1.3] This script used to call tools/geometry.py with --field, --shape, --face,
# --value, --like and --out, none of which exist, and it described writing a
# separate biomass0.dat that the solver has no tag to read. What the solver
# actually does is seed biomass on a MATERIAL NUMBER: <material_numbers><microbeN>
# names a code, and every voxel carrying it starts with that microbe's
# <initial_densities>. `geometry.py seed` writes that code into a box, which is
# the real form of this step.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"

GEOM=${GEOM:-geometry.dat}
OUT=${OUT:-geometry_seeded.dat}
NX=${NX:-128}; NY=${NY:-64}; NZ=${NZ:-64}
CODE=${CODE:-3}                                  # the <material_numbers><microbe0> code
BOX=${BOX:-"1 $((NX-2)) 1 4 1 $((NZ-2))"}        # x0 x1 y0 y1 z0 z1: a film on the y=0 wall

python3 "$ROOT/tools/geometry.py" seed "$GEOM" \
        --nx "$NX" --ny "$NY" --nz "$NZ" \
        --code "$CODE" --box $BOX -o "$OUT"

python3 "$ROOT/tools/geometry.py" inspect "$OUT" --nx "$NX" --ny "$NY" --nz "$NZ"
