#!/usr/bin/env bash
# A1 - build the pore space.
#   in : nothing (or a segmented image stack)
#   out: geometry.dat, one integer per voxel: 0 solid phase, 1 wall, 2 pore, 3+ microbe seed
#
# [v1.3] This script used to invoke tools/geometry.py with --type, --width, --out
# and --inspect, none of which that program defines, and without the subcommand it
# requires -- so it exited 2 before doing anything. The calls below are the real
# interface: `geometry.py create <kind> ...` and `geometry.py inspect <dat>`.
# Run `python3 tools/geometry.py create --help` for the full option list.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"

NX=${NX:-128}; NY=${NY:-64}; NZ=${NZ:-64}
KIND=${KIND:-channel}          # channel cylinders fracture layered random spheres
APERTURE=${APERTURE:-20}       # channel/fracture opening, in voxels
WALLS=${WALLS:-y}              # which faces get a confining wall: none x y z xy xz yz all
OUT=${OUT:-geometry.dat}

python3 "$ROOT/tools/geometry.py" create "$KIND" \
        --nx "$NX" --ny "$NY" --nz "$NZ" \
        --aperture "$APERTURE" --walls "$WALLS" -o "$OUT"

# never hand a geometry to the solver without looking at it first
python3 "$ROOT/tools/geometry.py" inspect "$OUT" --nx "$NX" --ny "$NY" --nz "$NZ"
