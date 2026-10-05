#!/usr/bin/env bash
#
# Assemble a buildable CompLaB3D tree carrying THIS package's chemistry.
#
#   ./tools/setup_build.sh <path to the v1.3.2 tree> [destination]
#
# Example, from inside this package:
#
#   ./tools/setup_build.sh ".." ~/abc_build
#   cd ~/abc_build && cmake -B build -S . -DPALABOS_ROOT=/path/to/palabos-v2.3.0
#   cmake --build build -j
#
# WHAT IT DOES, AND WHY EACH PIECE IS NEEDED
#
#   1. The solver sources, whole. complab.cpp includes forty headers out of src/
#      and there is no partial copy that builds.
#
#   2. THE PATCH. patch/apply_damkohler_feed_reference.py is applied to the COPY,
#      never to your original tree. Stock v1.3.2 builds the Damkohler banner's
#      reference composition from <initial_concentration>, and every one of those
#      is zero in this sweep, so without it the solver reports no Damkohler at
#      all for any case. The patch is idempotent and keeps a .orig beside the
#      file it changes.
#
#   3. THE THREE FILES THE SOLVER #INCLUDES BY NAME, at the tree root:
#          defineKinetics.hh          inert here, and still mandatory
#          defineAbioticKinetics.hh   A + B -> C, this package's copy
#          surrogateModel.hh          unused here, and still #included
#      complab3d_processors_part1.hh includes the first two unconditionally as
#      "../defineKinetics.hh", so they have to sit beside src/ rather than in it.
#      A tree assembled without any one of them does not compile; one assembled
#      without the KineticsStats namespace compiles and fails to link.
#
#   4. A reminder of which chemistry went in, written to CHEMISTRY.txt, because
#      the single most expensive mistake available here is building against the
#      shipped defaults and discovering weeks later that every case in a
#      Damkohler sweep ran identical chemistry.
#
# It does NOT copy a CompLaB.xml. Each case brings its own, written by
# tools/make_campaign.py, and the binary takes it as argv[1].

set -euo pipefail

usage() {
    echo "usage: setup_build.sh <path to CompLaB3D v1.3.2 tree> [destination]" >&2
    echo >&2
    echo "  the v1.3.2 tree is the folder holding src/, config/ and CMakeLists.txt" >&2
    echo "  destination defaults to ./build_tree" >&2
}

[[ $# -lt 1 || $1 == -h || $1 == --help ]] && { usage; exit 1; }

SRC="$(cd "$1" && pwd)"
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${2:-$PWD/build_tree}"

for f in src/complab.cpp src/complab_functions.hh CMakeLists.txt; do
    [[ -f "$SRC/$f" ]] || { echo "not a CompLaB3D tree: $SRC has no $f" >&2; exit 1; }
done
for f in kinetics/defineKinetics.hh kinetics/defineAbioticKinetics.hh \
         patch/apply_damkohler_feed_reference.py; do
    [[ -f "$PKG/$f" ]] || { echo "this package is incomplete: no $f" >&2; exit 1; }
done

mkdir -p "$DEST/src" "$DEST/output"
echo "assembling $DEST"

cp "$SRC"/src/*.hh "$SRC"/src/*.cpp "$DEST/src/"
cp "$SRC/CMakeLists.txt" "$DEST/"
[[ -f "$SRC/tools/complab3d_cobrapy.py" ]] && cp "$SRC/tools/complab3d_cobrapy.py" "$DEST/src/" || true
echo "  solver        $(ls "$DEST/src" | wc -l | tr -d ' ') files"

python3 "$PKG/patch/apply_damkohler_feed_reference.py" "$DEST/src/complab.cpp" \
    | sed 's/^/  patch         /'

cp "$PKG/kinetics/defineKinetics.hh"        "$DEST/defineKinetics.hh"
cp "$PKG/kinetics/defineAbioticKinetics.hh" "$DEST/defineAbioticKinetics.hh"
cp "$SRC/src/surrogateModel.hh"             "$DEST/surrogateModel.hh"
echo "  chemistry     defineKinetics.hh (inert), defineAbioticKinetics.hh (A + B -> C)"
echo "  also          surrogateModel.hh, unused and still required"

cat > "$DEST/CHEMISTRY.txt" <<EOF
This tree was assembled by $PKG/tools/setup_build.sh
from $SRC on $(date -u +%Y-%m-%dT%H:%M:%SZ).

Its chemistry is the ABC sweep's, NOT the CompLaB defaults:

  defineAbioticKinetics.hh   A + B -> C, R = k[A][B], k from PRT_KABIO
  defineKinetics.hh          every rate zero; <biotic_mode>false</biotic_mode>

The binary built here reads PRT_KABIO, PRT_DT and PRT_MAXFRAC from the
environment and prints a [KIN] line naming all three on its first reaction step.
If that line is missing from a run's log, the binary is NOT this one.

Prove it before spending cluster time:
  python3 $PKG/tools/check_campaign.py check --campaign <campaign> --complab $DEST/complab
EOF

echo
echo "assembled. Next:"
echo "  cd $DEST"
echo "  cmake -B build -S . -DPALABOS_ROOT=/path/to/palabos-v2.3.0"
echo "  cmake --build build -j"
echo
echo "The binary lands at $DEST/complab (CMAKE_RUNTIME_OUTPUT_DIRECTORY is '../')."
echo "Check the first run's log for the [KIN] line before trusting anything."
