#!/usr/bin/env bash
#
# Assemble a runnable working directory for one example.
#
# An example directory holds EVERYTHING that case needs: its configuration, its
# pore space, its chemistry, the metabolic model it reads, the training code
# that produced whatever it was trained on, and its own pre- and
# post-processing.  Nothing has to be fetched from elsewhere in the tree.
#
# So this script does very little.  It copies the case, adds the shared kinetics
# defaults underneath it where the case does not override them, adds the solver
# sources so the case builds on its own, and gets out of the way:
#
#   ./scripts/setup_case.sh 13_precipitation  run/precip
#   cd run/precip && ./pipeline.sh
#
# You can also just work inside examples/13_precipitation/ directly.  Assembling
# a copy keeps the example clean so the next run starts from a known state.

set -euo pipefail

usage() {
    cat <<'EOF'
usage: setup_case.sh <example> [destination]

  <example>      a directory name under examples/, with or without the number:
                 13_precipitation, or just precipitation
  [destination]  where to assemble it; default run/<example>

Then:
  cd <destination>
  ./pipeline.sh                                  # the whole chain

or step by step:
  python3 preprocess.py                          # the pore space
  ./offline.sh                                   # where the case has one
  cmake -B build -S . && cmake --build build -j
  ./complab CompLaB.xml 2>&1 | tee output/run.log
  python3 postprocess.py                         # is it worth believing
EOF
}

[[ $# -lt 1 || $1 == -h || $1 == --help ]] && { usage; exit 1; }

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WANT="$1"

CASE_DIR=""
for d in "$ROOT"/examples/[0-9]*; do
    b="$(basename "$d")"
    if [[ "$b" == "$WANT" || "${b#*_}" == "$WANT" ]]; then CASE_DIR="$d"; break; fi
done
if [[ -z "$CASE_DIR" ]]; then
    echo "setup_case.sh: no example matches '$WANT'." >&2
    echo "Available:" >&2
    for d in "$ROOT"/examples/[0-9]*; do echo "  $(basename "$d")" >&2; done
    exit 1
fi

CASE="$(basename "$CASE_DIR")"
DEST="${2:-$ROOT/run/$CASE}"
mkdir -p "$DEST"

echo "Assembling $CASE"

# 1. the shared default chemistry, so the case always compiles
#
#    All three are needed by EVERY case, not only the ones that use them.  The
#    solver includes defineKinetics.hh, defineAbioticKinetics.hh and
#    surrogateModel.hh unconditionally -- they are user-editable files, and a
#    file that is #included cannot be optional -- so a case assembled without
#    any one of them does not compile.  A case that uses none of the three
#    still gets the defaults, which do nothing.
cp "$ROOT/config/kinetics/defineKinetics.default.hh"        "$DEST/defineKinetics.hh"
cp "$ROOT/config/kinetics/defineAbioticKinetics.default.hh" "$DEST/defineAbioticKinetics.hh"
cp "$ROOT/src/surrogateModel.hh"                            "$DEST/surrogateModel.hh"
echo "  defaults      defineKinetics.hh, defineAbioticKinetics.hh, surrogateModel.hh"

# 2. the case itself, whole: input, models, training, its own scripts
cp -r "$CASE_DIR/." "$DEST/"
rm -f "$DEST/README.md"
for f in preprocess.py postprocess.py offline.sh pipeline.sh; do
    [[ -f "$DEST/$f" ]] && chmod +x "$DEST/$f"
done
echo "  the case      $(cd "$CASE_DIR" && find . -type f | wc -l | tr -d ' ') files, copied whole"

# 3. this case's own chemistry sits on top of the defaults
if [[ -d "$DEST/kinetics" ]]; then
    for f in "$DEST/kinetics"/*.hh; do
        cp "$f" "$DEST/$(basename "$f")"
        echo "  own chemistry $(basename "$f")"
    done
    rm -rf "$DEST/kinetics"
fi

# 4. the solver sources, so the case builds without leaving this directory
mkdir -p "$DEST/src" "$DEST/output"
cp "$ROOT"/src/*.hh "$ROOT"/src/*.cpp "$DEST/src/" 2>/dev/null || true
cp "$ROOT/CMakeLists.txt" "$DEST/" 2>/dev/null || true

# 4b. the COBRApy bridge module, which the embedded interpreter imports BY NAME
#     from <src_path> at start-up.  It is a .py and lives in tools/, so the *.hh
#     and *.cpp copy above does not pick it up -- and without it the COBRApy path
#     fails at start-up with "failed to import module complab3d_cobrapy", after
#     the build has succeeded and the geometry has loaded.  Example 10 shipped in
#     exactly that state: it assembled, it compiled, and it could not run.
#     Copied unconditionally: it is 20 kB and costs nothing in a case that never
#     asks for COBRApy.
cp "$ROOT/tools/complab3d_cobrapy.py" "$DEST/src/" 2>/dev/null || true

echo
echo "Assembled $CASE in $DEST"
echo
echo "The whole pipeline in one command:"
echo "  cd $DEST && ./pipeline.sh"
echo
echo "or step by step:"
echo "  cd $DEST"
echo "  python3 preprocess.py"
[[ -f "$DEST/offline.sh" ]] && echo "  ./offline.sh"
echo "  cmake -B build -S . && cmake --build build -j"
echo "  ./complab CompLaB.xml 2>&1 | tee output/run.log"
echo "  python3 postprocess.py"
echo
if [[ -f "$CASE_DIR/README.md" ]]; then
    echo "What this case shows: $CASE_DIR/README.md"
fi
