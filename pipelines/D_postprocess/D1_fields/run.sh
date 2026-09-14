#!/usr/bin/env bash
# D1 - fields and pictures.
#   in : a run folder, holding CompLaB.xml and output/
#   out: summary.csv, profiles.csv and PNG plots, in <run>/output
#
# [v1.3] This script used to pass --dir, --slice, --at, --out and --history, none
# of which tools/postprocess.py defines, and it never supplied the run folder the
# program requires as a positional argument -- so it exited 2. postprocess.py takes
# the run folder and does the slices, the histories and the plots in one pass.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"

DIR=${1:-$ROOT/run/13_precipitation}
OUT=${OUT:-}                       # defaults to <DIR>/output

if [ -n "$OUT" ]; then
    python3 "$ROOT/tools/postprocess.py" "$DIR" --output "$OUT"
    echo "wrote $OUT"
else
    python3 "$ROOT/tools/postprocess.py" "$DIR"
    echo "wrote $DIR/output"
fi
