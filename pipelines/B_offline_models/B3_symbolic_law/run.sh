#!/usr/bin/env bash
# B3 - search for a rate law.
#   in : a table, one column per variable plus one target column
#   out: a .sym text file, read at start-up, no rebuild needed
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"

DATA=${DATA:-samples.csv}
TARGET=${TARGET:-growth}
VARS=${VARS:-acetate}
POP=${POP:-600}
GENS=${GENS:-60}
OUT=${OUT:-discovered.sym}

python "$ROOT/tools/fit_symbolic.py" --data "$DATA" --target "$TARGET" \
       --inputs "$VARS" --pop "$POP" --gens "$GENS" --depth 6 --seed 1 \
       --out "$OUT"

cat <<'EOF'

What came back is a Pareto set, one expression per node count, not one answer.
Take the elbow: the shortest expression whose error is acceptable.

Then finish the file by hand:
  - write the other species as multiples of the fitted rate, so the
    stoichiometric ratios cannot drift
  - add one  range  line per variable.  It is enforced, not advisory.
EOF
