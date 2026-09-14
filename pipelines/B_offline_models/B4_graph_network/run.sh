#!/usr/bin/env bash
# B4 - train the graph network.
#   in : stoich.csv, samples.csv, one Damkohler number per reaction
#   out: a .gnn text file, read at start-up, no rebuild needed
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"

STOICH=${STOICH:-stoich.csv}
SAMPLES=${SAMPLES:-samples.csv}
DA=${DA:-1.0,1.0}
OUT=${OUT:-network.gnn}

python "$ROOT/tools/train_graphnet.py" --stoich "$STOICH" --data "$SAMPLES" \
       --da "$DA" --rounds 1 --width 8 --epochs 4000 --out "$OUT"

echo
echo "checking the trained network against the samples"
python "$ROOT/tests/xval_gnn.py" --net "$OUT" --samples "$SAMPLES"

cat <<'EOF'

The number that matters is not the correlation, it is the ratio between species
rates.  It should come out close to the stoichiometry without ever having been
a training target.  If it does not, the graph is wrong, not the fit.
EOF
