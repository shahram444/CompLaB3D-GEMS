#!/usr/bin/env bash
# B2 - fit the surrogate, then rebuild.  The longest preparation here, and the
# only one that ends in a recompile.
#   in : an SBML/BiGG model, and the range of uptake rates expected
#   out: src/surrogateModel.hh, compiled into the solver
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
T="$ROOT/tools/surrogate"

MODEL=${MODEL:-$ROOT/models/e_coli_core.xml.gz}
OBJECTIVE=${OBJECTIVE:-Biomass_Ecoli_core}

# generateTrainingData.py reads .xml, .mat and .json, not .gz.
if [[ "$MODEL" == *.gz ]]; then
    gunzip -kf "$MODEL"
    MODEL="${MODEL%.gz}"
fi
EX1=${EX1:-EX_glc__D_e};  LO1=${LO1:-0.001};   HI1=${HI1:-10}
EX2=${EX2:-EX_o2_e};      LO2=${LO2:-0.00003}; HI2=${HI2:-0.5}
GRID=${GRID:-141}
NAME=${NAME:-surrogate}

echo "1/4  sweeping the linear program over the grid (this is the slow part)"
python3 "$T/generateTrainingData.py" "$MODEL" --objective "$OBJECTIVE" \
       --exchange "$EX1" --range "$LO1" "$HI1" --log \
       --exchange "$EX2" --range "$LO2" "$HI2" --log \
       --grid "$GRID" -o training_data.csv

echo "2/4  fitting the network"
python3 "$T/trainSurrogate.py" training_data.csv --name "$NAME" \
       --layers 10 10 10 10 --restarts 5 -o "surrogate_weights_${NAME}.hh"

echo "3/4  checking the export reproduces the trainer"
python3 "$T/verifyExport.py" "surrogate_weights_${NAME}.hh"

echo "4/4  installing the header"
cp "surrogate_weights_${NAME}.hh" "$ROOT/src/surrogateModel.hh"

cat <<'EOF'

Now rebuild the solver:   cmake --build build -j

Before you trust the fit, evaluate it across the range you will actually visit:

    python tools/surrogate/inspectSurrogate.py src/surrogateModel.hh --eval 9.0 0.45

A large part of a fitted box returning zero growth is common and invisible in
the weights.  Find that out now rather than after a week of simulation.
EOF
