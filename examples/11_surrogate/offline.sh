#!/usr/bin/env bash
#
# 11_surrogate  --  OFFLINE STEP
#
# Sweep the linear program, fit the network to what it returned, verify the
# exported header reproduces the trainer, install it, and look at the response
# surface.
#
# Run before the solver. pipeline.sh does it for you.

set -euo pipefail

# THE LONGEST PREPARATION IN THE REPOSITORY, and the only one ending in a
# recompile. It solves a linear program at every point of a grid, so the cost is
# GRID^2 solves. The repository ships a fitted src/surrogateModel.hh, so this
# case runs without any of the below: this is how that header was made, and how
# you would make one for a model of your own.

MODEL=${MODEL:-models/e_coli_core.xml.gz}
OBJECTIVE=${OBJECTIVE:-Biomass_Ecoli_core}

# 41 x 41 = 1681 solves, a couple of minutes. The shipped header used 141, which
# is 19881 solves; go back to that for anything you intend to publish.
GRID=${GRID:-41}

# generateTrainingData.py reads .xml, .mat and .json, not .gz.
if [[ "$MODEL" == *.gz ]]; then
    echo "expanding $MODEL"
    gunzip -kf "$MODEL"
    MODEL="${MODEL%.gz}"
fi

echo
echo "1/4  sweeping the linear program over a ${GRID} x ${GRID} grid (the slow part)"
python3 training/generateTrainingData.py "$MODEL" \
        --objective "$OBJECTIVE" \
        --exchange EX_glc__D_e --range 0.001 10   --log \
        --exchange EX_o2_e     --range 3e-5 0.5   --log \
        --grid "$GRID" -o training_data.csv

echo
echo "2/4  fitting the network to that sweep"
python3 training/trainSurrogate.py training_data.csv --name surrogate \
        --layers 10 10 10 10 --restarts 5 -o surrogate_weights.hh

echo
echo "3/4  checking the exported header reproduces the trainer"
# Do this EVERY time. A transcription error between the fitted weights and the
# C++ header is invisible in the numbers and fatal in the results.
python3 training/verifyExport.py surrogate_weights.hh

echo
echo "4/4  looking at the response surface BEFORE trusting it"
# This path enforces no training range at all. Outside the fitted box it returns
# a confident number and no warning, and a large region returning zero growth
# looks exactly like a region that does not.
python3 training/inspectSurrogate.py surrogate_weights.hh --eval 9.0 0.45 || true

cat <<'EOF'

The fitted header is surrogate_weights.hh. It is NOT installed, because
installing it changes what the solver computes and requires a rebuild. When you
are satisfied with the response surface:

    cp surrogate_weights.hh src/surrogateModel.hh
    cmake --build build -j

Until then the case runs on the header that ships with the repository.
EOF
