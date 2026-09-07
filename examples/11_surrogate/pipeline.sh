#!/usr/bin/env bash
#
# 11_surrogate  --  THE WHOLE PIPELINE
#
# A fitted network standing in for the linear program.
#
# Everything this case needs is in this directory. Run it from an assembled
# case:
#
#     ./scripts/setup_case.sh 11_surrogate run/mycase
#     cd run/mycase
#     ./pipeline.sh
#
# Each step is also runnable on its own, and each is a file you can read:
#
#     preprocess.py    builds the pore space
#     offline.sh       sweep the linear program, fit the network, verify it. It does
#                      NOT install the fitted header: the case keeps running on the
#                      one that ships with the repository until you copy it over.
#     CompLaB.xml      what the solver reads
#     postprocess.py   says whether the run is worth believing
#
# Read them before you run them. The point of keeping them here rather than in
# a shared tools directory is that the whole chain for THIS case fits on one
# screen and needs no cross-referencing.

set -euo pipefail

# ---------------------------------------------------------------------------
# 1. PRE-PROCESS   ->  input/geometry.dat
# ---------------------------------------------------------------------------
echo "== pre-processing"
python3 preprocess.py

# ---------------------------------------------------------------------------
# 2. OFFLINE       ->  what has to exist before the solver starts
# ---------------------------------------------------------------------------
echo
# [v1.3] This step used to run unconditionally. offline.sh sweeps the linear
# program through training/generateTrainingData.py, which exits with an error
# when COBRApy is not installed, and with `set -euo pipefail` above that aborted
# the pipeline before the build. The case was therefore unrunnable without
# COBRApy, contradicting this case's own README, which says of offline.sh "Not
# needed to run the case, a fitted header ships with the repository". The step
# is now opt in: set RETRAIN=1 in the environment, or pass --retrain as the
# first argument, to refit the network.
RETRAIN="${RETRAIN:-0}"
if [ "${1:-}" = "--retrain" ]; then RETRAIN=1; fi

if [ "$RETRAIN" = "1" ]; then
    echo "== offline: sweep the linear program, fit the network, verify it. The fitted"
    echo "   header is not installed; the run below uses the shipped one."
    ./offline.sh
else
    echo "== offline: skipped. Using the fitted network that ships with the case."
    echo "   Set RETRAIN=1 (or pass --retrain) to refit it, which needs COBRApy."
fi

# ---------------------------------------------------------------------------
# 3. BUILD
# ---------------------------------------------------------------------------
echo
echo "== building"
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# ---------------------------------------------------------------------------
# 4. RUN           ->  output/
# ---------------------------------------------------------------------------
echo
echo "== running. READ THE START-UP LINES: the geometry, the enabled features"
echo "   and each organism's rate path are echoed before the first step. If any"
echo "   of it is not what you meant, stop now rather than in a fortnight."
mkdir -p output
./complab CompLaB.xml 2>&1 | tee output/run.log

# ---------------------------------------------------------------------------
# 5. POST-PROCESS  ->  a verdict on the run
# ---------------------------------------------------------------------------
echo
echo "== checking the run before believing any picture of it"
python3 postprocess.py
