#!/usr/bin/env bash
#
# 12_mixed_reaction_types  --  THE WHOLE PIPELINE
#
# Flux balance analysis beside compiled kinetics, in one run.
#
# Everything this case needs is in this directory. Run it from an assembled
# case:
#
#     ./scripts/setup_case.sh 12_mixed_reaction_types run/mycase
#     cd run/mycase
#     ./pipeline.sh
#
# Each step is also runnable on its own, and each is a file you can read:
#
#     preprocess.py    builds the pore space
#     offline.sh       export the model for the FBA organism only
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
echo "== offline: export the model for the FBA organism only"
./offline.sh

# ---------------------------------------------------------------------------
# 3. BUILD
# ---------------------------------------------------------------------------
echo
echo "== building"
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DENABLE_GLPK=ON
cmake --build build -j

# ---------------------------------------------------------------------------
# 4. RUN           ->  output/
# ---------------------------------------------------------------------------
echo
echo "== running. READ THE START-UP LINES: the geometry, the enabled features"
echo "   and each organism's rate path are echoed before the first step. If any"
echo "   of it is not what you meant, stop now rather than in a fortnight."
mkdir -p output
./build/complab CompLaB.xml 2>&1 | tee output/run.log

# ---------------------------------------------------------------------------
# 5. POST-PROCESS  ->  a verdict on the run
# ---------------------------------------------------------------------------
echo
echo "== checking the run before believing any picture of it"
python3 postprocess.py
