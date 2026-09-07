#!/usr/bin/env bash
#
# 22_cybernetic_switching  --  THE WHOLE PIPELINE
#
# Growth as a competition between one option per carbon source, blended by how much of the limiting element each brings in.
#
# Everything this case needs is in this directory. Run it from an assembled
# case:
#
#     ./scripts/setup_case.sh 22_cybernetic_switching run/mycase
#     cd run/mycase
#     ./pipeline.sh
#
# Each step is also runnable on its own:
#
#     preprocess.py    builds the pore space
#     CompLaB.xml      what the solver reads
#     postprocess.py   says whether the run is worth believing
#
set -euo pipefail

# ---------------------------------------------------------------------------
# 1. PRE-PROCESS   ->  input/geometry.dat
# ---------------------------------------------------------------------------
echo "== pre-processing"
python3 preprocess.py

# ---------------------------------------------------------------------------
# 2. BUILD
#    This case solves a linear program in every voxel, so GLPK is required.
#    The CMake option is OFF by default; without it the build succeeds and the
#    solver then stops at start-up saying it was built without GLPK.
# ---------------------------------------------------------------------------
echo "== building"
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DENABLE_GLPK=ON
cmake --build build -j"$(nproc 2>/dev/null || echo 4)"

# ---------------------------------------------------------------------------
# 3. RUN
# ---------------------------------------------------------------------------
echo "== running"
mkdir -p output
./complab CompLaB.xml 2>&1 | tee output/run.log

# ---------------------------------------------------------------------------
# 4. POST-PROCESS
# ---------------------------------------------------------------------------
echo "== post-processing"
python3 postprocess.py
