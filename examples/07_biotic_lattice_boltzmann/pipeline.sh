#!/usr/bin/env bash
#
# 07_biotic_lattice_boltzmann  --  THE WHOLE PIPELINE
#
# The same case again, with the biomass moved by the lattice Boltzmann solver
# instead of the cellular automaton of 05 or the finite difference of 06.
#
# [v1.3] Earlier text here said "with the population planktonic". That is wrong:
# CompLaB.xml gives this microbe an entry in <material_numbers>, which makes it an
# ATTACHED BIOFILM, meaning a population seeded on its own material number. What
# changes between 05, 06 and 07 is the biomass solver, not whether the population
# is attached.
#
# Everything this case needs is in this directory. Run it from an assembled
# case:
#
#     ./scripts/setup_case.sh 07_biotic_lattice_boltzmann run/mycase
#     cd run/mycase
#     ./pipeline.sh
#
# Each step is also runnable on its own, and each is a file you can read:
#
#     preprocess.py    builds the pore space
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
# 2. BUILD
# ---------------------------------------------------------------------------
echo
echo "== building"
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# ---------------------------------------------------------------------------
# 3. RUN           ->  output/
# ---------------------------------------------------------------------------
echo
echo "== running. READ THE START-UP LINES: the geometry, the enabled features"
echo "   and each organism's rate path are echoed before the first step. If any"
echo "   of it is not what you meant, stop now rather than in a fortnight."
mkdir -p output
./complab CompLaB.xml 2>&1 | tee output/run.log

# ---------------------------------------------------------------------------
# 4. POST-PROCESS  ->  a verdict on the run
# ---------------------------------------------------------------------------
echo
echo "== checking the run before believing any picture of it"
python3 postprocess.py
