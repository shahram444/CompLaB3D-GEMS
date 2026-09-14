#!/usr/bin/env bash
#
# 19_thermodynamic_gate  --  THE WHOLE PIPELINE
#
# A reaction rate multiplied by the free energy available for it, per voxel.
#
# Everything this case needs is in this directory. Run it from an assembled
# case:
#
#     ./scripts/setup_case.sh 19_thermodynamic_gate run/mycase
#     cd run/mycase
#     ./pipeline.sh
#
# Each step is also runnable on its own, and each is a file you can read:
#
#     preprocess.py           builds the aggregate
#     offline.sh              checks the energetics before the solver reads them
#     input/aom.thm           the energetics themselves
#     defineKinetics.hh       the ungated rate law
#     CompLaB.xml             what the solver reads
#     postprocess.py          says whether the gate did anything
#
# Read them before you run them.

set -euo pipefail

# ---------------------------------------------------------------------------
# 1. PRE-PROCESS   ->  input/geometry.dat
# ---------------------------------------------------------------------------
echo "== pre-processing"
python3 preprocess.py

# ---------------------------------------------------------------------------
# 2. OFFLINE       ->  what has to be right before the solver starts
# ---------------------------------------------------------------------------
echo
echo "== offline: the gate these energetics describe"
./offline.sh

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
echo "== running. READ THE [THM] START-UP LINES: the temperature, the reaction as"
echo "   parsed, dG0 and the energy threshold are echoed before the first step."
mkdir -p output
./complab CompLaB.xml 2>&1 | tee output/run.log

# ---------------------------------------------------------------------------
# 5. RUN AGAIN, UNGATED  ->  the only honest comparison
# ---------------------------------------------------------------------------
# A gated run on its own says nothing about what the gate did. The same case
# with the gate switched off and NOTHING else changed is the control.
echo
echo "== running again with the gate off, as a control"
# The sed range is anchored to the tag at the start of a line, because the
# comment at the top of CompLaB.xml also contains the word <thermodynamics>,
# and an unanchored range would run from that comment down through the
# <diagnostics> block and switch THAT off instead. The control would then write
# no summary.csv, the comparison below would silently read the gated run's file,
# and the two runs would look identical for a reason that has nothing to do
# with the gate.
sed '/^    <thermodynamics>$/,/^    <\/thermodynamics>$/ s|<enabled>true</enabled>|<enabled>false</enabled>|' \
    CompLaB.xml > nogate.xml
diff CompLaB.xml nogate.xml || true          # exactly one line must differ

mv -f output/summary.csv output/summary_gated.csv 2>/dev/null || true
./complab nogate.xml 2>&1 | tee output/run_nogate.log
mv -f output/summary.csv output/summary_nogate.csv 2>/dev/null || true
mv -f output/summary_gated.csv output/summary.csv 2>/dev/null || true

# ---------------------------------------------------------------------------
# 6. POST-PROCESS  ->  a verdict on the run
# ---------------------------------------------------------------------------
echo
echo "== checking the run before believing any picture of it"
python3 postprocess.py
