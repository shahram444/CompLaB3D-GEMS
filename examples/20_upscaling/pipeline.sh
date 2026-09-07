#!/usr/bin/env bash
#
# 20_upscaling  --  THE WHOLE PIPELINE
#
# One resolved aggregate, reduced to the effectiveness factor and Thiele modulus
# a continuum model needs in its place.
#
# Everything this case needs is in this directory. Run it from an assembled case:
#
#     ./scripts/setup_case.sh 20_upscaling run/mycase
#     cd run/mycase
#     ./pipeline.sh
#
# Each step is also runnable on its own, and each is a file you can read:
#
#     preprocess.py           builds the aggregate and the box around it
#     input/aom.thm           the energetics: what makes the core stop
#     kinetics/defineKinetics.hh   the ungated dual-Monod rate law
#     CompLaB.xml             what the solver reads
#     postprocess.py          says whether the number is worth quoting
#     offline/upscale.py      the sweep: one run is one point, this is the curve
#
# EXPECT ROUGHLY FOUR MINUTES for the single run on one core. The sweep at the
# end is optional and takes considerably longer -- it is a separate command, not
# part of this script, so that this script stays something you can run and watch.

set -euo pipefail

# ---------------------------------------------------------------------------
# 1. PRE-PROCESS   ->  input/geometry.dat
# ---------------------------------------------------------------------------
echo "== pre-processing: the aggregate, and the wall that closes the box"
python3 preprocess.py

# ---------------------------------------------------------------------------
# 2. BUILD
# ---------------------------------------------------------------------------
echo
echo "== building"
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# ---------------------------------------------------------------------------
# 3. RUN           ->  output/upscaling.csv
# ---------------------------------------------------------------------------
# READ THE [UPSCALE] LINES AT THE END. They carry the whole result, and the
# steadiness verdict that says whether to believe it.
echo
echo "== running"
mkdir -p output
./complab CompLaB.xml 2>&1 | tee output/run.log

# ---------------------------------------------------------------------------
# 4. POST-PROCESS  ->  a verdict
# ---------------------------------------------------------------------------
echo
echo "== checking the run before believing any number in it"
python3 postprocess.py

# ---------------------------------------------------------------------------
# 5. THE SWEEP, WHICH IS THE ACTUAL POINT
# ---------------------------------------------------------------------------
# Not run here, deliberately. One run gives one (phi, F_T, eta). A continuum
# model needs the curve, and that means one solver run per point.
cat <<'NEXT'

== next: one run is one point

   The run above measured the effectiveness factor for ONE aggregate size at ONE
   bulk composition. A continuum model will be asked for a rate at sizes and
   compositions this run did not visit, so it needs the curve:

       python3 offline/upscale.py --radii 3,4,5 --bulk 1e-5,3e-4,1e-3

   That is nine solver runs. Each is sized by its own diffusion time, R^2/D, so
   the larger radii take substantially longer -- doubling R quadruples the run.
   Start with --radii 3,4 while you are checking it does what you expect.

   It writes upscale_sweep.csv and fits the correction that the classical Thiele
   curve needs once the reaction has an energy limit.
NEXT
