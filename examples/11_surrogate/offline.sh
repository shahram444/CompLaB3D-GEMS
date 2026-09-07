#!/usr/bin/env bash
#
# 11_surrogate  --  OFFLINE STEP
#
# Sweep the linear program, fit the network to what it returned, verify the
# exported header reproduces the trainer, and look at the response surface. The
# fitted header is not installed: copying it over src/surrogateModel.hh is a
# separate step, printed at the end.
#
# Run before the solver. pipeline.sh does it for you.

set -euo pipefail

# THE LONGEST PREPARATION IN THE REPOSITORY, and the only one ending in a
# recompile. It solves a linear program at every point of a grid, so the cost is
# GRID^2 solves. The repository ships a fitted src/surrogateModel.hh, so this
# case runs without any of the below: this is the procedure that made it, shown
# on the bundled E. coli core model, and how you would make one for a model of
# your own. The shipped header itself came from a sweep of the Geobacter
# metallireducens model iAF987, which is not bundled.

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

# WHAT THE SWEEP RECORDS. Not just the objective: the linear program returns a
# whole flux vector, and the exchange fluxes in it are the ones the solver would
# otherwise have to guess with a Monod term. That guess is exact only where the
# swept bound was the binding constraint, and wrong everywhere else.
#
# So every swept exchange is recorded as an output alongside growth, and --also
# adds a product. Pass --growth-only for the old single-output behaviour.
ALSO=${ALSO:-EX_ac_e}

echo
echo "1/4  sweeping the linear program over a ${GRID} x ${GRID} grid (the slow part)"
python3 training/generateTrainingData.py "$MODEL" \
        --objective "$OBJECTIVE" \
        --exchange EX_glc__D_e --range 0.001 10   --log \
        --exchange EX_o2_e     --range 3e-5 0.5   --log \
        --also "$ALSO" \
        --grid "$GRID" -o training_data.csv

echo
echo "2/4  fitting the network to that sweep"
# TWO FILES FROM ONE FIT, holding the same numbers:
#   surrogate_weights.hh   compiled into the solver. Fastest, needs a rebuild.
#   surrogate_weights.srg  read at start-up by <weights_file>. Swapping network
#                          is then an edit to CompLaB.xml, no rebuild -- which
#                          is how a growth-only and a growth-and-fluxes network
#                          get compared on one case.
python3 training/trainSurrogate.py training_data.csv --name surrogate \
        --layers 10 10 10 10 --restarts 5 \
        -o surrogate_weights.hh --srg surrogate_weights.srg

echo
echo "3/4  checking the exported header reproduces the trainer, OUTPUT BY OUTPUT"
# Do this EVERY time. A transcription error between the fitted weights and the
# C++ header is invisible in the numbers and fatal in the results.
python3 training/verifyExport.py surrogate_weights.hh

echo
echo "4/4  looking at the response surface BEFORE trusting it"
# This path enforces no training range at all. Outside the fitted box it returns
# a confident number and no warning, and a large region returning zero growth
# looks exactly like a region that does not. With several outputs, check each:
# growth is a smooth surface and fits easily, while the flux columns are
# piecewise linear with kinks where the binding constraint changes, and a
# network smooths exactly those.
python3 training/inspectSurrogate.py surrogate_weights.hh --eval 9.0 0.45 || true

cat <<'EOF'

The fitted header is surrogate_weights.hh. It is NOT installed, because
installing it changes what the solver computes and requires a rebuild. When you
are satisfied with the response surface:

    cp surrogate_weights.hh surrogateModel.hh
    cmake --build build -j

[v1.3] The destination is surrogateModel.hh at the CASE ROOT, beside
CMakeLists.txt. setup_case.sh lays a second copy down in src/, but that copy is
not the one compiled: src/complab3d_processors_surrogate.hh includes the header
as "../surrogateModel.hh". Installing into src/ therefore updates an unused
duplicate, and the rebuild succeeds and gives bit identical results from the old
network, with no error to tell you the new fit was ignored.

Until then the case runs on the header that ships with the repository.

Or skip the rebuild entirely and use the run-time copy:

    <surrogate>
        <enabled>true</enabled>
        <weights_file>surrogate_weights.srg</weights_file>
    </surrogate>

Both files came from this one fit and hold the same numbers.
tests/test_surrogate_parity.cpp checks that the two evaluators agree on them.
EOF
