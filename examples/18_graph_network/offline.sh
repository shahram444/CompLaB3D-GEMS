#!/usr/bin/env bash
#
# 18_graph_network  --  OFFLINE STEP
#
# Train your own graph network (OPTIONAL)
#
# Run before the solver. pipeline.sh does it for you.

set -euo pipefail

# The case ships input/aom.gnn, so this is optional. It is how that file was
# made, and how you would make one for your own reaction network.
#
# You supply two CSVs: the stoichiometric matrix, which is given to the
# trainer as STRUCTURE rather than learned, and the samples. The Damkohler
# number enters as an edge input, so one network can span reaction-limited
# and transport-limited regimes.

# Both tables ship with this case. aom_stoich.csv is the stoichiometric matrix
# -- the STRUCTURE the trainer is given rather than asked to learn -- and
# aom_samples.csv is the 400 samples input/aom.gnn was fitted from. Point
# STOICH and SAMPLES at your own to train a network for your own reactions.
STOICH=${STOICH:-training/aom_stoich.csv}
SAMPLES=${SAMPLES:-training/aom_samples.csv}

# epochs is deliberately small so this finishes in about a minute. The shipped
# input/aom.gnn was trained with the same settings; more epochs buy very little
# here because the network is small and the samples are clean.
EPOCHS=${EPOCHS:-250}

python3 training/train_graphnet.py --stoich "$STOICH" --data "$SAMPLES" \
        --da 1.0 --rounds 2 --width 6 --epochs "$EPOCHS" --units per_hour \
        --out input/aom_retrained.gnn

# The number that matters is NOT the correlation. It is the ratio between
# species rates, which should come out close to the stoichiometry without
# ever having been a training target. If it does not, the graph is wrong.
# The trainer writes its own fit quality into the file's provenance lines.
echo
echo "what the trainer recorded about the fit it just made:"
grep '^provenance' input/aom_retrained.gnn

cat <<'EOF'

Compare the R and the scaled mse above with the shipped network's:

    grep '^provenance' input/aom.gnn

The retrained file is NOT installed. CompLaB.xml still points at input/aom.gnn.
Copy it over only once you are satisfied, and remember that the species names in
the file must match <name_of_substrates> exactly or the run stops at start-up.
EOF
