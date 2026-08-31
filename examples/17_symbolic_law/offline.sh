#!/usr/bin/env bash
#
# 17_symbolic_law  --  OFFLINE STEP
#
# Search for your own rate law (OPTIONAL)
#
# Run before the solver. pipeline.sh does it for you.

set -euo pipefail

# The case ships input/growth.sym, so this is optional. It is how you would
# find a law from data of your own.
#
# What comes back is a PARETO SET, one expression per node count, not one
# answer. Take the elbow: the shortest expression whose error is acceptable.
# Then finish the file by hand -- write the other species as multiples of the
# fitted rate so the stoichiometry cannot drift, and add one range line per
# variable, because ranges are enforced rather than advisory.

# training/growth_samples.csv ships with this case: 400 points drawn from the
# same dual-Monod law input/growth.sym describes, with 2% noise on the growth
# column so the search has something realistic to work on. Point DATA at your
# own table to fit your own law.
DATA=${DATA:-training/growth_samples.csv}

# pop and gens are deliberately small so this finishes in under a minute. A
# real search wants --pop 600 --gens 60 or more; the Pareto front barely moves
# after that, but it does move.
POP=${POP:-200}
GENS=${GENS:-15}

python3 training/fit_symbolic.py --data "$DATA" --target growth \
        --inputs acetate,o2 --units per_hour \
        --pop "$POP" --gens "$GENS" --depth 6 --seed 1 \
        --out input/growth_discovered.sym

cat <<'EOF'

The search wrote input/growth_discovered.sym. It is NOT installed: compare it
with the shipped input/growth.sym first, and copy it over only when you are
satisfied with the trade-off you picked off the Pareto front.
EOF
