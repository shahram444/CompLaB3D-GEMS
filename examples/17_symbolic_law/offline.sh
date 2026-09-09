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

# pop and gens are deliberately small so this finishes in under a minute, and
# that is enough to see the machinery work but NOT enough to find the law.
#
# What the small search returns on this data is the bilinear product
# 647 x acetate x o2, at about 4% error. That is not a failure: over the
# concentrations in the samples, both Monod terms are in their linear part and
# the product IS the right answer to two significant figures. Recovering the
# saturation, and with it the two half-saturation constants, needs a longer
# search, which is what the values in the comment below are for:
#
#     POP=600 GENS=60 ./offline.sh
#
# At that size the list reaches expressions that carry 0.010 and 0.050 inside
# them, which are the true oxygen and acetate half-saturation constants, found
# without the search ever being told the law is Monod.
POP=${POP:-200}
GENS=${GENS:-15}

# --verify runs the whole search a second time and says whether the two agree.
# It doubles the run time and it is worth it before a number goes in a paper:
# the constant fitting sits on scipy's optimiser, and on an over-parameterised
# expression that optimiser is not bit-reproducible even on one thread, so the
# only honest claim is one that has been checked on the machine making it.
VERIFY=${VERIFY:-1}
[ "$VERIFY" = "1" ] && VERIFY_FLAG=--verify || VERIFY_FLAG=

python3 training/fit_symbolic.py --data "$DATA" --target growth \
        --inputs acetate,o2 --units per_hour \
        --pop "$POP" --gens "$GENS" --depth 6 --seed 1 $VERIFY_FLAG \
        --out input/growth_discovered.sym

cat <<'EOF'

The search wrote input/growth_discovered.sym. It is NOT installed: compare it
with the shipped input/growth.sym first, and copy it over only when you are
satisfied with the trade-off you picked off the Pareto front.

READ THE LIST, DO NOT TRUST THE AUTOMATIC PICK. It takes the steepest gain per
node, which on this data lands on the crude bilinear expression and skips the
one that recovers the half-saturation constants. --pick <nodes> takes any other
row; the node counts differ from run to run, so read them off the list you just
got rather than one from yesterday.

The file carries its own provenance: the command, the seed, the thread setting
and the library versions are written into its header, because a fitted law whose
command has been lost is not a result.
EOF
