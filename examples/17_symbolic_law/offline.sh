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
# What comes back is a LIST, one formula per length, not one answer. Length is
# how many pieces a formula has, where a number, a name or an operator is one
# piece. Read down the list and take the last formula where the error still
# fell by something worth having.
#
# The file the search writes is complete and runs as it stands. It carries the
# fitted rate line, the reaction, the yield and the biomass name, one range
# line per variable, and its own provenance header. Nothing is left for you to
# finish by hand.
#
# If you would rather answer questions than remember flags:
#
#     python3 training/make_rate_law.py
#
# reads any training table, asks what its columns mean, and runs the same
# search. It prints the equivalent command at the end so a run can be repeated.

# training/growth_samples.csv ships with this case: 400 points drawn from the
# same dual-Monod law input/growth.sym describes, with 2% noise on the growth
# column so the search has something realistic to work on. Point DATA at your
# own table to fit your own law.
#
# That table is not a fixed asset. training/make_training_data.py writes it, from
# a law you can read, and writes the law beside it as growth_samples_truth.json.
# That is what lets this case ask whether the search found the RIGHT law rather
# than merely a law that fits:
#
#     python3 training/make_training_data.py --n 2000 --noise 0.05
#
# It also prints how far the sampling reaches past each half-saturation constant,
# and warns if it stops short. That ratio decides whether the constants can be
# recovered at all: below its constant a Monod term is a straight line, and every
# value of the constant fits a straight line equally well.
DATA=${DATA:-training/growth_samples.csv}

# pop and gens are deliberately small so this finishes in under a minute, and
# that is enough to see the machinery work but NOT enough to find the law.
#
# The small search stops around 26% error and every formula on its list is
# visibly the wrong shape. That is the honest signal, and it is worth saying
# that it was not always so: on the narrower training table this case used to
# ship, the small search reached 4% with a plain product of the two
# concentrations, because the data never left the straight-line part of either
# Monod term. A short run looked like a success. It is better that it fails
# where it should.
#
# Finding the law needs a longer search:
#
#     POP=600 GENS=60 ./offline.sh
#
# Expect that to take minutes rather than seconds, and longer than it did on the
# narrow table: fitting the constants of a curve is real work, whereas fitting
# the slope of a straight line is not.
POP=${POP:-200}
GENS=${GENS:-15}

# --verify runs the whole search a second time and says whether the two agree.
# It doubles the run time and it is worth it before a number goes in a paper:
# the constant fitting sits on scipy's optimiser, and on an over-parameterised
# formula that optimiser is not bit-reproducible even on one thread, so the
# only honest claim is one that has been checked on the machine making it.
#
# It compares what the formulas COMPUTE, not how they are spelled, and it
# judges the gap between the two runs against each formula's own error against
# the data. Expect the short and middle lengths to agree and the longest not
# to: a long formula carries more numbers than 400 rows can pin down, so there
# is no single best answer for it, and those are the same formulas nobody
# should be quoting.
VERIFY=${VERIFY:-1}
[ "$VERIFY" = "1" ] && VERIFY_FLAG=--verify || VERIFY_FLAG=

# --reaction, --yield and --biomass are what make the written file complete.
# The solver derives one substrate line per species from them at start-up:
#
#     acetate:  (-1 / 1) / 0.4  =  -2.5   x growth x Bug
#     o2:       (-2 / 1) / 0.4  =  -5.0   x growth x Bug
#
# so the ratio between the two is arithmetic on the reaction rather than two
# numbers typed by hand that nothing checks. For a reaction with no organism,
# drop --yield and --biomass and pass --rate-name extent instead.
python3 training/fit_symbolic.py --data "$DATA" --target growth \
        --inputs acetate,o2 --units per_hour \
        --reaction "acetate -1  o2 -2" --yield "acetate 0.4" --biomass Bug \
        --pop "$POP" --gens "$GENS" --depth 6 --seed 1 $VERIFY_FLAG \
        --out input/growth_discovered.sym

cat <<'EOF'

The search wrote input/growth_discovered.sym. It is NOT installed: compare it
with the shipped input/growth.sym first, and copy it over only when you are
satisfied with the length you picked off the list.

READ THE LIST, DO NOT TRUST THE AUTOMATIC PICK. It takes the biggest drop in
error per extra piece, which on this data lands on the crude bilinear formula
and skips the one that recovers the half-saturation constants. --pick <length>
takes any other row; the lengths differ from run to run, so read them off the
list you just got rather than one from yesterday.

The file carries its own provenance: the command, the seed, the thread setting
and the library versions are written into its header, because a fitted law whose
command has been lost is not a result.
EOF
