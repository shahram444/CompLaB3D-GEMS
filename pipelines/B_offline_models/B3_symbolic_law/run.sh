#!/usr/bin/env bash
# B3 - search for a rate law.
#   in : a table, one column per variable plus one target column
#   out: a .sym text file, read at start-up, no rebuild needed
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"

DATA=${DATA:-samples.csv}
TARGET=${TARGET:-growth}
VARS=${VARS:-acetate}
POP=${POP:-600}
GENS=${GENS:-60}
OUT=${OUT:-discovered.sym}

# The chemistry, so the file that comes out is complete and runs as it stands.
# REACTION is the balanced reaction, negative for consumed. For a law that
# belongs to an organism, give YIELD and BIOMASS as well and the rate is a
# specific growth rate. For a reaction with no organism, leave them empty and
# set RATE_NAME=extent, which is how fast the reaction itself turns.
#
# Leave REACTION empty and the search writes only the fitted line, which is the
# older behaviour: you then write the other species by hand.
REACTION=${REACTION:-}
YIELD=${YIELD:-}
BIOMASS=${BIOMASS:-}
RATE_NAME=${RATE_NAME:-$TARGET}
UNITS=${UNITS:-per_hour}

ARGS=(--data "$DATA" --target "$TARGET" --inputs "$VARS"
      --units "$UNITS" --rate-name "$RATE_NAME"
      --pop "$POP" --gens "$GENS" --depth 6 --seed 1 --out "$OUT")
[ -n "$REACTION" ] && ARGS+=(--reaction "$REACTION")
[ -n "$YIELD" ]    && ARGS+=(--yield "$YIELD")
[ -n "$BIOMASS" ]  && ARGS+=(--biomass "$BIOMASS")

python "$ROOT/tools/fit_symbolic.py" "${ARGS[@]}"

cat <<'EOF'

What came back is a list, one formula per length, not one answer. Length is how
many pieces a formula has: a number, a name or an operator is one piece. Read
down the list and take the last length where the error still fell by something
worth having. --pick <length> writes that one instead of the automatic choice.

If REACTION was given, the file is complete and runs as it stands: the solver
derives one substrate line per species from the reaction at start-up and prints
them in the log. If it was not, the file holds only the fitted line and you
write the other species yourself, as multiples of it, so the ratios cannot
drift.

Either way, the range lines are already there, one per variable, and they are
enforced at every evaluation rather than advisory.

If you would rather answer questions than remember flags:
    python "$ROOT/tools/make_rate_law.py"
EOF
