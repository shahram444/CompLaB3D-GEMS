#!/usr/bin/env bash
# D2 - check the run before believing any picture of it.
#   in : a run folder, holding CompLaB.xml and output/
#   out: massbalance.txt, and a non-zero exit if a declared sum drifted
#
# [v1.3] This script used to pass --dir and --report, neither of which exists, and
# omitted the positional run folder -- so it exited 2 and the mass-balance report
# this stage exists to produce was unreachable. --conserve is what asks for the
# check; name a sum the chemistry cannot change, and NOT one that is fed from a
# Dirichlet boundary, which is supplied from outside the domain and will always
# "drift". Repeat the flag for several sums.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"

DIR=${1:-$ROOT/run/13_precipitation}
TOL=${TOL:-1e-6}
CONSERVE=${CONSERVE:-}             # e.g. CONSERVE="Ca+calcite"

args=("$DIR" --tol "$TOL")
for c in $CONSERVE; do args+=(--conserve "$c"); done

python3 "$ROOT/tools/postprocess.py" "${args[@]}"
