#!/usr/bin/env bash
# C - assemble a case, build it, run it.
#   in : everything A and B produced
#   out: VTI fields, the log, the reports
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

CASE=${1:-13_precipitation}
DEST=${2:-$ROOT/run/$CASE}
NP=${NP:-1}

"$ROOT/scripts/setup_case.sh" "$CASE" "$DEST"

cd "$DEST"
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

if [[ "$NP" -gt 1 ]]; then
    mpirun -np "$NP" ./complab CompLaB.xml
else
    ./complab CompLaB.xml
fi

echo
echo "Read the start-up lines above.  If the geometry, the enabled features or"
echo "the rate path is not what you meant, stop now rather than after the run."
