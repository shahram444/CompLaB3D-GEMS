#!/usr/bin/env bash
# B1 - prepare the metabolic model.
#   in : an SBML or BiGG model
#   out: the flat XML the GLPK path reads (COBRApy reads the SBML directly)
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"

MODEL=${MODEL:-$ROOT/models/e_coli_core.xml.gz}

# extractMM.py reads .xml, .mat and .json, not .gz.
if [[ "$MODEL" == *.gz ]]; then
    gunzip -kf "$MODEL"
    MODEL="${MODEL%.gz}"
fi
OBJ=${OBJ:-Biomass_Ecoli_core}
OUT=${OUT:-model_mm.xml}

python3 "$ROOT/tools/extractMM.py" "$MODEL" --objective "$OBJ" -o "$OUT" -f

# extractMM.py prints the exchange-reaction table as it converts. Fill
# <exchange_reaction_names> from it: a name either resolves or stops the run,
# whereas a positional index that has moved because the model was revised does
# neither, and the run does not complain.
#
# If the objective name is wrong, extractMM.py lists every reaction whose id
# contains "biomass" -- which is the fastest way to find the right one for a
# model you have not used before.
