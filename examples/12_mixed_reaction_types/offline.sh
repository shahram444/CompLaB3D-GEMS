#!/usr/bin/env bash
#
# 12_mixed_reaction_types  --  OFFLINE STEP
#
# Prepare the metabolic model for the ONE population that uses the linear program.
#
# The other population takes its rate from defineKinetics.hh and needs nothing
# prepared at all.
#
# Run before the solver. pipeline.sh does it for you.

set -euo pipefail

# input/toy_model.xml IS ALREADY THE FLAT FORM the solver reads: four reactions,
# small enough to open and check by hand. Nothing has to be done for this case
# to run, so this script does nothing by default and exists to show you the step
# you WILL need with a model of your own.
#
# extractMM.py converts an SBML, MATLAB or JSON genome-scale model into that
# flat form, and prints the exchange-reaction table you fill
# <exchange_reaction_indices> in from. Give it your model:
#
#     MODEL=iAF987.xml OBJECTIVE=BIOMASS_Gm_GS15_core_79p20M ./offline.sh
#
# A gzipped model has to be expanded first: extractMM.py reads .xml, .mat and
# .json, not .gz.

MODEL=${MODEL:-}

if [[ -z "$MODEL" ]]; then
    cat <<'EOF'
Nothing to prepare: this case ships input/toy_model.xml, already in the flat
form the solver reads.

    grep -c '<reaction' input/toy_model.xml     # four reactions, readable in full

To convert a model of your own instead:

    MODEL=yourmodel.xml OBJECTIVE=<biomass reaction id> ./offline.sh

Then put the output name into <model_filename> in CompLaB.xml, without the
.xml, and fill <exchange_reaction_indices> from the table extractMM.py prints.
Better still, use <exchange_reaction_names>: an index that has silently moved
because the model was revised is a whole afternoon.
EOF
    exit 0
fi

case "$MODEL" in
    *.gz) echo "expanding $MODEL"; gunzip -kf "$MODEL"; MODEL="${MODEL%.gz}" ;;
esac

OBJECTIVE=${OBJECTIVE:-}
if [[ -n "$OBJECTIVE" ]]; then
    python3 training/extractMM.py "$MODEL" -o input/converted_model.xml --objective "$OBJECTIVE" -f
else
    python3 training/extractMM.py "$MODEL" -o input/converted_model.xml -f
fi

echo
echo "Wrote input/converted_model.xml. Set <model_filename>converted_model</model_filename>"
echo "in CompLaB.xml, and fill the exchange entries from the table above."
