#!/usr/bin/env bash
#
# 16_complete_pipeline  --  OFFLINE STEP
#
# Expand and convert the genome-scale model, so the run does not depend on a
# remote database still serving the same revision.
#
# Run before the solver. pipeline.sh does it for you.

set -euo pipefail

# This case resolves its model through <model_source>bigg:e_coli_core</model_source>,
# which can reach the network. models/e_coli_core.xml.gz ships here so it does
# not have to: converting it now pins the model to a file you hold, and a result
# from six months ago stays reproducible.

MODEL=${MODEL:-models/e_coli_core.xml.gz}
OBJECTIVE=${OBJECTIVE:-Biomass_Ecoli_core}
OUT=${OUT:-input/e_coli_core_mm.xml}

# extractMM.py reads .xml, .mat and .json, not .gz.
if [[ "$MODEL" == *.gz ]]; then
    echo "expanding $MODEL"
    gunzip -kf "$MODEL"
    MODEL="${MODEL%.gz}"
fi

# NOTE ON THE OBJECTIVE NAME. It is Biomass_Ecoli_core in this file, not the
# BIOMASS_Ecoli_core_w_GAM you will see quoted for other releases of the same
# model. extractMM.py lists every reaction whose id contains "biomass" when the
# name you gave is not one of them, which is the fastest way to find the right
# one for a model you have not used before.
python3 training/extractMM.py "$MODEL" -o "$OUT" --objective "$OBJECTIVE" -f

cat <<EOF

Wrote $OUT, and printed the exchange-reaction table above.

To use it instead of the network lookup, replace <model_source> in CompLaB.xml
with a per-microbe <model_filename>$(basename "${OUT%.xml}")</model_filename>,
and set <exchange_reaction_names> from that table.
EOF
