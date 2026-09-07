#!/usr/bin/env bash
#
# 19_thermodynamic_gate  --  OFFLINE STEP
#
# Check the energetics before the solver reads them.
#
# Run before the solver. pipeline.sh does it for you.

set -euo pipefail

# Every other rate path in CompLaB3D has an offline stage: a metabolic model to
# fetch, a network to train, an expression to fit. This one has a file of
# measured energetics, and the equivalent question is whether the threshold
# those energetics describe sits anywhere near the concentrations the run will
# actually visit.
#
# A gate that is shut everywhere looks in the output exactly like a reaction
# that never happened; a gate that is open everywhere looks exactly like no
# gate at all. Both are impossible to tell apart after the fact and trivial to
# see beforehand, which is what this step is for.

python3 offline/thermo_curve.py

cat <<'MSG'

Read the two lines above the table before running the solver:

  - where F_T falls through 0.5, and where it reaches zero
  - whether both are inside the composition range the case covers

If the gate never shuts along that path, the run will finish and mean nothing.
Fix input/aom.thm, not the solver.
MSG
