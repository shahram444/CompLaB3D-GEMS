#!/usr/bin/env bash
#
# 10_fba_cobrapy  --  OFFLINE STEP
#
# Check the model loads in the interpreter COBRApy will use
#
# Run before the solver. pipeline.sh does it for you.

set -euo pipefail

# This path reads the SBML directly, so there is no export step. What there IS
# is a failure mode worth ruling out first: the embedded interpreter finding a
# different Python from the one you installed cobra into. The solver prints the
# interpreter path in its first few lines; this checks the model in the shell's
# Python, so a disagreement between the two is visible immediately.

python3 -c "import cobra, sys; print('cobra', cobra.__version__, 'in', sys.executable)"
python3 -c "import cobra; m = cobra.io.read_sbml_model('input/toy_model.xml'); print(m.optimize())"
