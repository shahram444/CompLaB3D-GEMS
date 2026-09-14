#!/usr/bin/env bash
#
# 10_fba_cobrapy  --  OFFLINE STEP
#
# Check that COBRApy is importable from the same Python the embedded
# interpreter will load, and that the model file parses.
#
# Run before the solver. pipeline.sh does it for you.

set -euo pipefail

# There is no export step on this path: the C++ side reads input/toy_model.xml
# in the flat <Metabolic_Model> format and hands the stoichiometry to Python as
# plain arrays, so COBRApy never opens the file itself. What IS worth ruling out
# before a run is the failure that costs the most time: the embedded interpreter
# finding a different Python from the one you installed cobra into. The solver
# prints its interpreter path in the first few lines; this prints the shell's,
# so a disagreement between the two is visible immediately.

python3 -c "import cobra, sys; print('cobra', cobra.__version__, 'in', sys.executable)"

# And that the model file is the flat format the solver expects, with a
# stoichiometric matrix whose length is metabolites x reactions.
python3 - <<'PY'
import xml.etree.ElementTree as ET
r = ET.parse('input/toy_model.xml').getroot()
nm = int(r.findtext('nmet'))
nr = int(r.findtext('nrxn'))
S = r.findtext('S').split()
assert len(S) == nm * nr, "S has %d entries, expected %d x %d" % (len(S), nm, nr)
print("model ok: %d metabolites, %d reactions, S is %d entries" % (nm, nr, len(S)))
PY
