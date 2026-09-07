#!/usr/bin/env bash
#
# Structural checks on the repository itself.
#
# These are not simulation tests — ./tests/run_tests.sh does that.  These check
# the properties the tree is supposed to have and that are easy to break by
# accident: no file duplicated, one reference configuration, every example
# assembling, every internal link resolving.
#
#   ./tests/check_repo.sh
#
# Exit status is the number of failed checks.

set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

FAIL=0
pass() { printf '  ok    %s\n' "$1"; }
fail() { printf '  FAIL  %s\n' "$1"; FAIL=$((FAIL+1)); }

echo "Checking the repository tree"

# --- 1. nothing exists twice OUTSIDE the examples ----------------------------
#
# Inside examples/ duplication is deliberate. An example directory carries
# everything that case needs -- its pore space, its metabolic model, the
# training code that produced what it was trained on -- so that the folder is
# the whole procedure and nothing has to be fetched from elsewhere in the tree.
# Several cases legitimately hold the same geometry or the same trainer.
#
# THE COST OF THAT, stated plainly: a fix to a trainer has to be applied to
# every example that carries it. `make sync-examples` is not provided on
# purpose -- see CONTRIBUTING.md, which lists what to update together.
#
# Everywhere else the rule still holds: one copy, and no second one.
DUPES="$(find . -path ./run -prune -o -path ./.git -prune -o -path ./examples -prune -o \
         -path ./build -prune -o -name '*.o' -prune -o -name '*.a' -prune -o \
         -type f -size +200c -print \
         | xargs -r md5sum | sort | awk '{c[$1]=c[$1]" "$2; n[$1]++} END {for (k in n) if (n[k]>1) print c[k]}')"
if [[ -z "$DUPES" ]]; then
    pass "no two identical files over 200 bytes outside examples/"
else
    fail "identical files found:"
    echo "$DUPES" | sed 's/^/          /'
fi

# --- 1b. every copy inside examples/ still MATCHES its source ----------------
# Duplication is allowed; drift is not. A trainer or a geometry that has been
# fixed in one place and not the others is exactly what this catches.
DRIFT=""
check_copy() {   # check_copy <source> <path under examples/*/>
    local src="$1" rel="$2" f
    [[ -f "$src" ]] || return 0
    for f in examples/*/"$rel"; do
        [[ -f "$f" ]] || continue
        cmp -s "$src" "$f" || DRIFT="$DRIFT $f"
    done
}
check_copy tools/extractMM.py                          training/extractMM.py
check_copy tools/fit_symbolic.py                       training/fit_symbolic.py
check_copy tools/train_graphnet.py                     training/train_graphnet.py
check_copy tools/makeEquilibrium.py                    training/makeEquilibrium.py
check_copy tools/makeKinetics.py                       training/makeKinetics.py
check_copy tools/complab3d_cobrapy.py                  training/complab3d_cobrapy.py
check_copy tests/xval_gnn.py                           training/xval_gnn.py
check_copy tools/upscale_sweep.py                      offline/upscale.py
check_copy tools/surrogate/generateTrainingData.py     training/generateTrainingData.py
check_copy tools/surrogate/trainSurrogate.py           training/trainSurrogate.py
check_copy tools/surrogate/verifyExport.py             training/verifyExport.py
check_copy tools/surrogate/inspectSurrogate.py         training/inspectSurrogate.py
check_copy models/e_coli_core.xml.gz                   models/e_coli_core.xml.gz
check_copy models/toy_model.xml                        input/toy_model.xml
[[ -z "$DRIFT" ]] && pass "every copy inside examples/ matches its source" \
                  || fail "copies have drifted from their source:$DRIFT"

# --- 2. exactly one annotated reference configuration ------------------------
N=$(find . -name '*reference*.xml' -not -path './run/*' | wc -l)
[[ "$N" == 1 ]] && pass "exactly one reference configuration" \
                || fail "expected 1 reference configuration, found $N"

# --- 3. the reference configuration parses, and covers all six rate paths -----
if python3 -c "import xml.etree.ElementTree as E; E.parse('config/CompLaB.reference.xml')" 2>/dev/null; then
    MISSING=""
    for p in kinetics fba cobrapy surrogate symbolic graphnet thermodynamics; do
        grep -q "$p" config/CompLaB.reference.xml || MISSING="$MISSING $p"
    done
    [[ -z "$MISSING" ]] && pass "reference configuration parses, names all six rate paths and the gate" \
                        || fail "reference configuration does not mention:$MISSING"
else
    fail "config/CompLaB.reference.xml does not parse"
fi

# --- 2b. the everything file lists every tag the solver actually reads --------
#   config/CompLaB.everything.xml claims to have one line per ability. A claim
#   like that decays the first time someone adds a tag and forgets. This walks
#   the source for every XML path the solver opens and fails on any that the
#   file does not mention, so the claim stays true by construction rather than
#   by diligence.
UNLISTED=$(python3 tests/check_everything_xml.py)
if python3 -c "import xml.etree.ElementTree as E; E.parse('config/CompLaB.everything.xml')" 2>/dev/null; then
    [[ -z "$UNLISTED" ]] && pass "the everything file lists every tag the solver reads" \
                         || fail "read by the solver but missing from CompLaB.everything.xml:$UNLISTED"
else
    fail "config/CompLaB.everything.xml does not parse"
fi

# --- 3b. every rate path is actually CALLED from the solver -------------------
#   This is the check that was missing.  Documenting a path in the reference XML,
#   implementing its processor and unit-testing that processor all passed while
#   nothing in complab.cpp ever constructed it: the path was correct and dead.
#   A processor that is never instantiated is not a feature.
DEAD=""
for proc in run_kinetics runFBA_glpk3D runFBA_cobrapy3D run_surrogate3D \
            run_symbolic3D run_graphnet3D run_symbolic_abiotic3D run_graphnet_abiotic3D; do
    grep -q "new ${proc}<" src/complab.cpp || DEAD="$DEAD $proc"
done
[[ -z "$DEAD" ]] && pass "every rate processor is constructed in src/complab.cpp" \
                 || fail "implemented but never called from complab.cpp:$DEAD"

# --- 3b2. the thermodynamic gate is applied in EVERY biotic rate path ---------
#   complab3d_thermo.hh is reachable the moment one processor includes it, so 3c
#   below would pass with the gate wired into one path and missing from five.
#   The claim in the documentation is that the gate multiplies whichever path an
#   organism already uses; this is that claim, checked.
UNGATED=""
for f in complab3d_processors_part1.hh complab3d_processors_fba.hh \
         complab3d_processors_surrogate.hh complab3d_processors_symbolic.hh \
         complab3d_processors_graphnet.hh; do
    grep -q "complab_thermo::gateFor" "src/$f" || UNGATED="$UNGATED $f"
done
#   complab3d_processors_fba.hh holds two paths, so one grep is not enough there.
N=$(grep -c "complab_thermo::gateFor" src/complab3d_processors_fba.hh || true)
[[ "$N" -ge 2 ]] || UNGATED="$UNGATED complab3d_processors_fba.hh(only-one-of-two)"
[[ -z "$UNGATED" ]] && pass "the thermodynamic gate is applied in every biotic rate path" \
                    || fail "rate paths that never consult the gate:$UNGATED"

# --- 3b3. the two FBA back ends get separate organism lists -------------------
#   A GLPK organism and a COBRApy organism may share a run. What makes that safe
#   is that each processor is handed only its OWN organisms: complab.cpp builds
#   glpk_globalId and cpy_globalId separately, and each FBA processor gets its own
#   lattice vector. Hand both the same list again -- the single usesFBA() list the
#   older code built -- and each back end starts solving for the other's organisms,
#   where cfg.vec_lp[gM] is null for a COBRApy microbe. That is a segfault a
#   thousand iterations in, not a compile error, so it is checked here.
BAD=""
grep -q "glpk_globalId" src/complab.cpp || BAD="$BAD (no glpk_globalId)"
grep -q "cpy_globalId"  src/complab.cpp || BAD="$BAD (no cpy_globalId)"
#   Each dispatch must name its own list and its own lattice vector.
grep -q "runFBA_glpk3D<T,RXNDES>(nx, num_of_substrates, (plint) glpk_globalId.size()" src/complab.cpp \
    || BAD="$BAD (glpk dispatch does not use glpk_globalId)"
grep -q "runFBA_cobrapy3D<T,RXNDES>(nx, num_of_substrates, (plint) cpy_globalId.size()" src/complab.cpp \
    || BAD="$BAD (cobrapy dispatch does not use cpy_globalId)"
grep -q "ptr_glpk_lattices)" src/complab.cpp || BAD="$BAD (glpk lattice vector)"
grep -q "ptr_cpy_lattices)"  src/complab.cpp || BAD="$BAD (cobrapy lattice vector)"
#   And the ban itself must stay gone: mixing the back ends is allowed now.
grep -q "cannot be mixed in one simulation" src/complab3d_metabolic.hh \
    && BAD="$BAD (the mixing ban is back)"
[[ -z "$BAD" ]] && pass "GLPK and COBRApy organisms may share a run, each with its own list" \
                || fail "FBA back-end separation:$BAD"

# --- 3c. no header in src/ is unreachable from the program --------------------
#   The same failure one level up: a header nothing includes is a header nothing
#   compiles, so it can rot without any test noticing.
UNREACH=$(python3 - <<'PYEOF'
import os, re, collections
src = "src"
heads = [f for f in os.listdir(src) if f.endswith(".hh")]
inc = collections.defaultdict(set)
for f in os.listdir(src):
    if not f.endswith((".hh", ".cpp")):
        continue
    t = open(os.path.join(src, f), encoding="utf8", errors="replace").read()
    for m in re.finditer(r'#include\s+"([^"]+)"', t):
        inc[f].add(os.path.basename(m.group(1)))
seen, stack = set(), ["complab.cpp"]
while stack:
    cur = stack.pop()
    if cur in seen:
        continue
    seen.add(cur)
    stack.extend(n for n in inc.get(cur, ()) if n in heads)
# _reserve is a deliberately kept alternative implementation, not dead code
skip = {"complab3d_processors_part4_eqsolver_NR_reserve.hh"}
print(" ".join(sorted(set(heads) - seen - skip)))
PYEOF
)
[[ -z "$UNREACH" ]] && pass "every header in src/ is reachable from complab.cpp" \
                    || fail "unreachable from complab.cpp:$UNREACH"

# --- 4. every example carries what it needs ----------------------------------
BAD=""
for d in examples/[0-9]*; do
    c="$(basename "$d")"
    [[ -f "$d/CompLaB.xml"     ]] || BAD="$BAD $c:CompLaB.xml"
    [[ -f "$d/README.md"       ]] || BAD="$BAD $c:README.md"
    [[ -f "$d/preprocess.py"   ]] || BAD="$BAD $c:preprocess.py"
    [[ -f "$d/postprocess.py"  ]] || BAD="$BAD $c:postprocess.py"
    [[ -f "$d/pipeline.sh"     ]] || BAD="$BAD $c:pipeline.sh"
    [[ -f "$d/input/geometry.dat" ]] || BAD="$BAD $c:input/geometry.dat"
    # nothing in an example may reach back into the shared tree
    if grep -rn "\.\./\.\./tools\|\.\./\.\./models" "$d" --include=*.sh --include=*.py >/dev/null 2>&1; then
        BAD="$BAD $c:(reaches into ../../)"
    fi
done
[[ -z "$BAD" ]] && pass "every example is self-contained" || fail "incomplete:$BAD"

# --- 4b. the shipped geometry is the one the generator writes ----------------
#   Every case ships input/geometry.dat so it runs before you have run anything,
#   and a preprocess.py that rebuilds it. Nothing checked that the two agree, and
#   they had drifted: examples 13, 14 and 15 shipped a pore space with 696 mineral
#   voxels while their own generator wrote 72. Both ran. Which one a published
#   number came from was unknowable.
#
#   This also catches the other half of the same problem: <nx>/<ny>/<nz> in
#   CompLaB.xml must be the dimensions the generator actually wrote, or the solver
#   reads the file into a domain of the wrong shape.
BAD="$(python3 - <<'PYEOF'
import glob, os, re, subprocess, filecmp, shutil, tempfile
bad = []
for d in sorted(glob.glob("examples/[0-9]*")):
    gp = os.path.join(d, "input", "geometry.dat")
    if not os.path.exists(gp):
        bad.append(f"{os.path.basename(d)}: no geometry.dat"); continue
    keep = tempfile.mktemp()
    shutil.copy(gp, keep)
    try:
        out = subprocess.run(["python3", "preprocess.py"], cwd=d,
                             capture_output=True, text=True)
        if out.returncode != 0:
            bad.append(f"{os.path.basename(d)}: preprocess.py failed"); continue
        if not filecmp.cmp(gp, keep, shallow=False):
            bad.append(f"{os.path.basename(d)}: geometry.dat differs from preprocess.py output")
        m = re.search(r"(\d+) x (\d+) x (\d+)", out.stdout)
        x = open(os.path.join(d, "CompLaB.xml")).read()
        dims = re.search(r"<nx>(\d+)</nx>\s*\n\s*<ny>(\d+)</ny>\s*\n\s*<nz>(\d+)</nz>", x)
        if m and dims and m.groups() != dims.groups():
            bad.append(f"{os.path.basename(d)}: CompLaB.xml says {dims.groups()}, "
                       f"the generator writes {m.groups()}")
    finally:
        shutil.copy(keep, gp)
        os.unlink(keep)
print("\n".join(bad))
PYEOF
)"
[[ -z "$BAD" ]] && pass "every shipped geometry is the one its generator writes" \
               || { fail "geometry out of step with its generator:"; echo "$BAD" | sed 's/^/          /'; }

# --- 5. every example assembles ----------------------------------------------
TMP="$(mktemp -d)"
BAD=""
for d in examples/[0-9]*; do
    c="$(basename "$d")"
    if ./scripts/setup_case.sh "$c" "$TMP/$c" >/dev/null 2>&1; then
        for need in CompLaB.xml defineKinetics.hh defineAbioticKinetics.hh input/geometry.dat; do
            [[ -f "$TMP/$c/$need" ]] || BAD="$BAD $c:$need"
        done
    else
        BAD="$BAD $c:(setup_case.sh failed)"
    fi
done
rm -rf "$TMP"
[[ -z "$BAD" ]] && pass "all $(ls -d examples/[0-9]* | wc -l) examples assemble and are complete" \
                || fail "incomplete:$BAD"

# --- 6. every shell script parses, and is executable -------------------------
BAD=""
while read -r s; do
    bash -n "$s" 2>/dev/null || BAD="$BAD $s(syntax)"
    [[ -x "$s" ]]            || BAD="$BAD $s(not executable)"
done < <(find . -name '*.sh' -not -path './run/*')
[[ -z "$BAD" ]] && pass "every shell script parses and is executable" || fail "$BAD"

# --- 7. every internal markdown link resolves --------------------------------
BAD="$(python3 - <<'PY'
import re, os, pathlib
bad = []
for md in pathlib.Path('.').rglob('*.md'):
    if 'run/' in str(md): continue
    for target in re.findall(r'\]\(([^)#]+)\)', md.read_text(encoding='utf-8', errors='replace')):
        if target.startswith(('http', 'mailto:', '#')): continue
        if not (md.parent / target).exists():
            bad.append(f'{md}: {target}')
print('\n'.join(bad))
PY
)"
[[ -z "$BAD" ]] && pass "every internal markdown link resolves" \
               || { fail "broken links:"; echo "$BAD" | sed 's/^/          /'; }

# --- 8. every pipeline stage documents itself --------------------------------
BAD=""
for d in pipelines/*/ pipelines/*/*/; do
    [[ -f "$d/README.md" ]] || BAD="$BAD $d"
done
[[ -z "$BAD" ]] && pass "every pipeline stage has a README" || fail "no README in:$BAD"

echo
if [[ "$FAIL" == 0 ]]; then
    echo "All structural checks passed."
else
    echo "$FAIL check(s) failed."
fi
exit "$FAIL"
