#!/usr/bin/env bash
#
# Run the whole CompLB3D example suite and report pass or fail per case.
#
# Part of the CompLaB program.  GNU Affero General Public License v3 or later.
#
# ---------------------------------------------------------------------------
# WHY THIS REBUILDS BETWEEN CASES
#
# CompLaB compiles its rate laws in. defineKinetics.hh and
# defineAbioticKinetics.hh are #included, not read at run time, so two cases
# with different chemistry need different binaries. This script copies each
# case's headers into the source tree and rebuilds. That is not this script
# being clumsy; it is how the code is built, and the 2D suite worked the same
# way.
#
# Three of the cases also need a different cmake configuration: 09 and 12 need
# GLPK, 10 needs Python. The script groups them so the number of full rebuilds
# stays down.
#
# ---------------------------------------------------------------------------
# USAGE
#
#   ./runAllExamples.sh <path to CompLB3D>            everything it can build
#   ./runAllExamples.sh <path to CompLB3D> 01 05 13   only those cases
#
#   SKIP_BUILD=1 ./runAllExamples.sh ...   reuse the existing binary, for when
#                                          you are iterating on one case
#
# It restores your original defineKinetics.hh and defineAbioticKinetics.hh on
# exit, including on Ctrl-C. Check that it did before you commit anything.
# ---------------------------------------------------------------------------

set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE="${1:-}"
if [ -z "$CODE" ] || [ ! -f "$CODE/CMakeLists.txt" ]; then
    echo "usage: $0 <path to CompLB3D> [case numbers]"
    echo "       the path must contain CMakeLists.txt"
    exit 2
fi
CODE="$(cd "$CODE" && pwd)"
shift || true
WANT=("$@")

RESULTS="$HERE/results"
mkdir -p "$RESULTS"

# ---------------------------------------------------------------------------
# Put the user's own kinetics headers back, whatever happens.
# ---------------------------------------------------------------------------
BACKUP="$(mktemp -d)"
for f in defineKinetics.hh defineAbioticKinetics.hh; do
    [ -f "$CODE/$f" ] && cp "$CODE/$f" "$BACKUP/$f"
done
restore() {
    for f in defineKinetics.hh defineAbioticKinetics.hh; do
        [ -f "$BACKUP/$f" ] && cp "$BACKUP/$f" "$CODE/$f"
    done
    rm -rf "$BACKUP"
}
trap restore EXIT INT TERM

# ---------------------------------------------------------------------------
# Which cmake flags each case needs. Cases sharing a configuration are run
# together so the tree is reconfigured three times, not fifteen.
# ---------------------------------------------------------------------------
cmake_flags_for() {
    case "$1" in
        09|12) echo "-DENABLE_GLPK=ON" ;;
        10)    echo "-DENABLE_COBRAPY=ON" ;;
        *)     echo "" ;;
    esac
}

CURRENT_FLAGS="__none__"
configure_if_needed() {
    local flags="$1"
    [ "$flags" = "$CURRENT_FLAGS" ] && return 0
    echo "    configuring: cmake $flags .."
    mkdir -p "$CODE/build"
    ( cd "$CODE/build" && cmake $flags .. ) > "$RESULTS/cmake.log" 2>&1 || {
        echo "    cmake FAILED, see $RESULTS/cmake.log"
        return 1
    }
    CURRENT_FLAGS="$flags"
    return 0
}

PASS=0; FAIL=0; SKIP=0
FAILED_CASES=""

for dir in "$HERE"/[0-9][0-9]_*; do
    [ -d "$dir" ] || continue
    name="$(basename "$dir")"
    num="${name:0:2}"

    if [ ${#WANT[@]} -gt 0 ]; then
        found=0
        for w in "${WANT[@]}"; do [ "$w" = "$num" ] && found=1; done
        [ $found -eq 1 ] || continue
    fi

    echo ""
    echo "==================================================================="
    echo "  $name"
    echo "==================================================================="

    flags="$(cmake_flags_for "$num")"

    if [ "${SKIP_BUILD:-0}" != "1" ]; then
        cp "$dir/defineKinetics.hh"        "$CODE/defineKinetics.hh"
        cp "$dir/defineAbioticKinetics.hh" "$CODE/defineAbioticKinetics.hh"

        if ! configure_if_needed "$flags"; then
            echo "  SKIP: this configuration could not be prepared."
            [ -n "$flags" ] && echo "        $flags usually means a missing optional dependency."
            SKIP=$((SKIP+1)); continue
        fi
        echo "    building ..."
        if ! ( cd "$CODE/build" && make -j4 ) > "$RESULTS/$name.build.log" 2>&1; then
            echo "  FAIL: build failed, see results/$name.build.log"
            tail -5 "$RESULTS/$name.build.log" | sed 's/^/        /'
            FAIL=$((FAIL+1)); FAILED_CASES="$FAILED_CASES $name"; continue
        fi
    fi

    BIN="$CODE/build/complab"
    [ -x "$BIN" ] || BIN="$CODE/complab"
    if [ ! -x "$BIN" ]; then
        echo "  SKIP: no complab executable found."
        SKIP=$((SKIP+1)); continue
    fi

    # Each case runs in its own scratch directory. CompLaB.xml must be in the
    # working directory: the code reads it by that exact name.
    RUN="$RESULTS/$name"
    rm -rf "$RUN"; mkdir -p "$RUN"
    cp "$dir/CompLaB.xml" "$RUN/"
    cp -r "$dir/input"    "$RUN/"
    mkdir -p "$RUN/output" "$RUN/src"
    # COBRApy imports its module by name at run time from <src_path>
    [ -f "$CODE/src/complab3d_cobrapy.py" ] && cp "$CODE/src/complab3d_cobrapy.py" "$RUN/src/"

    # Under Slurm, launch through srun so the job's ranks are used. Outside
    # it, run serially: these domains are far too small for MPI to help, and a
    # login-node mpirun is usually unwelcome anyway.
    LAUNCH=""
    if [ -n "${SLURM_JOB_ID:-}" ]; then LAUNCH="${MPIRUN_LAUNCHER:-srun}"; fi

    echo "    running ${LAUNCH:+with $LAUNCH }..."
    START=$(date +%s)
    ( cd "$RUN" && ${LAUNCH} "$BIN" ) > "$RUN/run.log" 2>&1
    RC=$?
    ELAPSED=$(( $(date +%s) - START ))

    # ---- pass or fail -----------------------------------------------------
    # A zero exit is necessary but not sufficient: this code can finish
    # cleanly having produced nonsense. So also insist that nothing reported a
    # negative concentration and that no line says Terminating.
    PROBLEM=""
    [ $RC -ne 0 ] && PROBLEM="exit code $RC"
    if grep -q "Terminating" "$RUN/run.log" 2>/dev/null; then
        PROBLEM="${PROBLEM:+$PROBLEM; }input rejected at start-up"
    fi
    NEG=$(grep -c "\[NEG!\]" "$RUN/run.log" 2>/dev/null || echo 0)
    if [ "$NEG" -gt 0 ]; then
        PROBLEM="${PROBLEM:+$PROBLEM; }$NEG negative-concentration warning(s)"
    fi

    if [ -z "$PROBLEM" ]; then
        echo "  PASS  (${ELAPSED}s)"
        PASS=$((PASS+1))
    else
        echo "  FAIL  (${ELAPSED}s): $PROBLEM"
        echo "        log: results/$name/run.log"
        tail -12 "$RUN/run.log" 2>/dev/null | sed 's/^/        /'
        FAIL=$((FAIL+1)); FAILED_CASES="$FAILED_CASES $name"
    fi
done

echo ""
echo "==================================================================="
echo "  $PASS passed, $FAIL failed, $SKIP skipped"
[ -n "$FAILED_CASES" ] && echo "  failed:$FAILED_CASES"
echo "  logs and output in $RESULTS"
echo "==================================================================="

# ---------------------------------------------------------------------------
# What this script does NOT check
#
# It checks that each case runs and stays physical. It does not check the
# NUMBERS. The per-case README says what to expect, and several of them have a
# right answer worth confirming by hand:
#
#   02  the steady profile is a straight line from 1 to 0
#   03  total A + total C constant, and total B + total C constant
#   09  growth settles near 8 mmol/gDW/h, from the stoichiometry
#   10  must agree with 09
#   13  porosity falls in steps;  14  porosity rises in steps
#   15  porosity falls then rises, and settles rather than oscillating
#
# Automating those means parsing the VTK output, which is worth doing once the
# suite has run against real Palabos and the log format is settled.
# ---------------------------------------------------------------------------

exit $([ $FAIL -eq 0 ] && echo 0 || echo 1)
