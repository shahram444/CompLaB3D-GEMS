#!/usr/bin/env bash
# One case. Never exits non-zero: it writes status.json either way, so a Slurm
# array task that dies cannot take the collector's triage down with it.
set -u
CASE="${1:?usage: run_one.sh <case directory> [complab binary]}"
BIN="${2:-./complab}"
cd "$CASE" || exit 0
mkdir -p output
. ./env.sh
START=$(date +%s)
"$BIN" CompLaB.xml > output/run.log 2>&1
RC=$?
END=$(date +%s)
STATE=ok
REASON=""
if [ $RC -ne 0 ]; then
  STATE=failed; REASON="exit code $RC"
elif ! grep -q '\[KIN\] abiotic A+B->C' output/run.log; then
  STATE=failed; REASON="no [KIN] line: the binary is not reading PRT_KABIO"
fi
cat > status.json <<EOF
{"state": "$STATE", "reason": "$REASON", "exit_code": $RC, "wall_s": $((END-START))}
EOF
exit 0
