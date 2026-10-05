#!/usr/bin/env bash
#SBATCH --job-name=abc_sweep
#SBATCH --partition=batch
#SBATCH --array=0-59%20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/slurm_%A_%a.out

# One process per case, not one big job. Array tasks are separate processes, so
# one crash cannot take the rest down, and run_one.sh never propagates a failure.
#
#   sbatch submit.sbatch                 the whole sweep
#   sbatch --array=3,7,12 submit.sbatch  three cases again
#
# make_campaign.py retry prints the --array list for whatever failed.

set -u
cd "$(dirname "$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}' | awk '{print $1}')")" 2>/dev/null || cd "$SLURM_SUBMIT_DIR"
mkdir -p logs
CASE=$(awk -F, -v n="$SLURM_ARRAY_TASK_ID" 'NR>1 && $1==n {print $2}' runs.csv)
if [ -z "$CASE" ]; then echo "no case for array index $SLURM_ARRAY_TASK_ID"; exit 0; fi
./run_one.sh "runs/$CASE" "./complab"
