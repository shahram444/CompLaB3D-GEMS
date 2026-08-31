#!/bin/bash
#
# Run the CompLB3D example suite on Tahoma (or any Slurm cluster).
#
#   sbatch runExamples.slurm                 all 15 cases
#   sbatch runExamples.slurm 09 10 12        only these
#
# FILL IN --account AND --mail-user BELOW, or the job is rejected.
#
# Part of the CompLaB program.  GNU Affero General Public License v3 or later.

#SBATCH --account=PROJECT                 # from 'gbalance'. Required.
#SBATCH --job-name=complab_examples
#SBATCH --partition=normal

# ---------------------------------------------------------------------------
# ONE NODE, FOUR RANKS, ON PURPOSE.
#
# Every example is 24 x 24 x 6 = 3456 voxels. comp.sh asks for 36 tasks because
# it is written for production runs; here that would give each rank about 96
# voxels, almost all of them envelope, and the halo exchange would cost more
# than the work. Four ranks is already generous. Use 1 if you want the timings
# to mean anything.
#
# These cases are a correctness check, not a scaling study. Do not read
# performance numbers off them.
# ---------------------------------------------------------------------------
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=01:00:00

#SBATCH --output=complab_examples_%j.out
#SBATCH --error=complab_examples_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@pnnl.gov

set -u
cd "$SLURM_SUBMIT_DIR"

module purge
module load gcc openmpi cmake

# Cases 09 and 12 need GLPK; case 10 needs Python with cobrapy importable at
# RUN time, not just at build time. If either module is missing, the runner
# reports those cases as skipped and the rest still go through.
module load glpk        2>/dev/null || echo "NOTE: no glpk module; cases 09 and 12 will be skipped."
module load python      2>/dev/null || true
python3 -c "import cobra" 2>/dev/null || echo "NOTE: cobrapy not importable; case 10 will fail at start-up."

# Point this at your CompLB3D tree. CMakeLists.txt hard-codes
#   /scratch/sa01687/Tahoma/versionControl/palabos-v2.3.0
# so if Palabos has moved, edit lines 85 to 93 there before submitting.
CODE="${CODE:-$SLURM_SUBMIT_DIR/..}"

# The runner builds each case in turn. Only complab.cpp recompiles between
# cases; the Palabos static library is built once and cached, so after the
# first case each rebuild is fast. Switching cmake configuration (three times
# across the suite) costs a reconfigure but still reuses that library.
export MPIRUN_LAUNCHER="srun"
./runAllExamples.sh "$CODE" "$@"

echo
echo "Results in $SLURM_SUBMIT_DIR/results/"
echo "Per-case logs: results/<case>/run.log"
