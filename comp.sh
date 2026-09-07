#!/bin/bash
# ============================================================================
#  CompLB3D  --  Slurm submit script for Sapelo2 (GACRC, University of Georgia)
# ============================================================================
#
#  BUILD FIRST, then submit.  See the block at the bottom of this file for the
#  exact build commands.
#
#  Submit from the directory that holds ./complab and CompLaB.xml:
#
#      sbatch comp.sh
#
#  CompLaB.xml is read from the WORKING DIRECTORY, by that exact name.  The
#  directory named in <output_path> must already exist -- the code does not
#  create it, and a run that cannot write its output dies late rather than
#  early.  This script makes it for you.
#
#  Part of the CompLaB program.  GNU Affero General Public License v3 or later.
#  Meile Lab, University of Georgia.  shahram.asgari@uga.edu
# ============================================================================

#SBATCH --job-name=complab3d           # shown in squeue
#SBATCH --partition=batch              # Sapelo2's standard queue; check with: sinfo -s
#SBATCH --nodes=1                      # one node
#SBATCH --ntasks=32                    # MPI ranks; check a node's core count with: sinfo -o "%n %c"
#SBATCH --mem=64gb                     # memory for the whole job
#SBATCH --time=04:00:00                # wall clock; the job is killed past this
#SBATCH --output=complab_%j.out        # stdout, %j = job id
#SBATCH --error=complab_%j.err         # stderr kept separate, which makes failures readable
#SBATCH --mail-type=END,FAIL           # mail on finish or failure
#SBATCH --mail-user=shahram.asgari@uga.edu

# ---------------------------------------------------------------------------
# Working directory
# ---------------------------------------------------------------------------
cd "$SLURM_SUBMIT_DIR" || { echo "cannot cd to $SLURM_SUBMIT_DIR"; exit 1; }
mkdir -p output

# ---------------------------------------------------------------------------
# Modules
#
#   THESE VERSION STRINGS ARE A STARTING POINT, NOT A GUARANTEE.  Confirm them
#   once with
#       module spider foss
#       module spider CMake
#       module spider GLPK
#   and edit the three lines below to whatever Sapelo2 actually offers.
#
#   The compiler toolchain and GLPK must come from the SAME generation --
#   GCCcore-11.3.0 here.  Mixing them links, then fails at run time with
#   missing symbols, which is a miserable thing to debug.
#
#   Load exactly what you BUILT with.  An executable built against one MPI and
#   launched under another will start, print nothing useful, and hang.
# ---------------------------------------------------------------------------
module purge
module load foss/2022a                      # GCC + OpenMPI + FlexiBLAS
module load GLPK/5.0-GCCcore-11.3.0         # only needed if built with -DENABLE_GLPK=ON
# module load SciPy-bundle/2022.05-foss-2022a   # only for -DENABLE_COBRAPY=ON

# ---------------------------------------------------------------------------
# A short record at the top of the log, so a result can be traced to a run
# ---------------------------------------------------------------------------
echo "============================================================"
echo " job id     : $SLURM_JOB_ID"
echo " host       : $(hostname)"
echo " nodes      : $SLURM_JOB_NUM_NODES"
echo " tasks      : $SLURM_NTASKS"
echo " directory  : $(pwd)"
echo " started    : $(date)"
echo "============================================================"
module list 2>&1
echo "============================================================"

if [ ! -x ./complab ]; then
    echo "ERROR: ./complab is not here or is not executable. Build it first."
    exit 1
fi
if [ ! -f ./CompLaB.xml ]; then
    echo "ERROR: CompLaB.xml is not in this directory. The code reads it by that exact name."
    exit 1
fi

# ---------------------------------------------------------------------------
# Run
#
#   srun, not mpirun: it inherits the allocation from Slurm, so the rank count
#   always matches --ntasks and there is nothing to keep in step by hand.
#
#   For a SERIAL run -- which is the right choice for the example cases, they
#   are 24x24x6 and MPI only adds overhead -- set --ntasks=1 above and replace
#   the line below with:   ./complab
# ---------------------------------------------------------------------------
srun ./complab
status=$?

echo "============================================================"
echo " finished   : $(date)"
echo " exit code  : $status"
echo "============================================================"
exit $status

# ============================================================================
#  BUILDING, for reference
# ============================================================================
#
#   cd /scratch/sa01687/complab3Dsurrogate
#
#   module purge
#   module load foss/2022a
#   module load CMake/3.24.3-GCCcore-11.3.0
#   module load GLPK/5.0-GCCcore-11.3.0
#
#   mkdir -p build && cd build
#   cmake -DENABLE_GLPK=ON ..
#   make -j8
#   cd ..
#
#   CMakeLists.txt finds Palabos at versionControl/palabos-v2.3.0 next to
#   itself.  If your copy is named or placed differently:
#
#       cmake -DENABLE_GLPK=ON -DPALABOS_ROOT=/full/path/to/palabos-v2.3.0 ..
#
#   -DENABLE_GLPK=ON is REQUIRED if any run will use <reaction_type>glpk, or
#   <surrogate><train_if_missing>true -- training sweeps a linear program even
#   when no organism runs FBA at every voxel.
#
#   The executable lands at ../complab, next to CMakeLists.txt, because
#   CMAKE_RUNTIME_OUTPUT_DIRECTORY is set to "../" at the top of that file.
#
#   THE BUILD ALSO NEEDS THESE AT THE TOP LEVEL, beside CMakeLists.txt:
#       defineKinetics.hh
#       defineAbioticKinetics.hh
#       surrogateModel.hh
#   src/complab.cpp includes them as "../defineKinetics.hh" and friends.  The
#   build stops on the first one that is missing.
# ============================================================================
