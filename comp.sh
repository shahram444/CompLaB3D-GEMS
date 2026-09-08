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
#SBATCH --ntasks=1                     # SERIAL. See "HOW MANY RANKS" below before raising it.
#SBATCH --mem=16gb                     # memory for the whole job
#SBATCH --time=02:00:00                # wall clock; the job is killed past this
#SBATCH --output=complab_%j.out        # stdout, %j = job id
#SBATCH --error=complab_%j.err         # stderr kept separate, which makes failures readable
#SBATCH --mail-type=END,FAIL           # mail on finish or failure
#SBATCH --mail-user=shahram.asgari@uga.edu

# ---------------------------------------------------------------------------
# HOW MANY RANKS
#
#   ONE, for every example case in this repository. They are 24x26x8 = 4992
#   voxels. Split across 32 ranks that is 150 voxels each, and the ranks spend
#   longer talking to each other than computing. A serial run finishes them in
#   seconds.
#
#   Serial also avoids the interconnect entirely, which matters here: UCX 1.12.1
#   in the foss/2022a stack aborts during connection setup on some Sapelo2
#   nodes, with
#
#       ib_iface.c:742  Assertion `gid->global.interface_id != 0' failed
#
#   That is a bug in the cluster's InfiniBand layer, not in this code, and it
#   kills the job before the solver has read anything. If you hit it on a run
#   that genuinely needs several ranks, uncomment ONE of the two blocks in the
#   RUN section below to route around it, and report the backtrace to GACRC.
#
#   Raise --ntasks only for a domain large enough to pay for the communication,
#   and read the MPI note in pipelines/C_run/README.md first: dissolution loses
#   mineral mass across block boundaries, so cases 14 and 15 must stay serial
#   whatever their size.
# ---------------------------------------------------------------------------

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
#   and edit the lines below to whatever Sapelo2 actually offers.
#
#   The compiler toolchain and GLPK must come from the SAME generation --
#   GCCcore-11.3.0 here.  Mixing them links, then fails at run time with
#   missing symbols, which is a miserable thing to debug.
#
#   Load exactly what you BUILT with.  An executable built against one MPI and
#   launched under another will start, print nothing useful, and hang.
# ---------------------------------------------------------------------------
module purge
module load foss/2022a                        # GCC + OpenMPI + FlexiBLAS
# module load GLPK/5.0-GCCcore-11.3.0         # only if built with -DENABLE_GLPK=ON
# module load SciPy-bundle/2022.05-foss-2022a # only if built with -DENABLE_COBRAPY=ON

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
#   One rank runs the executable directly. No srun, no mpirun: neither adds
#   anything to a single process, and going through the launcher is what drags
#   the interconnect into a job that has no use for it.
#
#   FOR A PARALLEL RUN, set --ntasks above and use srun instead:
#
#       srun ./complab
#
#   srun rather than mpirun, because it inherits the allocation from Slurm and
#   the rank count always matches --ntasks with nothing to keep in step by hand.
#
#   If that aborts inside UCX (see HOW MANY RANKS above), uncomment ONE of these
#   before the srun line. The first keeps OpenMPI's UCX path but drops the
#   InfiniBand transport; the second bypasses UCX altogether. Both are slower
#   than a working interconnect and both complete.
#
#       export UCX_TLS=tcp,self,sm
#
#       export OMPI_MCA_pml=ob1
#       export OMPI_MCA_btl=self,vader,tcp
# ---------------------------------------------------------------------------
if [ "${SLURM_NTASKS:-1}" -gt 1 ]; then
    srun ./complab
else
    ./complab
fi
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
#   cd /scratch/<you>/CompLaB3D-GEMS
#
#   module purge
#   module load foss/2022a
#   module load CMake/3.23.1-GCCcore-11.3.0
#   module load GLPK/5.0-GCCcore-11.3.0        # only for the flux-balance cases
#
#   mkdir -p build && cd build
#   cmake ..                                   # add -DENABLE_GLPK=ON if needed
#   make -j8
#   cd ..
#
#   Build in an interactive job, not on a login node: this compiles all of
#   Palabos into a static library, which takes 10 to 20 minutes and several GB.
#
#       interact -c 8 --mem 40gb --time 2:00:00
#
#   CMakeLists.txt finds Palabos at versionControl/palabos-v2.3.0 next to
#   itself.  If your copy is named or placed differently:
#
#       cmake -DPALABOS_ROOT=/full/path/to/palabos-v2.3.0 ..
#
#   or set PALABOS_ROOT in your environment, which is what you want if you build
#   inside a case assembled by scripts/setup_case.sh -- that directory has no
#   versionControl/ beside it, so the built-in fallback does not resolve.
#
#   -DENABLE_GLPK=ON is REQUIRED if any run will use <reaction_type>glpk, or
#   <surrogate><train_if_missing>true -- training sweeps a linear program even
#   when no organism runs FBA at every voxel.
#
#   Do NOT build with -DENABLE_MPI=OFF. The option exists, but the solver calls
#   MPI reductions unguarded in about twenty places, so a serial-only build does
#   not compile. Leave it ON and control the rank count with --ntasks, which is
#   what this script does.
#
#   The executable lands at ../complab, next to CMakeLists.txt, because
#   CMAKE_RUNTIME_OUTPUT_DIRECTORY is set to "../" at the top of that file.
#
#   THE BUILD ALSO NEEDS THESE AT THE TOP LEVEL, beside CMakeLists.txt:
#       defineKinetics.hh
#       defineAbioticKinetics.hh
#       surrogateModel.hh
#   src/complab.cpp includes them as "../defineKinetics.hh" and friends.  The
#   build stops on the first one that is missing.  scripts/setup_case.sh lays
#   all three down for you; building at the top of the repository instead means
#   copying them yourself:
#
#       cp config/kinetics/defineKinetics.default.hh        defineKinetics.hh
#       cp config/kinetics/defineAbioticKinetics.default.hh defineAbioticKinetics.hh
#       cp src/surrogateModel.hh                            surrogateModel.hh
#
#   A case that ships its own kinetics/ directory (03, 05-08, 12, 13-15, 19, 20)
#   overrides the matching default, and changing it means rebuilding.  Palabos is
#   already compiled by then, so that rebuild is about a minute.
# ============================================================================
