#!/bin/bash
# ============================================================================
#  CompLB3D  --  Slurm submit script for the COBRApy cases (example 10, 16)
#  Sapelo2, GACRC, University of Georgia
# ============================================================================
#
#  This differs from comp.sh in exactly two ways, and both matter:
#
#    1. it loads a Python module, because the executable is linked against
#       libpython and will not even start without it;
#    2. it runs on few ranks, because every rank starts its own embedded
#       interpreter and COBRApy is one to two orders of magnitude slower
#       than GLPK on the same linear program.
#
#      sbatch comp_cobrapy.sh
#
#  Part of the CompLaB program.  GNU Affero General Public License v3 or later.
#  Meile Lab, University of Georgia.  shahram.asgari@uga.edu
# ============================================================================

#SBATCH --job-name=complab3d_cobrapy
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=4                     # keep this small; see the note above
#SBATCH --mem=32gb
#SBATCH --time=04:00:00
#SBATCH --output=complab_%j.out
#SBATCH --error=complab_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shahram.asgari@uga.edu

cd "$SLURM_SUBMIT_DIR" || { echo "cannot cd to $SLURM_SUBMIT_DIR"; exit 1; }
mkdir -p output

# ---------------------------------------------------------------------------
# Modules
#
#   THE PYTHON MODULE MUST BE THE ONE YOU BUILT AGAINST.  If cmake reported
#   "found Python 3.10.4", load Python/3.10.4-GCCcore-11.3.0 here.  A different
#   minor version gives "libpython3.X.so.1.0: cannot open shared object file"
#   before main() runs.
# ---------------------------------------------------------------------------
module purge
module load foss/2022a
module load Python/3.10.4-GCCcore-11.3.0
module load SciPy-bundle/2022.05-foss-2022a
module load GLPK/5.0-GCCcore-11.3.0        # harmless, and needed if built with both

# pip install --user put cobra in ~/.local; make sure it is on the path the
# EMBEDDED interpreter searches, which is the same PYTHONPATH the module sets.
export PYTHONNOUSERSITE=              # unset it: we WANT ~/.local on sys.path

echo "============================================================"
echo " job id     : $SLURM_JOB_ID"
echo " host       : $(hostname)"
echo " tasks      : $SLURM_NTASKS"
echo " directory  : $(pwd)"
echo " started    : $(date)"
echo "============================================================"
module list 2>&1
echo "--- cobra visible to the module python? ---"
python3 -c "import sys, cobra; print(' python ', sys.version.split()[0]); print(' cobra  ', cobra.__version__)" \
    || echo " WARNING: the module python cannot import cobra. The embedded one probably cannot either."
echo "============================================================"

if [ ! -x ./complab ]; then
    echo "ERROR: ./complab is not here or is not executable. Build it first."
    exit 1
fi
if [ ! -f ./CompLaB.xml ]; then
    echo "ERROR: CompLaB.xml is not in this directory."
    exit 1
fi
if [ ! -f ./src/complab3d_cobrapy.py ]; then
    echo "ERROR: src/complab3d_cobrapy.py is missing. <src_path> in CompLaB.xml"
    echo "       says 'src', and the module is imported by name at run time."
    exit 1
fi

srun ./complab
status=$?

echo "============================================================"
echo " finished   : $(date)"
echo " exit code  : $status"
echo "============================================================"
exit $status

# ============================================================================
#  BUILDING for COBRApy, for reference
# ============================================================================
#
#   ON THE LOGIN NODE (compute nodes have no outbound network):
#
#     module purge
#     module load foss/2022a
#     module load Python/3.10.4-GCCcore-11.3.0
#     module load SciPy-bundle/2022.05-foss-2022a
#     python3 -m pip install --user cobra
#     python3 -c "import cobra; print(cobra.__version__)"
#
#   THEN BUILD, with the SAME modules loaded:
#
#     module load CMake/3.24.3-GCCcore-11.3.0
#     module load GLPK/5.0-GCCcore-11.3.0
#     cd /scratch/sa01687/complab3Dsurrogate
#     rm -rf build && mkdir build && cd build
#     cmake -DENABLE_GLPK=ON -DENABLE_COBRAPY=ON \
#           -DPython3_EXECUTABLE=$(which python3) ..
#     make -j8
#     cd ..
#
#   cmake must print   "found Python 3.10.4 headers in ..."   pointing at the
#   module, NOT at /usr/include/python3.6m.  If it points at the system Python
#   the run will fail to import cobra, because cobra was installed for 3.10.
#
#   -DENABLE_GLPK=ON is kept on so one binary runs cases 9, 10, 12 and 16.
# ============================================================================
