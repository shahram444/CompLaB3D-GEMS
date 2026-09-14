#!/bin/bash
#
# 07_biotic_lattice_boltzmann
# The same organism, biomass spread by lattice Boltzmann.
#
# THIS SCRIPT ONLY RUNS THE SOLVER. It does not compile anything.
# Build ./complab first in an interactive session: see COMPILE.txt in this
# folder, which lists the commands line by line.
#
# Then, from this folder:     sbatch run.sh
# Watch it:                   squeue -u $USER
# Read the log while it runs:  tail -f output/run.log
#
# ---------------------------------------------------------------------------
# THE SLURM HEADER. Every line below is read by sbatch, not by bash, so each
# comment sits on its own line above the directive it explains.
# ---------------------------------------------------------------------------

# what this job is called in squeue, and the stem of the two log files below
#SBATCH --job-name=complab_07_biotic_lattice_boltzmann
# the lab partition
#SBATCH --partition=meile_p
# mail on start, on finish and on failure
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sa01687@uga.edu
# one node is enough: this is a 24 x 26 x 8 grid
#SBATCH --nodes=1
# ONE rank, deliberately. Read the note at the bottom of this file before raising it.
#SBATCH --ntasks=1
# memory for the whole job
#SBATCH --mem=16gb
# wall clock ceiling; 1000 steps on this grid finishes well inside it
#SBATCH --time=01:00:00
# %x is the job name set above, %j is the job id SLURM assigns
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# ---------------------------------------------------------------------------
# THE JOB
# ---------------------------------------------------------------------------

cd $SLURM_SUBMIT_DIR                       # SLURM starts you in $HOME; the case lives where you submitted from

module purge                               # drop whatever your login shell had loaded
module load foss/2022a                     # GCC 11.3 and OpenMPI: the toolchain complab was built against
module load SciPy-bundle/2022.05-foss-2022a  # provides python3 for preprocess.py and postprocess.py

python3 preprocess.py                      # rebuild input/geometry.dat, and report porosity

mkdir -p output                            # the solver writes here and will not create it

# One rank, so no srun. READ THE START-UP LINES in the log: the geometry, every
# enabled feature and each reaction path are echoed before the first step.
./complab CompLaB.xml 2>&1 | tee output/run.log

python3 postprocess.py 2>&1 | tee output/checks.log   # the verdict: is this run worth believing

# ---------------------------------------------------------------------------
# WHAT TO OPEN IN PARAVIEW
# ---------------------------------------------------------------------------
#
#   donor_*.vti
#   product_*.vti
#   rate_<species>_*.vti      the reaction rate as a field, mol/L/s
#   Bug_*.vti
#
# ---------------------------------------------------------------------------
# WHY --ntasks=1 AND NO srun
# ---------------------------------------------------------------------------
#
# A multi-rank run on this cluster aborts inside UCX before the first step:
#
#   ib_iface.c:742  Assertion `gid->global.interface_id != 0` failed
#
# It is an InfiniBand device-address bug in UCX 1.12.1, not a CompLaB fault.
# These example cases are small enough that one rank is quick anyway. If you do
# need more ranks on a production geometry, raise --ntasks, put srun back in
# front of ./complab, and expect to work around UCX first.
