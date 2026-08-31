#!/bin/bash                            # Run this script with the bash shell.

#SBATCH --account=PROJECT             # Charge the compute time to this project account (the name from 'gbalance'). Required, or the job is rejected.
#SBATCH --job-name=complab3d          # A label for the job, shown in 'squeue' so it's easy to spot.
#SBATCH --partition=normal            # Which queue to run in. 'normal' is the standard compute queue (up to 100 nodes, 2-day max).
                                      
#SBATCH --nodes=1                     # Ask for 1 compute node.
#SBATCH --ntasks-per-node=36          # Run 36 MPI processes on that node (a Tahoma node has 36 cores).
                                                                    
#SBATCH --time=04:00:00               # Max wall-clock time (4 hours). The job is killed if it runs longer, so set it a bit above the expected runtime.
                                                                  
#SBATCH --output=complab_%j.out       # Send normal program output to this file (%j = job ID,  so each run gets its own file).
                                      
#SBATCH --error=complab_%j.err        # Send error messages to a separate file (helps debugging).
#SBATCH --mail-type=END,FAIL          # Email us when the job finishes or fails.
#SBATCH --mail-user=your_email@pnnl.gov   # Where those emails go (put your real address here).

cd $SLURM_SUBMIT_DIR                  # Move into the folder we submitted from, the one holding
                                      # 'complab' and 'CompLaB.xml', since the program reads
                                      # CompLaB.xml from its current directory.
module purge                          # Clear default modules for a clean, conflict-free environment.
module load gcc openmpi               # Load the SAME compiler + MPI we built with, so the
                                      # executable starts correctly (must match Step 6).

srun ./complab                        # Launch 'complab' in parallel across the 36 tasks Slurm
                                      # allocated. 'srun' is what actually starts the MPI run.
