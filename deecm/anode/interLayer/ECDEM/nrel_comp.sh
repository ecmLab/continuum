#!/bin/bash

#SBATCH --job-name="COM_1M_5C_Ag_100MPA"
#SBATCH --account=interlayer
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=jmv8431@rit.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Set the path to your LIGGGHTS executable
export PATH=/kfs3/scratch/jmv8431/LIGGGHTS-PUBLIC/src:$PATH

# --- SIMULATION WORKFLOW ---
# The simulations will run in the order listed below.

# This first step must complete and write a restart file before others can run.
srun -n $SLURM_CPUS_ON_NODE lmp_auto -in in.compaction