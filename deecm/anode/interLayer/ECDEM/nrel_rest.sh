#!/bin/bash
#SBATCH --job-name="REST_Ag_M1C5_vr15_fp250_sp14"
#SBATCH --account=interlayer
#SBATCH --time=4:00:00
#SBATCH --ntasks-per-node=94
#SBATCH --nodes=1
#SBATCH --mail-user=jmv8431@rit.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Set the path to your LIGGGHTS executable
export PATH=/kfs3/scratch/jmv8431/LIGGGHTS-PUBLIC/src:$PATH

# --- SIMULATION WORKFLOW ---
# The simulations will run in the order listed below.

# This first step must complete and write a restart file before others can run.
srun lmp_auto -in rest.in