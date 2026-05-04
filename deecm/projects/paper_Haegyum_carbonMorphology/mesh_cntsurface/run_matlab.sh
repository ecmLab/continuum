#!/bin/bash
#SBATCH --job-name="mat_10"
#SBATCH --account=interlayer
#SBATCH --time=1:00:00
#SBATCH --ntasks-per-node=94
#SBATCH --nodes=1
#SBATCH --mail-user=jmv8431@rit.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load the MATLAB module
module load matlab

matlab -batch mdl_run