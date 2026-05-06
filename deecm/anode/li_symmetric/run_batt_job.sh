#!/bin/bash
#SBATCH --job-name=batt_run
#SBATCH --partition=RM
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --time=24:00:00
#SBATCH --account=mat250014p
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=vazquezm

# Echo commands before execution for debugging
set -x

# Change to the directory where the job was submitted
cd $SLURM_SUBMIT_DIR

# Load any required modules (adjust as needed for your environment)
module load gcc/10.2.0 openmpi/5.0.3-gcc13.2.1 cuda/11.7.1

# Run LIGGGHTS with MPI
# Using 256 total MPI tasks (2 nodes * 128 cores/node)
mpiexec -np 256 /ocean/projects/mat250014p/shared/software/LIGGGHTS-PUBLIC/src/lmp_auto -in in.batt_sym_cyc_pack