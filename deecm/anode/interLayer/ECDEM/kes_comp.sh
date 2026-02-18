#!/bin/bash
#SBATCH --job-name=COM_Ag_M1C5_vr15_fp250_sp14
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=03:00:00
#SBATCH --account=mat250014p
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vazquezm

# Load required modules LIGGGHTS
module load gcc/10.2.0 openmpi/5.0.3-gcc13.2.1 cuda/11.7.1

# The simulations will run in the order listed below.
mpiexec -np $SLURM_NTASKS /ocean/projects/mat250014p/shared/software/LIGGGHTS-PUBLIC/src/lmp_auto -in pck.in