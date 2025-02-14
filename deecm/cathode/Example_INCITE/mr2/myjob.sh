#!/bin/bash
# Job name:
#SBATCH --job-name=NMC_LPS
#
# Partition:
#SBATCH --partition=tier3
#
# Processors:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#
# Wall clock limit:
#SBATCH --time=60:30:00
#SBATCH --mem=30g


## Commands to run:
spack load lammps@20211027
srun lmp -log none -in mr.in
