#!/bin/bash
# Job name:
#SBATCH --job-name=NMC_LPS
#
# Partition:
#SBATCH --partition=tier3
#
# Processors:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#
# Wall clock limit:
#SBATCH --time=120:30:00
#SBATCH --mem=20g


## Commands to run:
spack load lammps@20211027
srun lmp -log none -in mr.in
