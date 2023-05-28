#!/bin/bash
# Job name:
#SBATCH --job-name=NMC_LPS
#
# Partition:
#SBATCH --partition=tier3
#
# Processors:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=00:30:00
#SBATCH --mem=10g


## Commands to run:

srun lmp -log none -in mr.in
