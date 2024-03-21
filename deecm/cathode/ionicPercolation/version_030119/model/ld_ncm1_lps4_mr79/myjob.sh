#!/bin/bash
# Job name:
#SBATCH --job-name=NMC_LPS
#
# Partition:
#SBATCH --partition=savio2
#
# Processors:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=24:00:00
#

## Commands to run:
module load matlab
mpirun -np 24 matlab -nodisplay -nosplash < mdl_all.m
