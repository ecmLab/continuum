#!/bin/bash
# Job name:
#SBATCH --job-name=Ag_aC
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
#SBATCH --time=01:50:00
#SBATCH --mem=30g

## Commands to run:
srun matlab -nodisplay -nosplash < mdl_all.m
