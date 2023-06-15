#!/bin/bash
# Job name:
#SBATCH --job-name=Ag_aC
#
# Partition:
#SBATCH --partition=tier3
#
# Processors:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=10:00:00
#SBATCH --mem=6g

## Commands to run:
srun matlab -nodisplay -nosplash < mdl_run.m
