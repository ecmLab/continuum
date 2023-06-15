#!/bin/bash
# Job name:
#SBATCH --job-name=Ag_aC
#
# Partition:
#SBATCH --partition=debug
#
# Processors:
##SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#
# Wall clock limit:
#SBATCH --time=20:00:00
#SBATCH --mem=160g

## Commands to run:
srun lmp -log none -in pck_run.in
## srun --nodes=1 --ntasks=2 --cpus-per-task=2 lmp -in lmp_mr.in 
