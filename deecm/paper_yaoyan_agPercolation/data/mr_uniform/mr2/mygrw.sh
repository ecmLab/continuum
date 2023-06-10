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
#SBATCH --cpus-per-task=4
#
# Wall clock limit:
#SBATCH --time=40:00:00
#SBATCH --mem=6g

## Commands to run:
srun lmp -log none -in grw.in
## srun --nodes=1 --ntasks=2 --cpus-per-task=2 lmp -in lmp_mr.in 
