#!/bin/bash
#SBATCH --job-name=Ag_aC
#SBATCH --account=membrane
#SBATCH --partition=tier3
#
# Processors:
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#
# Wall clock limit:
#SBATCH --time=02:00:00
#SBATCH --mem=10g

## Commands to run:
spack load lammps@20220107
srun lmp -log none -in mr.in
