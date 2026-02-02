#!/bin/bash
# Job name:
#SBATCH --job-name=bonds
#
# Partition:
#SBATCH --partition=tier3
#
# Processors:
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=12
#
# Wall clock limit:
#SBATCH --time=30:00:00
#SBATCH --mem=20g


## Commands to run:
spack env activate purewater-24013001
spack load liggghts-flexible-fibers@2023-03-28
srun liggghts < in.carbonbond
#srun instead of mpirun
