#!/bin/bash
# Job name:
#SBATCH --job-name=PrimaryParticle
#
# Partition:
#SBATCH --partition=savio
#
# Processors:
#SBATCH --nodes=4
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=24:00:00
#

lmpl=/global/home/users/howietu/software/LIGGGHTS-PUBLIC/src/lmp_auto

## load or unload modules:
#module unload gcc
#module unload cuda

## Commands to run:
mpirun -np 80 $lmpl < cohesion.in
