#!/bin/bash
# Job name:
#SBATCH --job-name=PrimaryParticle
#
# Partition:
#SBATCH --partition=savio2
#
# Processors:
#SBATCH --nodes=3
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=24:00:00
#

lmpl=/global/home/users/howietu/software/LIGGGHTS-PUBLIC/src/lmp_auto
#gmx=/global/software/sl-7.x86_64/modules/apps/ms/gromacs/5.1.4-cuda/bin/gmx_mpi

## load or unload modules:
#module unload gcc
#module unload cuda

## Commands to run:
mpirun -np 72 $lmpl < in.particle_particle
