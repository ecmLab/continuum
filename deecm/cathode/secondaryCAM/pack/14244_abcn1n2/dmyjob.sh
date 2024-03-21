#!/bin/bash
# Job name:
#SBATCH --job-name=PrimaryParticle
#
# Partition:
#SBATCH --partition=savio
#
#QoS:
#SBATCH --qos=savio_debug
#
# Processors:
#SBATCH --nodes=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=00:30:00
#

lmpl=/global/home/users/howietu/software/LIGGGHTS-PUBLIC/src/lmp_auto

## load or unload modules:
#module unload gcc
#module unload cuda

## Commands to run:
mpirun -np 24 $lmpl < in.lmp
