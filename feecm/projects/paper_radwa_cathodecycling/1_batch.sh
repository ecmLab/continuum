#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH -t 08:00:00
#SBATCH --ntasks-per-node=116
#SBATCH --array=1-3
#SBATCH --job-name=MOOSE_LPSC_Sweep
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vazquezm
#SBATCH -A mat250014p
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err

# ---------------------------------------------------------------
# Parameter Sweep: 3 Young's Moduli (ymod_lpsc) values = 3 tasks
# Array task IDs 1-3 are mapped to the ymod_lpsc_vals array.
#
# For sync_times.i run the following:
# python -c "print(\"output_times_str = '\" + ' '.join([str(round(i*0.1, 1)) for i in range(10001)]) + \"'\")" > sync_times.i
#
# How to run:
#   Make sure to run "mkdir -p logs rst runs"
#   sbatch 1_batch.sh
# ---------------------------------------------------------------

set -x

# ---- Paths (SET THESE for your HPC staging location) ----------
BASE=/ocean/projects/mat250014p/shared/projects/paper_Radwa_MOOSECathodeCycling/contact_loss
EXE="$BASE/contact_loss-opt"
cd "$BASE"

# ---- Modules -------------------------------------------------
module purge
module load gcc/10.2.0
module load openmpi/4.0.5-gcc10.2.0
module load anaconda3/2022.10

# ---- MPI compiler wrappers -----------------------------------
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F90=mpif90
export F77=mpif77

# ---- Parameter arrays ---------------------------------------
# Array of the 3 requested ymod_lpsc values
ymod_lpsc_vals=(370 550 22000)

# Map SLURM array ID (1-3) to bash array index (0-2)
TASK_ID=$(( SLURM_ARRAY_TASK_ID - 1 ))
YMOD=${ymod_lpsc_vals[$TASK_ID]}

TAG="E${YMOD}_LPSC"

echo "=============================================="
echo "Task $SLURM_ARRAY_TASK_ID -> ymod_lpsc=$YMOD MPa (tag=$TAG)"
echo "=============================================="

# Ensure concurrent runs do not write to the same files
RUNDIR="$BASE/runs/$TAG"
rm -rf "$RUNDIR" 
mkdir -p "$RUNDIR" logs
cd "$RUNDIR"

# Use SLURM's actual core count
NP=${SLURM_NTASKS:-$SLURM_NTASKS_PER_NODE}

# Run MOOSE, overriding the ymod_lpsc variable from 1_contact.i
mpirun -np "$NP" "$EXE" \
    -i "$BASE/1_contact.i" \
    "ymod_lpsc=$YMOD" \
    "Outputs/file_base=${TAG}_out"

EXIT_CODE=$?
echo "MOOSE exit code: $EXIT_CODE"
exit $EXIT_CODE