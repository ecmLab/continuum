#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 02:00:00
#SBATCH --ntasks-per-node=48
#SBATCH --array=1-100
#SBATCH --job-name=MOOSE_ParamSweep
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vazquezm
#SBATCH -A mat250014p
#SBATCH --output=logs/job_%A_%a.out   # %A = job ID, %a = array task ID
#SBATCH --error=logs/job_%A_%a.err

# ---------------------------------------------------------------
# Parameter Sweep: 10 Young's Moduli x 10 Hardness = 100 tasks
# Array task IDs 1-100 are mapped to (i_ymod, i_hv) index pairs
# using row-major order:
#
#   TASK_ID (0-indexed) = i_ymod * 10 + i_hv
#
#   i_ymod = (TASK_ID) / 10   → selects Young's Modulus row
#   i_hv   = (TASK_ID) % 10   → selects Hardness column
#
# Example:
#   SLURM_ARRAY_TASK_ID=1  → TASK_ID=0  → ymod=100, Hv=3
#   SLURM_ARRAY_TASK_ID=11 → TASK_ID=10 → ymod=200, Hv=3
#   SLURM_ARRAY_TASK_ID=100→ TASK_ID=99 → ymod=1000,Hv=30
# How to run:
#   sbatch batch_sweep.sh
#   sbatch --array=1,100 1_sweep.sh (For just the first and last tasks)
#   sbatch --array=1-100:2 submit_sweep.sh (For every other task)
#   sbatch --array=1-10 submit_sweep.sh (For range of tasks)
# ---------------------------------------------------------------

set -x

# ---- Working directory ----------------------------------------
cd /ocean/projects/mat250014p/shared/projects/paper_Zhao_NACSCathode/contact_loss

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

# ---- Parameter arrays (0-indexed, 10 values each) ------------
ymod_se=(100 200 300 400 500 600 700 800 900 1000)  # MPa
Hv_se=(3 6 9 12 15 18 21 24 27 30)                  # MPa

# ---- Map SLURM_ARRAY_TASK_ID (1–100) to array indices --------
TASK_ID=$(( SLURM_ARRAY_TASK_ID - 1 ))   # shift to 0-based
i_ymod=$(( TASK_ID / 10 ))
i_hv=$(( TASK_ID % 10 ))

YMOD=${ymod_se[$i_ymod]}
HV=${Hv_se[$i_hv]}

echo "=============================================="
echo "SLURM Array Task ID : $SLURM_ARRAY_TASK_ID"
echo "Parameter indices   : i_ymod=$i_ymod, i_hv=$i_hv"
echo "ymod_se             : $YMOD MPa"
echo "Hv_se               : $HV MPa"
echo "=============================================="

# ---- Ensure output directories exist ------------------------
mkdir -p rst
mkdir -p logs

# ---- Run MOOSE (override top-level variables via CLI) --------
# MOOSE resolves top-level key=value pairs before any block,
# so passing them on the CLI replaces the values in contact.i
# without editing the input file.
mpirun -np $SLURM_NTASKS_PER_NODE ./contact_loss-opt \
    -i 1_contact.i \
    "ymod_se=$YMOD" \
    "Hv_se=$HV"

EXIT_CODE=$?

echo "MOOSE exit code: $EXIT_CODE"
exit $EXIT_CODE