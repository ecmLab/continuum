#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 10:00:00
#SBATCH --ntasks-per-node=64
#SBATCH --array=1-1000%10
#SBATCH --job-name=MOOSE_CZMSweep
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vazquezm
#SBATCH -A mat250014p
#SBATCH --output=logs_czm/job_%A_%a.out
#SBATCH --error=logs_czm/job_%A_%a.err

# ---------------------------------------------------------------
# Parameter Sweep: 10 Young's Moduli x 10 Hardness x 10 czm_B
#                  = 1000 tasks
#
# Array task IDs 1-1000 are mapped to (i_ymod, i_hv, i_czm) index
# triplets using row-major (lexicographic) order:
#
#   TASK_ID (0-indexed) = i_ymod * 100 + i_hv * 10 + i_czm
#
#   i_ymod = (TASK_ID) / 100        → Young's Modulus  (outermost)
#   i_hv   = (TASK_ID / 10) % 10    → Hardness         (middle)
#   i_czm  = (TASK_ID) % 10         → czm_B            (innermost)
#
# Example:
#   SLURM_ARRAY_TASK_ID=1    → TASK_ID=0   → ymod=100,  Hv=3,  C=1
#   SLURM_ARRAY_TASK_ID=2    → TASK_ID=1   → ymod=100,  Hv=3,  C=2
#   SLURM_ARRAY_TASK_ID=10   → TASK_ID=9   → ymod=100,  Hv=3,  C=10
#   SLURM_ARRAY_TASK_ID=11   → TASK_ID=10  → ymod=100,  Hv=6,  C=1
#   SLURM_ARRAY_TASK_ID=101  → TASK_ID=100 → ymod=200,  Hv=3,  C=1
#   SLURM_ARRAY_TASK_ID=1000 → TASK_ID=999 → ymod=1000, Hv=30, C=10
#
# How to run:
#   Make sure to run "mkdir -p logs_czm rst_czm runs_czm"
#   sbatch 2_sweep.sh
#   sbatch --array=1,1000 2_sweep.sh   (just the first and last tasks)
#   sbatch --array=1-1000:2 2_sweep.sh (every other task)
#   sbatch --array=1-10 2_sweep.sh     (czm_B sweep at ymod=100, Hv=3)
#   sbatch --array=991-1000 2_sweep.sh (czm_B sweep at ymod=1000, Hv=30)
#
# Note: 1000 tasks throttled to 10 concurrent (%10). At ~1 hr/task
# this serializes to a long wall clock; raise the %N throttle (e.g.
# %20, %50) and/or bump -t if CZM convergence is slower than the
# non-CZM runs.
# ---------------------------------------------------------------

set -x

# ---- Working directory ----------------------------------------
BASE=/ocean/projects/mat250014p/shared/projects/paper_Zhao_NACSCathode/contact_loss
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
ymod_se=(100 147 195 242 289 337 384 432 479 526 574 621 668 716 763 811 858 905 953 1000)  # MPa
Hv_se=(20 22 23 25 26 28 29 31 33 34 36 37 39 41 42 44 45 47 48 50)                  # MPa
czm_B=(1 2 3 4 5 6 7 8 9 10)                         # Need to Find a Value that Works for All 20x20 Ymod-Hv Pairs (Could Just Test the Extremes of Hv and Ymod)

TASK_ID=$(( SLURM_ARRAY_TASK_ID - 1 ))
i_ymod=$(( TASK_ID / 100 ))
i_hv=$(( (TASK_ID / 10) % 10 ))
i_czm=$(( TASK_ID % 10 ))
YMOD=${ymod_se[$i_ymod]}
HV=${Hv_se[$i_hv]}
CZM=${czm_B[$i_czm]}

TAG="E${YMOD}_H${HV}_C${CZM}"

echo "=============================================="
echo "Task $SLURM_ARRAY_TASK_ID -> ymod_se=$YMOD MPa, Hv_se=$HV MPa, czm_B=$CZM  (tag=$TAG)"
echo "=============================================="

# so no two concurrent runs write to the same Exodus / restart files.
RUNDIR="$BASE/runs_czm/$TAG"
mkdir -p "$RUNDIR" logs_czm
cd "$RUNDIR"

# Use SLURM's actual core count rather than assuming the var is set.
NP=${SLURM_NTASKS:-$SLURM_NTASKS_PER_NODE}

# Run from RUNDIR; point MOOSE at the input by absolute path.
# Outputs/file_base is overridden too so the output name is unique
# even if your Outputs block hardcodes a base name.
mpirun -np "$NP" "$BASE/contact_loss-opt" \
    -i "$BASE/2_czm.i" \
    "ymod_se=$YMOD" \
    "Hv_se=$HV" \
    "czm_B=$CZM" \
    "Outputs/file_base=${TAG}_out"

EXIT_CODE=$?
echo "MOOSE exit code: $EXIT_CODE"
exit $EXIT_CODE