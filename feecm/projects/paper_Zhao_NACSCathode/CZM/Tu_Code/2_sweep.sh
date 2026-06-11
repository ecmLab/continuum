#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 2:00:00
#SBATCH --ntasks-per-node=48
#SBATCH --array=1-400%20
#SBATCH --job-name=MOOSE_CZMSweep
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vazquezm
#SBATCH -A mat250014p
#SBATCH --output=logs_czm/job_%A_%a.out
#SBATCH --error=logs_czm/job_%A_%a.err

# ---------------------------------------------------------------
# Parameter Sweep: 20 Young's Moduli x 20 Hardness x 1 czm_B
#                  = 400 tasks
#
# Array task IDs 1-400 are mapped to (i_ymod, i_hv) index
# pairs using row-major (lexicographic) order. czm_B is fixed
# at a single value, so i_czm is always 0:
#
#   TASK_ID (0-indexed) = i_ymod * 20 + i_hv
#
#   i_ymod = (TASK_ID) / 20         → Young's Modulus  (outermost)
#   i_hv   = (TASK_ID) % 20         → Hardness         (innermost)
#   i_czm  = 0                      → czm_B            (fixed)
#
# Example:
#   SLURM_ARRAY_TASK_ID=1    → TASK_ID=0   → ymod=100,  Hv=20, C=1
#   SLURM_ARRAY_TASK_ID=2    → TASK_ID=1   → ymod=100,  Hv=22, C=1
#   SLURM_ARRAY_TASK_ID=20   → TASK_ID=19  → ymod=100,  Hv=50, C=1
#   SLURM_ARRAY_TASK_ID=21   → TASK_ID=20  → ymod=147,  Hv=20, C=1
#   SLURM_ARRAY_TASK_ID=400  → TASK_ID=399 → ymod=1000, Hv=50, C=1
#
# How to run:
#   Make sure to run "mkdir -p logs_czm rst_czm runs_czm"
#   sbatch 2_sweep.sh
#   sbatch --array=1,400 2_sweep.sh      (just the first and last tasks)
#   sbatch --array=1-400:2 2_sweep.sh    (every other task)
#   sbatch --array=1-20 2_sweep.sh       (Hv sweep at ymod=100)
#   sbatch --array=381-400 2_sweep.sh    (Hv sweep at ymod=1000)
#   sbatch --array=1-381:20 2_sweep.sh   (ymod sweep at Hv=20)
#   sbatch --array=20-400:20 2_sweep.sh  (ymod sweep at Hv=50)
#
# Note: 400 tasks throttled to 10 concurrent (%10). At ~1 hr/task
# this serializes to a long wall clock; raise the %N throttle (e.g.
# %20, %50) and/or bump -t if CZM convergence is slower than the
# non-CZM runs.
# ---------------------------------------------------------------

set -x

# ---- Paths (SET THESE for your HPC staging location) ----------
# BASE: cluster directory holding this project's inputs, meshes and the
#   built executable. Update to wherever you stage
#   projects/ecmTest/2026_paper_yangzhao_catholyteMechanics/calculation.
BASE=/ocean/projects/mat250014p/shared/projects/paper_Zhao_NACSCathode/contact_loss
# EXE: the ecm_test optimized executable (app name 'ecm' -> ecm-opt).
#   Copy ecm-opt into $BASE, or set the full path to feecm/ecm_test/ecm-opt.
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
ymod_nacs=(100 147 195 242 289 337 384 432 479 526 574 621 668 716 763 811 858 905 953 1000)  # MPa
Hv_nacs=(20 22 23 25 26 28 29 31 33 34 36 37 39 41 42 44 45 47 48 50)                  # MPa
czm_B=(1)                         # Single value applied to all 20x20 Ymod-Hv pairs (validate at the Hv/Ymod extremes)

TASK_ID=$(( SLURM_ARRAY_TASK_ID - 1 ))
i_ymod=$(( TASK_ID / 20 ))
i_hv=$(( TASK_ID % 20 ))
i_czm=0
YMOD=${ymod_nacs[$i_ymod]}
HV=${Hv_nacs[$i_hv]}
CZM=${czm_B[$i_czm]}

TAG="E${YMOD}_H${HV}_C${CZM}"

echo "=============================================="
echo "Task $SLURM_ARRAY_TASK_ID -> ymod_nacs=$YMOD MPa, Hv_nacs=$HV MPa, czm_B=$CZM  (tag=$TAG)"
echo "=============================================="

# so no two concurrent runs write to the same Exodus / restart files.
RUNDIR="$BASE/runs_czm/$TAG"
rm -rf "$RUNDIR"  # <--- ADD THIS LINE TO AUTO-CLEAN BEFORE RUNNING
mkdir -p "$RUNDIR" logs_czm
cd "$RUNDIR"

# Use SLURM's actual core count rather than assuming the var is set.
NP=${SLURM_NTASKS:-$SLURM_NTASKS_PER_NODE}

# Run from RUNDIR; point MOOSE at the input by absolute path.
# Outputs/file_base is overridden too so the output name is unique
# even if your Outputs block hardcodes a base name.
mpirun -np "$NP" "$EXE" \
    -i "$BASE/2_czm.i" \
    "ymod_nacs=$YMOD" \
    "Hv_nacs=$HV" \
    "czm_B=$CZM" \
    "Outputs/file_base=${TAG}_out"

EXIT_CODE=$?
echo "MOOSE exit code: $EXIT_CODE"
exit $EXIT_CODE