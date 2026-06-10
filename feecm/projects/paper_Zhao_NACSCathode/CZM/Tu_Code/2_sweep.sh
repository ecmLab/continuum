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
#   sbatch --array=91-991:100 2_sweep.sh (czm_B=1, Hv=30, ymod=100-1000)
#
# Note: 1000 tasks throttled to 10 concurrent (%10). At ~1 hr/task
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
ymod_nacs=(100 200 300 400 500 600 700 800 900 1000)  # MPa
Hv_nacs=(3 6 9 12 15 18 21 24 27 30)                  # MPa
czm_B=(1 2 3 4 5 6 7 8 9 10)                         # 1=Baseline ... 10

TASK_ID=$(( SLURM_ARRAY_TASK_ID - 1 ))
i_ymod=$(( TASK_ID / 100 ))
i_hv=$(( (TASK_ID / 10) % 10 ))
i_czm=$(( TASK_ID % 10 ))
YMOD=${ymod_nacs[$i_ymod]}
HV=${Hv_nacs[$i_hv]}
CZM=${czm_B[$i_czm]}

TAG="E${YMOD}_H${HV}_C${CZM}"

echo "=============================================="
echo "Task $SLURM_ARRAY_TASK_ID -> ymod_nacs=$YMOD MPa, Hv_nacs=$HV MPa, czm_B=$CZM  (tag=$TAG)"
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
mpirun -np "$NP" "$EXE" \
    -i "$BASE/2_czm.i" \
    "ymod_nacs=$YMOD" \
    "Hv_nacs=$HV" \
    "czm_B=$CZM" \
    "Outputs/file_base=${TAG}_out"

EXIT_CODE=$?
echo "MOOSE exit code: $EXIT_CODE"
exit $EXIT_CODE