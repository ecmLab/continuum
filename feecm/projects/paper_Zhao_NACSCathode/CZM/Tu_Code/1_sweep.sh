#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 01:00:00
#SBATCH --ntasks-per-node=48
#SBATCH --array=1-400%20
#SBATCH --job-name=MOOSE_ParamSweep
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vazquezm
#SBATCH -A mat250014p
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err

# ---------------------------------------------------------------
# Parameter Sweep: 20 Young's Moduli x 20 Hardness = 400 tasks
# Array task IDs 1-400 are mapped to (i_ymod, i_hv) index pairs
# using row-major order:
#
#   TASK_ID (0-indexed) = i_ymod * 20 + i_hv
#
#   i_ymod = (TASK_ID) / 20   → selects Young's Modulus row
#   i_hv   = (TASK_ID) % 20   → selects Hardness column
#
# Example:
#   SLURM_ARRAY_TASK_ID=1   → TASK_ID=0   → ymod=100,   Hv=10
#   SLURM_ARRAY_TASK_ID=21  → TASK_ID=20  → ymod=147.4, Hv=10
#   SLURM_ARRAY_TASK_ID=400 → TASK_ID=399 → ymod=1000,  Hv=50
# How to run:
#   Make sure to run "mkdir -p logs rst runs"
#   sbatch 1_sweep.sh
#   sbatch --array=1,400 1_sweep.sh (For just the first and last tasks)
#   sbatch --array=1-400:2 1_sweep.sh (For every other task)
#   sbatch --array=1-20 1_sweep.sh (For one ymod row = 20 Hv points)
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
Hv_nacs=(20 22 23 25 26 28 29 31 33 34 36 37 39 41 42 44 45 47 48 50)         # MPa

TASK_ID=$(( SLURM_ARRAY_TASK_ID - 1 ))
i_ymod=$(( TASK_ID / 20 ))
i_hv=$(( TASK_ID % 20 ))
YMOD=${ymod_nacs[$i_ymod]}
HV=${Hv_nacs[$i_hv]}

TAG="E${YMOD}_H${HV}"

echo "=============================================="
echo "Task $SLURM_ARRAY_TASK_ID -> ymod_nacs=$YMOD MPa, Hv_nacs=$HV MPa  (tag=$TAG)"
echo "=============================================="

# so no two concurrent runs write to the same Exodus / restart files.
RUNDIR="$BASE/runs/$TAG"
rm -rf "$RUNDIR"  # <--- ADD THIS LINE TO AUTO-CLEAN BEFORE RUNNING
mkdir -p "$RUNDIR" logs
cd "$RUNDIR"

# Use SLURM's actual core count rather than assuming the var is set.
NP=${SLURM_NTASKS:-$SLURM_NTASKS_PER_NODE}

# Run from RUNDIR; point MOOSE at the input by absolute path.
# Outputs/file_base is overridden too so the output name is unique
# even if your Outputs block hardcodes a base name.
mpirun -np "$NP" "$EXE" \
    -i "$BASE/1_contact.i" \
    "ymod_nacs=$YMOD" \
    "Hv_nacs=$HV" \
    "Outputs/file_base=${TAG}_out"

EXIT_CODE=$?
echo "MOOSE exit code: $EXIT_CODE"
exit $EXIT_CODE