#!/bin/bash
#SBATCH --job-name=mat_60
#SBATCH --account=interlayer
#SBATCH --partition=shared
#SBATCH --qos=normal
#SBATCH --array=1-11
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --mail-user=jmv8431@rit.edu
#SBATCH --mail-type=END,FAIL,ARRAY_TASKS

# ------------------------------------------------------------
# Per-task AU math (Kestrel shared CPU, charge factor = 10):
#   fraction = max(1/104, 2/240) = 0.00962  (cores dominate)
#   AU / hr  = 0.0962
#
# Skewed example: 1 case x 2 h + 10 cases x 10 min:
#   AUs = 0.0962 * (2 + 10*(10/60)) = 0.36 AU
# Same workload as a 12-core parfor job: ~2.3 AU.  ~6x cheaper.
#
# Also: tasks run/fail independently. Re-run only the broken
# ones with e.g.  sbatch --array=3,7 run_matlab.sh
# ------------------------------------------------------------

mkdir -p logs massratio

module purge
module load matlab

cd "$SLURM_SUBMIT_DIR"

# Each MATLAB instance gets one core; pin internal threading to it.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Per-task tmp prefs dir keeps concurrent MATLABs from clobbering each other
export MATLAB_PREFDIR="${TMPDIR:-/tmp}/matlab_prefs_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$MATLAB_PREFDIR"

echo "[$(date)] Starting array task ${SLURM_ARRAY_TASK_ID} on $(hostname -s)"
matlab -nodisplay -nosplash -nodesktop -batch "mdl_run"
echo "[$(date)] Finished array task ${SLURM_ARRAY_TASK_ID}"