#!/bin/bash
#SBATCH --job-name=DRX_C
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=03:30:00
#SBATCH --account=mat250014p
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vazquezm

set -euo pipefail

# Load required modules LIGGGHTS
module load gcc/10.2.0 openmpi/5.0.3-gcc13.2.1 cuda/11.7.1

if [[ -z "${DATAFILE:-}" ]]; then
  echo "ERROR: Set DATAFILE to the mesh data file (e.g., /path/to/mr1/mdl.data)." >&2
  exit 1
fi

OUTDIR=${OUTDIR:-results/job}
LABEL=${LABEL:-job}

mkdir -p "$OUTDIR"

mpiexec -np $SLURM_NTASKS /ocean/projects/mat250014p/shared/software/LIGGGHTS-PUBLIC/src/lmp_auto -in lmp_mr.in \
  -var datafile "$DATAFILE" \
  -var outputdir "$OUTDIR" \
  -var label "$LABEL"