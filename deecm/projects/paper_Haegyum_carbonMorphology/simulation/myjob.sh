#!/bin/bash
#SBATCH --job-name=DRX_C
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=03:30:00
#SBATCH --account=mat250014p
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vazquezm

## Commands to run:
## spack load lammps@20211027

set -euo pipefail

module load openmpi/5.0.3-gcc13.2.1 cuda/11.7.1 mkl/2020.4.304 python/3.8.6 LAMMPS/29Aug24-gnu

if [[ -z "${DATAFILE:-}" ]]; then
  echo "ERROR: Set DATAFILE to the mesh data file (e.g., /path/to/mr1/mdl.data)." >&2
  exit 1
fi

OUTDIR=${OUTDIR:-results/job}
LABEL=${LABEL:-job}
PRESS_TGT=${PRESS_TGT:-2}
LMP_CMD=${LMP_CMD:-lmp}

mkdir -p "$OUTDIR"

srun "$LMP_CMD" -in lmp_mr.in \
  -var datafile "$DATAFILE" \
  -var outputdir "$OUTDIR" \
  -var label "$LABEL" \
  -var press_tgt "$PRESS_TGT"
