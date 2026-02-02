#!/bin/bash
# Job name:
#SBATCH --job-name=DRX_C
#
# Partition:
#SBATCH --partition=tier3
#
# Processors:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=03:30:00
#SBATCH --mem=10g


## Commands to run:
spack load lammps@20211027

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
