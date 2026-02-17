#!/bin/bash
set -euo pipefail

# Submit LIGGGHTS jobs for all morphologies and mass ratios to SLURM
# Usage: ./run_all.sh
# Dry Run: DRY_RUN=1 ./run_all.sh

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)
RESULT_ROOT="${SCRIPT_DIR}/results"
DRY_RUN=${DRY_RUN:-0}

MORPH_DIRS=(mesh_cnt10 mesh_cnt14 mesh_cnt16 mesh_cnt18 mesh_cnt22 mesh_cnt26 mesh_cnt30)
# MORPH_DIRS=(mesh_superp mesh_cnt mesh_graphene mesh_homo)

mkdir -p "$RESULT_ROOT"

echo "============================================"
echo "Submitting LIGGGHTS jobs to Bridges-2"
echo "  Project root: $PROJECT_ROOT"
echo "  Results:      $RESULT_ROOT"
echo "============================================"

job_count=0

for morph in "${MORPH_DIRS[@]}"; do
  mesh_base="${PROJECT_ROOT}/${morph}/massratio"
  
  if [[ ! -d "$mesh_base" ]]; then
    echo "[warn] Skipping ${morph}: ${mesh_base} not found" >&2
    continue
  fi
  
  for datafile in "${mesh_base}"/mr*/mdl.data; do
    [[ -e "$datafile" ]] || continue
    
    mr_dir=$(basename "$(dirname "$datafile")")
    outdir="${RESULT_ROOT}/${morph}/${mr_dir}"
    label="${morph}_${mr_dir}"
    
    mkdir -p "$outdir"
    
    sbatch_cmd=(
      sbatch
      --job-name="${label}"
      --account=mat250014p
      --partition=RM-shared
      --nodes=1
      --ntasks-per-node=48
      --time=03:30:00
      --output="${outdir}/%x_%j.out"
      --error="${outdir}/%x_%j.err"
      --export="DATAFILE=${datafile},OUTDIR=${outdir},LABEL=${label}"
      "${SCRIPT_DIR}/myjob.sh"
    )
    
    if [[ "$DRY_RUN" == "1" ]]; then
      echo "[dry-run] ${label}"
    else
      echo "[submit] ${label}"
      "${sbatch_cmd[@]}"
    fi
    
    # FIX: Use pre-increment or addition to avoid exit code 1 when job_count=0
    job_count=$((job_count + 1))
  done
done

echo "============================================"
echo "Submitted ${job_count} jobs"
echo "Monitor with: squeue -u \$USER"
echo "============================================"