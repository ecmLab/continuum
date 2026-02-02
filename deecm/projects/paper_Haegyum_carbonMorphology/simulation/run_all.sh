#!/bin/bash
set -euo pipefail

# Run LIGGGHTS for every generated mesh (SuperP, CNT, Graphene, Homogeneous)
# Usage: LMP_CMD=lmp ./run_all.sh

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)
RESULT_ROOT="${SCRIPT_DIR}/results"
PRESS_TGT=${PRESS_TGT:-2}
LMP_CMD=${LMP_CMD:-lmp}
MORPH_DIRS=(mesh_superp mesh_cnt mesh_graphene mesh_homo)

mkdir -p "$RESULT_ROOT"

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
    echo "[info] Running ${label} using ${datafile}"
    "$LMP_CMD" -in "${SCRIPT_DIR}/lmp_mr.in" \
      -var datafile "$datafile" \
      -var outputdir "$outdir" \
      -var label "$label" \
      -var press_tgt "$PRESS_TGT"
  done
done
