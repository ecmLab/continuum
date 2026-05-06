#!/bin/bash -l
#SBATCH --job-name=ActivePerc
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=03:00:00
#SBATCH --account=mat250014p
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vazquezm

set -euo pipefail

module load anaconda3/2024.10-1

# Run with sbatch active_perc_mod.sh
# This script processes dump files for different morphologies and mass ratios, computes the number of active Ag atoms, and summarizes the results in a text file.
# For (mesh_cnt10 mesh_cnt14 mesh_cnt18 mesh_cnt22 mesh_cnt26 mesh_cnt30) is time=06:00:00

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-.}"
RESULT_ROOT="${SCRIPT_DIR}/results"
SUMMARY="${SCRIPT_DIR}/active_summary.txt"

MORPH_DIRS=(mesh_cnt_surface mesh_homo)
# MORPH_DIRS=(mesh_cnt10 mesh_cnt14 mesh_cnt18 mesh_cnt22 mesh_cnt26 mesh_cnt30)

# Header for summary file
printf "%-25s %-8s  %8s  %8s  %10s\n" "label" "mr" "num_Ag" "active" "frac" > "$SUMMARY"

count=0

for morph in "${MORPH_DIRS[@]}"; do
  morph_results="${RESULT_ROOT}/${morph}"
  
  if [[ ! -d "$morph_results" ]]; then
    echo "[warn] Skipping ${morph}: ${morph_results} not found" >&2
    continue
  fi
  
  for mr_dir in "${morph_results}"/mr*/; do
    [[ -d "$mr_dir" ]] || continue
    
    mr=$(basename "$mr_dir")
    label="${morph}_${mr}"
    dump="${mr_dir}post/rst_pr35_${label}.dump"
    
    if [[ ! -f "$dump" ]]; then
      echo "[warn] No dump file: ${dump}" >&2
      continue
    fi
    
    out_json="${mr_dir}post/active_${label}.json"
    
    echo "[run] ${label}  ->  ${out_json}"
    result=$(python "${SCRIPT_DIR}/active_number_mod.py" "$dump" "$out_json")
    
    # result is: num_Ag  num_active  frac
    printf "%-25s %-8s  %s\n" "$label" "$mr" "$result" >> "$SUMMARY"
    
    count=$((count + 1))
  done
done

echo "============================================"
echo "Processed ${count} dump files"
echo "Summary:  ${SUMMARY}"
echo "============================================"