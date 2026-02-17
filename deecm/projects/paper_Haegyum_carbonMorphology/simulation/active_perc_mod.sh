#!/bin/bash -l
#SBATCH --job-name=ActivePerc
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=04:00:00
#SBATCH --account=mat250014p
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vazquezm

set -euo pipefail

module load python/3.8.6
module load numpy/1.19.4
# If scipy is not in the default module, pip install it to a venv or use:
# pip install --user scipy

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
RESULT_ROOT="${SCRIPT_DIR}/results"
SUMMARY="${SCRIPT_DIR}/active_summary.txt"

MORPH_DIRS=(mesh_cnt10 mesh_cnt14 mesh_cnt18 mesh_cnt22 mesh_cnt26 mesh_cnt30)

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
    dump="${mr_dir}post/rst_${label}.dump"
    
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