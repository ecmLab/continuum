#!/bin/bash
set -euo pipefail

# ============================================================
# Submit LIGGGHTS compaction jobs for ALL configurations
#
#   5 metals  × 10 M##  × 11 vr##  × 3 fp  =  1,650 jobs
#
# Usage:
#   ./nrel_run_all.sh            # submit all jobs
#   DRY_RUN=1 ./nrel_run_all.sh  # preview without submitting
# ============================================================

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
RESULT_ROOT="${SCRIPT_DIR}/results"
DRY_RUN=${DRY_RUN:-0}

# ---------- Configuration arrays -------------------------------------------

# Metal name → Young's modulus [kPa]
declare -A EME_MAP=(
  [Ag]=76.0e6
  [Mg]=44.0e6
  [Si]=112.4e6
  [Al]=68.0e6
  [Sn]=44.3e6
)

# Test Parameters
# METALS=(Ag Mg Si Al Sn)
# MNUMS=(1 10)
# VRNUMS=(1 15)
# FPS=(400)

METALS=(Ag)
MNUMS=(1 2 3 4 5 6 7 8 9 10)
VRNUMS=(1 2 3 4 5 6 7 8 9 10 15)
FPS=(100 250 400)

# ---------- Pre-flight checks ----------------------------------------------

mkdir -p "$RESULT_ROOT"

echo "============================================"
echo "Submitting LIGGGHTS compaction jobs"
echo "  Script dir:   $SCRIPT_DIR"
echo "  Results root:  $RESULT_ROOT"
if [[ "$DRY_RUN" == "1" ]]; then
  echo "  *** DRY RUN – no jobs will be submitted ***"
fi
echo "============================================"

job_count=0
skip_count=0

# ---------- Main loop -------------------------------------------------------

for mn in "${METALS[@]}"; do
  eme="${EME_MAP[$mn]}"

  for mnum in "${MNUMS[@]}"; do
    for vrnum in "${VRNUMS[@]}"; do

      # Verify the input .data file exists
      datafile="data${mn}/${mn}_M${mnum}C5_vr${vrnum}.data"
      if [[ ! -f "$datafile" ]]; then
        echo "[warn] Skipping – file not found: ${datafile}" >&2
        skip_count=$((skip_count + 1))
        continue
      fi

      for fp in "${FPS[@]}"; do

        label="${mn}_M${mnum}C5_vr${vrnum}_fp${fp}"
        outdir="${RESULT_ROOT}/${mn}/${label}"
        mkdir -p "$outdir"

        sbatch_cmd=(
          sbatch
          --job-name="${label}"
          --account=interlayer
          --time=2:00:00
          --ntasks-per-node=94
          --nodes=1
          --mail-user=jmv8431@rit.edu
          --mail-type=BEGIN,END,FAIL
          --output="${outdir}/%x_%j.out"
          --error="${outdir}/%x_%j.err"
          --export="ALL,MN=${mn},MNUM=${mnum},VRNUM=${vrnum},FP=${fp},EME=${eme},LABEL=${label},OUTDIR=${outdir}"
          "${SCRIPT_DIR}/nrel_comp.sh"
        )

        if [[ "$DRY_RUN" == "1" ]]; then
          echo "[dry-run] ${label}"
        else
          echo "[submit] ${label}"
          "${sbatch_cmd[@]}"
        fi

        job_count=$((job_count + 1))

      done   # fp
    done     # vrnum
  done       # mnum
done         # mn

echo "============================================"
echo "Submitted: ${job_count} jobs"
echo "Skipped:   ${skip_count} (missing .data files)"
echo "Monitor with: squeue -u \$USER"
echo "============================================"