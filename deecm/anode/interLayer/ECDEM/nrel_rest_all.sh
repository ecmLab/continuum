#!/bin/bash
set -euo pipefail

# ============================================================
# Submit LIGGGHTS rest/EC jobs for ALL configurations
#
#   5 metals × 10 M## × 11 vr## × 3 fp = 1,650 jobs
#
# Usage:
#   ./nrel_rest_all.sh            # submit all jobs
#   DRY_RUN=1 ./nrel_rest_all.sh  # preview without submitting
# ============================================================

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
RESULT_ROOT="${SCRIPT_DIR}/results_rest"
DRY_RUN=${DRY_RUN:-0}

# ---------- Fixed simulation parameters ------------------------------------
REST_TIME=30        # Rest duration [s]
DIFF_DT=0.0001      # Diffusion sub-timestep [s]

# ---------- Metal-specific lookup tables -----------------------------------

# Young's modulus [kPa]
declare -A EME_MAP=(
  [Ag]=76.0e6
  [Mg]=44.0e6
  [Si]=112.4e6
  [Al]=68.0e6
  [Sn]=44.3e6
)

# Poisson's ratio [-]
declare -A VME_MAP=(
  [Ag]=0.37
  [Mg]=0.35
  [Si]=0.28
  [Al]=0.36
  [Sn]=0.33
)

# AM self-diffusivity [um^2/s]
declare -A D_AM_MAP=(
  [Ag]=1.0e-10
  [Mg]=1.0e-15
  [Si]=1.0e-16
  [Al]=1.0e-13
  [Sn]=1.0e-18
)

# Cross-pair diffusivity = min(D_AM, 1.0e-14)
# Pre-computed since D_CB = D_LM = 1.0e-14 is constant
declare -A D_CROSS_MAP=(
  [Ag]=1.0e-14
  [Mg]=1.0e-15
  [Si]=1.0e-16
  [Al]=1.0e-14
  [Sn]=1.0e-18
)

# Max Li concentration in AM [mol/m^3]
declare -A C_LI_MAX_AM_MAP=(
  [Ag]=85155.85
  [Mg]=61993.32
  [Si]=90040.34
  [Al]=52631.58
  [Sn]=79199.37
)

# Max volume expansion AM [-]
declare -A V_EXP_MAX_AM_MAP=(
  [Ag]=10.291
  [Mg]=2.25
  [Si]=3.80
  [Al]=1.90
  [Sn]=3.44
)

# Partial molar volume AM
declare -A OMEGA_LI_AM_MAP=(
  [Ag]=1.0602e-05
  [Mg]=8.96e-6
  [Si]=8.18e-6
  [Al]=9e-6
  [Sn]=8.96e-6
)

# ---------- Sweep arrays ---------------------------------------------------

METALS=(Ag Mg Si Al Sn)
MNUMS=(1 2 3 4 5 6 7 8 9 10)
VRNUMS=(1 2 3 4 5 6 7 8 9 10 15)
FPS=(100 250 400)

# ---------- Pre-flight checks ----------------------------------------------

mkdir -p "$RESULT_ROOT"

echo "============================================"
echo "Submitting LIGGGHTS rest/EC jobs"
echo "  Script dir:    $SCRIPT_DIR"
echo "  Results root:  $RESULT_ROOT"
echo "  Rest time:     ${REST_TIME} s"
echo "  Diff dt:       ${DIFF_DT} s"
if [[ "$DRY_RUN" == "1" ]]; then
  echo "  *** DRY RUN – no jobs will be submitted ***"
fi
echo "============================================"

job_count=0
skip_count=0

# ---------- Main loop -------------------------------------------------------

for mn in "${METALS[@]}"; do
  eme="${EME_MAP[$mn]}"
  vme="${VME_MAP[$mn]}"
  d_am="${D_AM_MAP[$mn]}"
  d_cross="${D_CROSS_MAP[$mn]}"
  c_li_max_am="${C_LI_MAX_AM_MAP[$mn]}"
  v_exp_max_am="${V_EXP_MAX_AM_MAP[$mn]}"
  omega_li_am="${OMEGA_LI_AM_MAP[$mn]}"

  for mnum in "${MNUMS[@]}"; do
    for vrnum in "${VRNUMS[@]}"; do
      for fp in "${FPS[@]}"; do

        # The rest script reads the compacted data file produced by 1_pck.in
        datafile="data${mn}/${mn}_M${mnum}C5_vr${vrnum}_fp${fp}.data"
        if [[ ! -f "$datafile" ]]; then
          echo "[warn] Skipping – compacted file not found: ${datafile}" >&2
          skip_count=$((skip_count + 1))
          continue
        fi

        label="${mn}_M${mnum}C5_vr${vrnum}_fp${fp}_rest"
        outdir="${RESULT_ROOT}/${mn}/${label}"
        mkdir -p "$outdir"

        sbatch_cmd=(
          sbatch
          --job-name="${label}"
          --account=interlayer
          --time=8:00:00
          --ntasks-per-node=94
          --nodes=1
          --mail-user=jmv8431@rit.edu
          --mail-type=BEGIN,END,FAIL
          --output="${outdir}/%x_%j.out"
          --error="${outdir}/%x_%j.err"
          --export="ALL,MN=${mn},MNUM=${mnum},VRNUM=${vrnum},FP=${fp},EME=${eme},VME=${vme},D_AM=${d_am},D_CROSS=${d_cross},C_LI_MAX_AM=${c_li_max_am},V_EXP_MAX_AM=${v_exp_max_am},OMEGA_LI_AM=${omega_li_am},REST_TIME=${REST_TIME},DIFF_DT=${DIFF_DT},LABEL=${label},OUTDIR=${outdir}"
          "${SCRIPT_DIR}/nrel_rest.sh"
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
echo "Skipped:   ${skip_count} (missing compacted .data files)"
echo "Monitor with: squeue -u \$USER"
echo "============================================"