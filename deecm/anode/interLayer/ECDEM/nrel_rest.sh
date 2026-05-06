#!/bin/bash
#SBATCH --job-name="REST"
#SBATCH --account=interlayer
#SBATCH --time=2:00:00
#SBATCH --ntasks-per-node=94
#SBATCH --nodes=1
#SBATCH --mail-user=jmv8431@rit.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ============================================================
# SLURM job wrapper for parameterised LIGGGHTS rest/EC sim
#
# Expected environment variables (set by nrel_rest_all.sh):
#   MN           – metal name              (Ag, Mg, Si, Al, Sn)
#   MNUM         – M number                (1-10)
#   VRNUM        – volume ratio number     (1-10, 15)
#   FP           – fabrication pressure    (100, 250, 400)
#   EME          – Young's modulus         [kPa]
#   VME          – Poisson's ratio         [-]
#   D_AM         – AM self-diffusivity     [um^2/s]
#   D_CROSS      – cross-pair diffusivity  [um^2/s]
#   C_LI_MAX_AM  – max Li conc. in AM     [mol/m^3]
#   V_EXP_MAX_AM – max volume expansion AM [-]
#   OMEGA_LI_AM  – partial molar vol AM
#   REST_TIME    – rest duration           [s]
#   DIFF_DT      – diffusion sub-timestep  [s]
#   LABEL        – human-readable job label
#   OUTDIR       – output directory for logs
# ============================================================

set -euo pipefail

LMP=/kfs3/scratch/jmv8431/LIGGGHTS-PUBLIC/src/lmp_auto

echo "=========================================="
echo "Job:    ${LABEL}"
echo "Metal:  ${MN}  (E = ${EME} kPa, v = ${VME})"
echo "Config: M${MNUM} C5 vr${VRNUM}  fp${FP}"
echo "Rest:   ${REST_TIME} s,  diff_dt = ${DIFF_DT} s"
echo "D_AM:   ${D_AM},  D_cross: ${D_CROSS}"
echo "=========================================="

# Create dump output directory
mkdir -p "post/postr/rest_${MN}_M${MNUM}C5_vr${VRNUM}_fp${FP}_rest"

srun "${LMP}" \
  -var mn           "${MN}"           \
  -var mnum         "${MNUM}"         \
  -var vrnum        "${VRNUM}"        \
  -var fp           "${FP}"           \
  -var eme          "${EME}"          \
  -var vme          "${VME}"          \
  -var d_am         "${D_AM}"         \
  -var d_cross      "${D_CROSS}"      \
  -var c_li_max_am  "${C_LI_MAX_AM}"  \
  -var v_exp_max_am "${V_EXP_MAX_AM}" \
  -var omega_li_am  "${OMEGA_LI_AM}"  \
  -var rest_time    "${REST_TIME}"    \
  -var diff_dt      "${DIFF_DT}"      \
  -in 2_rest.in