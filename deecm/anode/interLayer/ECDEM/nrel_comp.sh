#!/bin/bash
#SBATCH --job-name="COMP"
#SBATCH --account=interlayer
#SBATCH --time=2:00:00
#SBATCH --ntasks-per-node=94
#SBATCH --nodes=1
#SBATCH --mail-user=jmv8431@rit.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ============================================================
# SLURM job wrapper for parameterised LIGGGHTS compaction
#
# Expected environment variables (set by nrel_run_all.sh):
#   MN    – metal name          (Ag, Mg, Si, Al, Sn)
#   MNUM  – M number            (1-10)
#   VRNUM – volume ratio number (1-10, 15)
#   FP    – fabrication pressure (100, 250, 400)
#   EME   – Young's modulus     [kPa]
#   LABEL – human-readable job label
#   OUTDIR– output directory for logs
# ============================================================

set -euo pipefail

# Set the path to your LIGGGHTS executable
LMP=/kfs3/scratch/jmv8431/LIGGGHTS-PUBLIC/src/lmp_auto

echo "=========================================="
echo "Job:    ${LABEL}"
echo "Metal:  ${MN}  (E = ${EME} kPa)"
echo "Config: M${MNUM} C5 vr${VRNUM}  fp${FP}"
echo "=========================================="

# Run LIGGGHTS, passing all parameters via -var
srun "${LMP}" \
  -var mn    "${MN}"    \
  -var mnum  "${MNUM}"  \
  -var vrnum "${VRNUM}" \
  -var fp    "${FP}"    \
  -var eme   "${EME}"   \
  -in 1_pck.in