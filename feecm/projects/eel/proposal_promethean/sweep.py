#!/usr/bin/env python3
"""
Design-map driver for pressure_free_stack.i

Sweeps the three primary design variables identified in the proposal --
the GPE polymerization shrink (eps_GPE), the package shrink (eps_pkg),
and the externally applied stack pressure (sigma_stack) -- and compiles
the post-cure (cycle_amplitude=0) contact-stability metrics into a
single design_map.csv.

Run from the example directory:
    python3 sweep.py

Outputs:
    sweeps/run_NNN.csv       per-run postprocessor history
    sweeps/design_map.csv    one row per sweep point with final-time PPs
"""
from __future__ import annotations

import csv
import itertools
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EEL = HERE / ".." / ".." / "eel-opt"
INPUT = HERE / "pressure_free_stack.i"
OUT_DIR = HERE / "sweeps"

# Sweep grids (proposal ranges).
EPS_GPE_GRID = [0.02, 0.05, 0.08, 0.12]    # GPE polymerization volumetric shrink
EPS_PKG_GRID = [0.02, 0.05, 0.08, 0.10]    # Package coating volumetric shrink
SIGMA_STACK_GRID = [0.0, 0.5, 1.0]         # External stack pressure [MPa]
E_SEP_GRID = [1.0, 10.0, 25.0, 50.0]       # GPE Young's modulus [MPa] (proposal: 1-50)

# Postprocessor columns of interest in the per-run CSV.
INTERFACE_PPS = [
    "pc_cu_anode",
    "pc_anode_sep",
    "pc_sep_cathode",
    "pc_cathode_al",
]
RC_PPS = [
    "Rc_cu_anode",
    "Rc_anode_sep",
    "Rc_sep_cathode",
    "Rc_cathode_al",
]
SUMMARY_PPS = ["pc_avg", "Pi_c"]


def run_one(run_idx: int, eps_GPE: float, eps_pkg: float,
            sigma_stack: float, E_sep: float) -> dict:
    """Invoke eel-opt for one design point and return the final-time PP row."""
    file_base = f"sweeps/run_{run_idx:04d}"
    cmd = [
        str(EEL),
        "-i", str(INPUT),
        f"eps_GPE={eps_GPE}",
        f"eps_pkg={eps_pkg}",
        f"sigma_stack={sigma_stack}",
        f"E_sep={E_sep}",
        "cycle_amplitude=0.0",
        "end_time=1.0",
        f"Outputs/file_base={file_base}",
        "Outputs/exodus=false",
    ]
    print(f"[{run_idx:04d}] eps_GPE={eps_GPE:.3f}  eps_pkg={eps_pkg:.3f}  "
          f"sigma_stack={sigma_stack:.3f}  E_sep={E_sep:.1f}", flush=True)
    proc = subprocess.run(cmd, cwd=HERE, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout[-2000:] + "\n")
        sys.stderr.write(proc.stderr[-2000:] + "\n")
        raise RuntimeError(f"run {run_idx} failed (returncode {proc.returncode})")

    csv_path = HERE / f"{file_base}.csv"
    with open(csv_path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise RuntimeError(f"empty CSV at {csv_path}")
    return rows[-1]


def main() -> int:
    if not EEL.exists():
        sys.stderr.write(f"eel-opt not found at {EEL}\n")
        return 1
    OUT_DIR.mkdir(exist_ok=True)

    grid = list(itertools.product(E_SEP_GRID, EPS_GPE_GRID, EPS_PKG_GRID, SIGMA_STACK_GRID))
    print(f"Sweeping {len(grid)} design points...")

    master_rows = []
    for idx, (E_sep, eps_GPE, eps_pkg, sigma_stack) in enumerate(grid):
        row = run_one(idx, eps_GPE, eps_pkg, sigma_stack, E_sep)

        pcs = [float(row[k]) for k in INTERFACE_PPS]
        rcs = [float(row[k]) for k in RC_PPS]
        record = {
            "run_idx": idx,
            "E_sep_MPa": E_sep,
            "eps_GPE": eps_GPE,
            "eps_pkg": eps_pkg,
            "sigma_stack_MPa": sigma_stack,
            **{k: float(row[k]) for k in INTERFACE_PPS},
            **{k: float(row[k]) for k in RC_PPS},
            "pc_worst_MPa": min(pcs),
            "Rc_worst_mOhm_cm2": max(rcs),
            "pc_avg_MPa": float(row["pc_avg"]),
            "Pi_c": float(row["Pi_c"]),
        }
        master_rows.append(record)
        print(f"      pc_avg={record['pc_avg_MPa']:.3f} MPa  "
              f"pc_worst={record['pc_worst_MPa']:.3f} MPa  "
              f"Rc_worst={record['Rc_worst_mOhm_cm2']:.1f} mOhm.cm2  "
              f"Pi_c={record['Pi_c']:.3f}")

    out_csv = OUT_DIR / "design_map.csv"
    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(master_rows[0].keys()))
        writer.writeheader()
        writer.writerows(master_rows)
    print(f"\nWrote {out_csv}")

    # Console pivot: best Pi_c at sigma_stack=0 for each (E_sep, eps_pkg)
    # at the most-favorable eps_GPE = min(EPS_GPE_GRID).
    print("\nPi_c at sigma_stack=0, eps_GPE={:.2f}:".format(min(EPS_GPE_GRID)))
    print(f"  {'E_sep [MPa]':>11s}  " +
          "  ".join(f"epkg={e:.2f}" for e in EPS_PKG_GRID))
    for E_sep in E_SEP_GRID:
        cells = []
        for eps_pkg in EPS_PKG_GRID:
            rec = next(r for r in master_rows
                       if abs(r["E_sep_MPa"] - E_sep) < 1e-9
                       and abs(r["eps_GPE"] - min(EPS_GPE_GRID)) < 1e-9
                       and abs(r["eps_pkg"] - eps_pkg) < 1e-9
                       and abs(r["sigma_stack_MPa"]) < 1e-9)
            cells.append(f"{rec['Pi_c']:.3f}")
        print(f"  {E_sep:>11.1f}  " + "  ".join(f"{c:>9s}" for c in cells))
    return 0


if __name__ == "__main__":
    sys.exit(main())
