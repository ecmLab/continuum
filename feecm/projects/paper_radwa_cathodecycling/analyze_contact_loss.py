#!/usr/bin/env python3
"""
analyze_contact_loss.py
========================

Post-process the 3-way ymod_lpsc sweep (see 1_batch.sh / 1_contact.i) run on
the NMC / LPSC quarter-circle contact model and plot:

  1. Interfacial contact loss [%] vs Time [h]
  2. End-of-cycle contact loss [%] vs Cycle Number [-]
  3. NMC volume (area) expansion [%] vs Time [h]
  4. Porosity [%] vs Time [h]
  5. End-of-cycle porosity [%] vs Cycle Number [-]

The 1_contact.i input only writes an Exodus file (no Postprocessors/CSV), so
every quantity below is derived directly from the Exodus nodal fields that
MOOSE's [Contact] action and mesh geometry already provide:

  * contact_pressure  (nodal, on the secondary surface 'block_LPSC_left')
        0  -> node NOT in contact  (used for contact loss %)
        >0 -> node in contact
  * penetration        (nodal, on the secondary surface)
        <0 -> local separation distance (gap) between NMC and LPSC
        >=0 -> overlap / in contact
        Used to reconstruct the local gap width, which is integrated along
        the (reference) interface arc-length to get the void ("blank
        space") area opened up between the two bodies.
  * coordx/coordy + disp_x/disp_y -> deformed nodal positions, used to
        (a) sum current element areas over block_NMC for the volume
            expansion metric, and
        (b) trace the outer domain boundary (block_bottom/right/top/left)
            with the shoelace formula to get the total available area for
            the porosity ratio.

Cycle boundaries are assumed at t = k * CYCLE_PERIOD_H (k = 1, 2, ...),
matching cycle_period=1.0 in 1_contact.i (one hour per charge/discharge
cycle); "end of cycle" = the synced output closest to that time.

Requires: numpy, matplotlib, netCDF4 (`pip install netCDF4`; already installed
in the `moose` conda env used to build/run this project).

Usage
-----
    python analyze_contact_loss.py [--root DIR]

By default DIR is this script's directory, and Exodus files are looked up
at  <DIR>/runs/<TAG>/<TAG>_out.e*  for TAG in {E370_LPSC, E550_LPSC,
E22000_LPSC}, matching the paths written by 1_batch.sh.
"""

import argparse
import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

# ---------------------------------------------------------------------------
# Run configuration (must match 1_batch.sh / 1_contact.i)
# ---------------------------------------------------------------------------
RUNS = [
    dict(ymod=370, tag="E370_LPSC", label="370 MPa", color="black"),
    dict(ymod=550, tag="E550_LPSC", label="550 MPa", color="red"),
    dict(ymod=22000, tag="E22000_LPSC", label="22 GPa", color="blue"),
]

CYCLE_PERIOD_H = 1.0  # cycle_period in 1_contact.i, interpreted as hours
CONTACT_PRESSURE_TOL = 1e-8  # MPa; <= this counts as "not in contact"

SIDESET_NMC_PRIMARY = "block_NMC_right"
SIDESET_LPSC_SECONDARY = "block_LPSC_left"
BOUNDARY_BOTTOM = "block_bottom"
BOUNDARY_RIGHT = "block_right"
BOUNDARY_TOP = "block_top"
BOUNDARY_LEFT = "block_left"
BLOCK_NMC = "block_NMC"


# ---------------------------------------------------------------------------
# Low-level Exodus (NetCDF) helpers
# ---------------------------------------------------------------------------
def _decode_name_array(var):
    """Decode an Exodus '...names' char array into a list of Python strings."""
    names = []
    for row in var[:]:
        raw = row.filled(b"\x00") if hasattr(row, "filled") else row
        names.append(bytes(raw).decode("utf-8", errors="ignore").strip("\x00").strip())
    return names


def find_exodus_files(root, tag):
    """Find Exodus file(s) for a run tag, in the layout written by 1_batch.sh."""
    pattern = os.path.join(root, "runs", tag, f"{tag}_out.e*")
    files = sorted(glob.glob(pattern))
    if not files:
        # Fall back to a recursive search anywhere under root.
        pattern = os.path.join(root, "**", f"{tag}_out.e*")
        files = sorted(glob.glob(pattern, recursive=True))
    return files


class RunData:
    """Container for the per-time-step series computed from one Exodus run."""

    def __init__(self, label, color):
        self.label = label
        self.color = color
        self.time_h = np.empty(0)
        self.contact_loss_pct = np.empty(0)
        self.porosity_pct = np.empty(0)
        self.nmc_expansion_pct = np.empty(0)


def _nodeset_lookup(ds, ns_names):
    """Return dict: boundary name -> 0-based node indices into coordx/coordy."""
    lookup = {}
    for i, name in enumerate(ns_names):
        if not name:
            continue
        if name in lookup:
            continue  # keep first occurrence (MOOSE may duplicate empty-named sets)
        lookup[name] = ds.variables[f"node_ns{i + 1}"][:].astype(np.int64) - 1
    return lookup


def _nodal_var_index(name_nod_var, name):
    return name_nod_var.index(name) + 1  # vals_nod_var<k> is 1-indexed


def _quad_areas(x, y):
    """Vectorized shoelace area for an (N,4) array of QUAD4 corner coords."""
    return 0.5 * np.abs(
        (x[:, 0] * y[:, 1] - x[:, 1] * y[:, 0])
        + (x[:, 1] * y[:, 2] - x[:, 2] * y[:, 1])
        + (x[:, 2] * y[:, 3] - x[:, 3] * y[:, 2])
        + (x[:, 3] * y[:, 0] - x[:, 0] * y[:, 3])
    )


def _polygon_area(x, y):
    """Shoelace area of a closed polygon given ordered vertex coordinates."""
    return 0.5 * np.abs(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _ordered_boundary_loop(coordx, coordy, ns_lookup):
    """Node index loop tracing the outer rectangle CCW using reference coords."""
    bottom = ns_lookup[BOUNDARY_BOTTOM]
    right = ns_lookup[BOUNDARY_RIGHT]
    top = ns_lookup[BOUNDARY_TOP]
    left = ns_lookup[BOUNDARY_LEFT]

    bottom = bottom[np.argsort(coordx[bottom])]  # (0,0) -> (L,0)
    right = right[np.argsort(coordy[right])]  # (L,0) -> (L,L)
    top = top[np.argsort(-coordx[top])]  # (L,L) -> (0,L)
    left = left[np.argsort(-coordy[left])]  # (0,L) -> (0,0)
    return np.concatenate([bottom, right, top, left])


def _arc_segment_weights(coordx, coordy, node_idx):
    """Sort interface nodes by reference angle and return (order, ds weights).

    ds weights are trapezoidal (half of each adjacent chord) reference
    arc-length weights that sum to the total reference interface length --
    used to turn a per-node contact/gap signal into an arc-length-weighted
    quantity.
    """
    ang = np.arctan2(coordy[node_idx], coordx[node_idx])
    order = np.argsort(ang)
    idx_sorted = node_idx[order]
    x = coordx[idx_sorted]
    y = coordy[idx_sorted]
    chord = np.hypot(np.diff(x), np.diff(y))  # length of each segment i -> i+1
    return idx_sorted, chord


def analyze_run(files, label, color):
    """Read one run's Exodus file(s) and compute the derived time series."""
    times_all = []
    contact_all = []
    porosity_all = []
    expansion_all = []

    a0_nmc = None  # reference (t=0) NMC area, carried across restart files

    for fpath in files:
        ds = nc.Dataset(fpath, "r")
        try:
            ns_names = _decode_name_array(ds.variables["ns_names"])
            eb_names = _decode_name_array(ds.variables["eb_names"])
            name_nod_var = _decode_name_array(ds.variables["name_nod_var"])

            ns_lookup = _nodeset_lookup(ds, ns_names)
            for req in (SIDESET_NMC_PRIMARY, SIDESET_LPSC_SECONDARY,
                        BOUNDARY_BOTTOM, BOUNDARY_RIGHT, BOUNDARY_TOP, BOUNDARY_LEFT):
                if req not in ns_lookup:
                    raise RuntimeError(f"{fpath}: node set '{req}' not found in Exodus file")

            idx_cp = _nodal_var_index(name_nod_var, "contact_pressure")
            idx_pen = _nodal_var_index(name_nod_var, "penetration")
            idx_dx = _nodal_var_index(name_nod_var, "disp_x")
            idx_dy = _nodal_var_index(name_nod_var, "disp_y")

            coordx = ds.variables["coordx"][:].astype(np.float64)
            coordy = ds.variables["coordy"][:].astype(np.float64)

            sec_nodes = ns_lookup[SIDESET_LPSC_SECONDARY]
            sec_idx_sorted, sec_ds = _arc_segment_weights(coordx, coordy, sec_nodes)
            total_arc_len = sec_ds.sum()
            # Trapezoidal node weight = half of each adjacent segment.
            node_w = np.zeros(len(sec_idx_sorted))
            node_w[:-1] += 0.5 * sec_ds
            node_w[1:] += 0.5 * sec_ds

            boundary_loop = _ordered_boundary_loop(coordx, coordy, ns_lookup)

            nmc_block_idx = eb_names.index(BLOCK_NMC) + 1
            conn_nmc = ds.variables[f"connect{nmc_block_idx}"][:].astype(np.int64) - 1

            cp_var = ds.variables[f"vals_nod_var{idx_cp}"]
            pen_var = ds.variables[f"vals_nod_var{idx_pen}"]
            dx_var = ds.variables[f"vals_nod_var{idx_dx}"]
            dy_var = ds.variables[f"vals_nod_var{idx_dy}"]
            time_whole = ds.variables["time_whole"][:].astype(np.float64)

            n_steps = len(time_whole)
            contact_loss = np.empty(n_steps)
            porosity = np.empty(n_steps)
            expansion = np.empty(n_steps)

            for t in range(n_steps):
                dx = dx_var[t, :]
                dy = dy_var[t, :]

                # --- contact loss % (arc-length weighted) ---
                cp = cp_var[t, sec_idx_sorted]
                not_in_contact = cp <= CONTACT_PRESSURE_TOL
                contact_loss[t] = 100.0 * node_w[not_in_contact].sum() / total_arc_len

                # --- porosity % = blank (void) area / total domain area ---
                pen = pen_var[t, sec_idx_sorted]
                gap = np.clip(-pen, 0.0, None)  # local separation distance
                blank_area = np.sum(0.5 * (gap[:-1] + gap[1:]) * sec_ds)
                bx = coordx[boundary_loop] + dx[boundary_loop]
                by = coordy[boundary_loop] + dy[boundary_loop]
                total_area = _polygon_area(bx, by)
                porosity[t] = 100.0 * blank_area / total_area

                # --- NMC area (volume) expansion % ---
                xe = coordx[conn_nmc] + dx[conn_nmc]
                ye = coordy[conn_nmc] + dy[conn_nmc]
                a_nmc = _quad_areas(xe, ye).sum()
                if a0_nmc is None:
                    a0_nmc = a_nmc
                expansion[t] = 100.0 * (a_nmc - a0_nmc) / a0_nmc

            times_all.append(time_whole)
            contact_all.append(contact_loss)
            porosity_all.append(porosity)
            expansion_all.append(expansion)
            print(f"  [{label}] {os.path.basename(fpath)}: {n_steps} time steps "
                  f"(t = {time_whole[0]:.3f} -> {time_whole[-1]:.3f} h)")
        finally:
            ds.close()

    time_h = np.concatenate(times_all)
    contact_loss_pct = np.concatenate(contact_all)
    porosity_pct = np.concatenate(porosity_all)
    nmc_expansion_pct = np.concatenate(expansion_all)

    # Merge restart files: sort by time and drop duplicate/overlapping times,
    # keeping the LAST occurrence (i.e. the most recently written value).
    order = np.argsort(time_h, kind="stable")
    time_h = time_h[order]
    contact_loss_pct = contact_loss_pct[order]
    porosity_pct = porosity_pct[order]
    nmc_expansion_pct = nmc_expansion_pct[order]

    _, keep_idx = np.unique(time_h[::-1], return_index=True)
    keep_idx = len(time_h) - 1 - keep_idx
    keep_idx.sort()

    run = RunData(label, color)
    run.time_h = time_h[keep_idx]
    run.contact_loss_pct = contact_loss_pct[keep_idx]
    run.porosity_pct = porosity_pct[keep_idx]
    run.nmc_expansion_pct = nmc_expansion_pct[keep_idx]
    return run


def end_of_cycle_series(run):
    """Sample a run's series at the synced output closest to each cycle end."""
    n_cycles = int(np.floor(run.time_h.max() / CYCLE_PERIOD_H + 1e-6))
    cycles = np.arange(1, n_cycles + 1)
    target_times = cycles * CYCLE_PERIOD_H
    idx = np.searchsorted(run.time_h, target_times)
    idx = np.clip(idx, 0, len(run.time_h) - 1)
    # searchsorted may land one step early/late; pick the closer neighbor.
    left = np.clip(idx - 1, 0, len(run.time_h) - 1)
    use_left = np.abs(run.time_h[left] - target_times) < np.abs(run.time_h[idx] - target_times)
    idx = np.where(use_left, left, idx)
    return cycles, idx


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def apply_style():
    plt.rcParams.update({
        "font.family": "Times New Roman",
        "font.size": 18,
        "axes.linewidth": 1.5,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "axes.labelsize": 24,
        "legend.fontsize": 16,
        "lines.linewidth": 1.5,
        "lines.linestyle": "-",
    })


def new_figure():
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_frame_on(True)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.5)
    return fig, ax


def finish_and_save(fig, ax, xlabel, ylabel, out_dir, basename):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="both", labelsize=18)
    legend = ax.legend(loc="best", fontsize=16, frameon=False)
    fig.tight_layout()
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(out_dir, f"{basename}.{ext}"), dpi=300)
    plt.close(fig)


def plot_vs_time(runs, attr, ylabel, out_dir, basename):
    fig, ax = new_figure()
    for run in runs:
        ax.plot(run.time_h, getattr(run, attr), color=run.color, label=run.label,
                linestyle="-", linewidth=1.5)
    finish_and_save(fig, ax, "Time [h]", ylabel, out_dir, basename)


def plot_end_of_cycle(runs, attr, ylabel, out_dir, basename):
    fig, ax = new_figure()
    for run in runs:
        cycles, idx = end_of_cycle_series(run)
        ax.plot(cycles, getattr(run, attr)[idx], color=run.color, label=run.label,
                linestyle="-", linewidth=1.5)
    finish_and_save(fig, ax, "Cycle Number [-]", ylabel, out_dir, basename)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", default=os.path.dirname(os.path.abspath(__file__)),
                         help="Project root containing the 'runs' directory (default: script dir)")
    args = parser.parse_args()

    out_dir = os.path.join(args.root, "plot")
    os.makedirs(out_dir, exist_ok=True)

    runs = []
    for cfg in RUNS:
        files = find_exodus_files(args.root, cfg["tag"])
        if not files:
            print(f"WARNING: no Exodus output found for tag '{cfg['tag']}' "
                  f"(expected under runs/{cfg['tag']}/{cfg['tag']}_out.e*) -- skipping.",
                  file=sys.stderr)
            continue
        print(f"Analyzing {cfg['label']} (tag={cfg['tag']}) ...")
        run = analyze_run(files, cfg["label"], cfg["color"])
        runs.append(run)

    if not runs:
        print("No runs found -- nothing to plot. Copy the HPC 'runs/' directory "
              "next to this script (or pass --root) and re-run.", file=sys.stderr)
        sys.exit(1)

    apply_style()

    plot_vs_time(runs, "contact_loss_pct", "Interfacial Contact Loss [%]",
                 out_dir, "1_contact_loss_vs_time")
    plot_end_of_cycle(runs, "contact_loss_pct", "End of Cycle Contact Loss [%]",
                       out_dir, "2_end_of_cycle_contact_loss")
    plot_vs_time(runs, "nmc_expansion_pct", "NMC Volume Expansion [%]",
                 out_dir, "3_nmc_volume_expansion_vs_time")
    plot_vs_time(runs, "porosity_pct", "Porosity [%]",
                 out_dir, "4_porosity_vs_time")
    plot_end_of_cycle(runs, "porosity_pct", "End of Cycle Porosity [%]",
                       out_dir, "5_end_of_cycle_porosity")

    print(f"\nSaved 5 figures (.svg + .png) to: {out_dir}")


if __name__ == "__main__":
    main()