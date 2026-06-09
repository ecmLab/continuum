#!/usr/bin/env python3
"""
Generate figures for design_principles.tex.

Reads:
  ../sweeps/design_map.csv     (from sweep.py)
  ../rst/cycling_demo.csv      (from running pressure_free_stack.i with cycle_amplitude=1)

Writes (vector PDF, into figures/):
  geometry_schematic.pdf       Labeled cross-section of the wrapped pouch
  design_map_heatmap.pdf       Pi_c heatmaps over (eps_GPE, eps_pkg) at three sigma_stack
  Pi_c_vs_eps_pkg.pdf          Pi_c vs eps_pkg, lines = eps_GPE, panels = sigma_stack
  Pi_c_vs_eps_GPE.pdf          Pi_c vs eps_GPE, lines = eps_pkg, panels = sigma_stack
  interface_pc_bars.pdf        Per-interface pc, baseline + 3 chosen design points
  cycling_trajectory.pdf       Pi_c(t) over the 2-cycle Li-breathing run
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

HERE = Path(__file__).resolve().parent
EXAMPLE_DIR = HERE.parent
FIG_DIR = HERE / "figures"
FIG_DIR.mkdir(exist_ok=True)

DESIGN_CSV = EXAMPLE_DIR / "sweeps" / "design_map.csv"
CYCLE_CSV = EXAMPLE_DIR / "rst" / "cycling_demo.csv"
T_GRID_CSV = EXAMPLE_DIR / "sweeps" / "T_grid.csv"


# -----------------------------------------------------------------------------
# Figure 1: geometry schematic
# -----------------------------------------------------------------------------
def fig_geometry():
    # Layer thicknesses [mm] (must match pressure_free_stack.i).
    t_pkg = 0.050
    t_cu = 0.010
    t_anode = 0.080
    t_sep = 0.025
    t_cathode = 0.080
    t_al = 0.015
    W_inner = 1.0
    W_outer = W_inner + 2 * t_pkg

    layers = [
        ("Cu collector", t_cu, "#d2691e"),
        ("Graphite anode", t_anode, "#404040"),
        ("Separator (GPE)", t_sep, "#7f9fbf"),
        ("LVPF cathode", t_cathode, "#1a4f8a"),
        ("Al collector", t_al, "#a8a8a8"),
    ]
    H_inner = sum(t for _, t, _ in layers)
    H_outer = H_inner + 2 * t_pkg
    pkg_color = "#d8b673"

    fig, ax = plt.subplots(figsize=(7.5, 4.6))

    # Outer package envelope (filled background)
    ax.add_patch(patches.Rectangle((0, 0), W_outer, H_outer,
                                   facecolor=pkg_color, edgecolor="black", linewidth=1.0))

    # Active layers (in inner column)
    y = t_pkg
    for label, t, color in layers:
        ax.add_patch(patches.Rectangle((t_pkg, y), W_inner, t,
                                       facecolor=color, edgecolor="black", linewidth=0.6))
        ax.text(t_pkg + W_inner / 2, y + t / 2, label,
                ha="center", va="center", fontsize=9, color="white", weight="bold")
        y += t

    # Hatch the inner column boundary so package wrap is obvious
    for x_left in (0, t_pkg + W_inner):
        ax.add_patch(patches.Rectangle((x_left, t_pkg), t_pkg, H_inner,
                                       facecolor=pkg_color, edgecolor="black",
                                       linewidth=0.6, hatch="///", alpha=0.7))
    # And top/bottom strips with same hatch
    for y_strip in (0, H_outer - t_pkg):
        ax.add_patch(patches.Rectangle((0, y_strip), W_outer, t_pkg,
                                       facecolor=pkg_color, edgecolor="black",
                                       linewidth=0.6, hatch="///", alpha=0.7))

    # Eigenstrain callouts (shrink arrows from outside)
    ax.annotate("$\\varepsilon_{\\rm pkg}$ shrink (envelope)",
                xy=(W_outer * 0.30, H_outer - t_pkg / 2),
                xytext=(W_outer * 0.30, H_outer + 0.06),
                ha="center", fontsize=10, color="firebrick",
                arrowprops=dict(arrowstyle="->", color="firebrick", lw=1.0))
    ax.annotate("$\\varepsilon_{\\rm GPE}$ shrink",
                xy=(t_pkg + W_inner * 0.85, t_pkg + t_cu + t_anode + t_sep / 2),
                xytext=(W_outer + 0.06, t_pkg + t_cu + t_anode + t_sep / 2),
                ha="left", va="center", fontsize=10, color="firebrick",
                arrowprops=dict(arrowstyle="->", color="firebrick", lw=1.0))

    # External stack pressure (top)
    for x_arrow in np.linspace(W_outer * 0.55, W_outer * 0.92, 4):
        ax.annotate("", xy=(x_arrow, H_outer), xytext=(x_arrow, H_outer + 0.03),
                    arrowprops=dict(arrowstyle="-|>", color="darkgreen", lw=0.9))
    ax.text(W_outer * 0.74, H_outer + 0.06,
            "$\\sigma_{\\rm stack}$ (applied)",
            ha="center", fontsize=10, color="darkgreen")

    # Left-edge roller
    for y_arrow in np.linspace(t_pkg * 0.5, H_outer - t_pkg * 0.5, 5):
        ax.annotate("", xy=(0, y_arrow), xytext=(-0.05, y_arrow),
                    arrowprops=dict(arrowstyle="-|>", color="steelblue", lw=0.9))
    ax.text(-0.07, H_outer / 2, "$u_x = 0$\n(symmetry roller)",
            ha="right", va="center", fontsize=8, color="steelblue")

    # Bottom roller
    for x_arrow in np.linspace(W_outer * 0.15, W_outer * 0.85, 5):
        ax.annotate("", xy=(x_arrow, 0), xytext=(x_arrow, -0.03),
                    arrowprops=dict(arrowstyle="-|>", color="steelblue", lw=0.9))
    ax.text(W_outer / 2, -0.06, "$u_y = 0$ (roller)",
            ha="center", va="top", fontsize=8, color="steelblue")

    ax.set_xlim(-0.32, W_outer + 0.42)
    ax.set_ylim(-0.13, H_outer + 0.13)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("Wrapped-pouch geometry (2D plane strain)", fontsize=11)

    fig.tight_layout()
    out = FIG_DIR / "geometry_schematic.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


# -----------------------------------------------------------------------------
# Read sweep CSV
# -----------------------------------------------------------------------------
DEFAULT_E_SEP = 10.0   # MPa, default GPE modulus for the (eps_GPE, eps_pkg, sigma) figures


def load_design_map(E_sep_filter: float | None = DEFAULT_E_SEP):
    """Load the master CSV and optionally filter to a single E_sep slice."""
    if not DESIGN_CSV.exists():
        sys.stderr.write(f"missing {DESIGN_CSV}; run sweep.py first\n")
        sys.exit(1)
    rows = []
    with open(DESIGN_CSV, newline="") as fh:
        for r in csv.DictReader(fh):
            rows.append({k: float(v) for k, v in r.items()})
    if E_sep_filter is not None and "E_sep_MPa" in rows[0]:
        rows = [r for r in rows if abs(r["E_sep_MPa"] - E_sep_filter) < 1e-9]
    eps_GPE = sorted({r["eps_GPE"] for r in rows})
    eps_pkg = sorted({r["eps_pkg"] for r in rows})
    sigma = sorted({r["sigma_stack_MPa"] for r in rows})
    return rows, eps_GPE, eps_pkg, sigma


def load_design_map_full():
    """Load the master CSV with no E_sep filter (for the E_sep lever figure)."""
    return load_design_map(E_sep_filter=None)


def pivot(rows, key, eps_GPE, eps_pkg, sigma_value):
    M = np.full((len(eps_GPE), len(eps_pkg)), np.nan)
    for r in rows:
        if abs(r["sigma_stack_MPa"] - sigma_value) > 1e-9:
            continue
        i = eps_GPE.index(r["eps_GPE"])
        j = eps_pkg.index(r["eps_pkg"])
        M[i, j] = r[key]
    return M


# -----------------------------------------------------------------------------
# Figure 2: heatmaps of Pi_c
# -----------------------------------------------------------------------------
def fig_heatmaps():
    rows, eps_GPE, eps_pkg, sigma = load_design_map()
    fig, axes = plt.subplots(1, len(sigma), figsize=(4.2 * len(sigma), 3.8),
                             sharey=True, constrained_layout=True)
    if len(sigma) == 1:
        axes = [axes]
    vmin = min(r["Pi_c"] for r in rows)
    vmax = max(r["Pi_c"] for r in rows)
    # Build cell-edge extents from data centers (half-step beyond first/last sample).
    def edges(xs):
        a = np.asarray(xs, dtype=float)
        d = np.diff(a)
        e = np.concatenate(([a[0] - d[0] / 2], a[:-1] + d / 2, [a[-1] + d[-1] / 2]))
        return float(e[0]), float(e[-1])
    x_lo, x_hi = edges(eps_pkg)
    y_lo, y_hi = edges(eps_GPE)
    for ax, s in zip(axes, sigma):
        M = pivot(rows, "Pi_c", eps_GPE, eps_pkg, s)
        im = ax.imshow(M, origin="lower", aspect="auto",
                       extent=[x_lo, x_hi, y_lo, y_hi],
                       vmin=vmin, vmax=vmax, cmap="viridis")
        # Pi_c = 1 boundary contour
        try:
            X, Y = np.meshgrid(eps_pkg, eps_GPE)
            cs = ax.contour(X, Y, M, levels=[1.0], colors="white",
                            linewidths=1.6, linestyles="-")
            ax.clabel(cs, inline=True, fontsize=7, fmt="$\\Pi_c=1$")
        except Exception:
            pass
        # overlay numbers at cell centers
        for i, eg in enumerate(eps_GPE):
            for j, ep in enumerate(eps_pkg):
                color = "white" if M[i, j] < (vmin + vmax) * 0.55 else "black"
                ax.text(ep, eg, f"{M[i,j]:.2f}",
                        ha="center", va="center", fontsize=7, color=color)
        ax.set_xlabel("$\\varepsilon_{\\rm pkg}$")
        ax.set_title(f"$\\sigma_{{\\rm stack}}={s:.2f}$ MPa", fontsize=10)
        ax.set_xticks(eps_pkg)
        ax.set_yticks(eps_GPE)
    axes[0].set_ylabel("$\\varepsilon_{\\rm GPE}$")
    cbar = fig.colorbar(im, ax=axes, shrink=0.85, label="$\\Pi_c$")
    fig.suptitle("Contact stability $\\Pi_c$ across the design grid", fontsize=11)
    out = FIG_DIR / "design_map_heatmap.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


# -----------------------------------------------------------------------------
# Figure 3: Pi_c vs eps_pkg (sigma panels, eps_GPE lines)
# -----------------------------------------------------------------------------
def fig_lines_eps_pkg():
    rows, eps_GPE, eps_pkg, sigma = load_design_map()
    fig, axes = plt.subplots(1, len(sigma), figsize=(4.0 * len(sigma), 3.6),
                             sharey=True, constrained_layout=True)
    if len(sigma) == 1:
        axes = [axes]
    cmap = plt.get_cmap("plasma")
    for ax, s in zip(axes, sigma):
        for k, eg in enumerate(eps_GPE):
            ys = [next(r["Pi_c"] for r in rows
                       if abs(r["sigma_stack_MPa"] - s) < 1e-9
                       and abs(r["eps_GPE"] - eg) < 1e-9
                       and abs(r["eps_pkg"] - ep) < 1e-9) for ep in eps_pkg]
            ax.plot(eps_pkg, ys, "-o", color=cmap(k / max(1, len(eps_GPE) - 1)),
                    label=f"$\\varepsilon_{{\\rm GPE}}={eg:.2f}$", markersize=4)
        ax.axhline(1.0, color="black", lw=0.8, ls="--", alpha=0.6)
        ax.set_xlabel("$\\varepsilon_{\\rm pkg}$")
        ax.set_title(f"$\\sigma_{{\\rm stack}}={s:.2f}$ MPa", fontsize=10)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("$\\Pi_c$")
    axes[-1].legend(fontsize=8, loc="upper left")
    fig.suptitle("$\\Pi_c$ scales nearly linearly with package shrink "
                 "(dashed line: $\\Pi_c=1$ target)",
                 fontsize=10)
    out = FIG_DIR / "Pi_c_vs_eps_pkg.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


# -----------------------------------------------------------------------------
# Figure 4: Pi_c vs eps_GPE (sigma panels, eps_pkg lines)
# -----------------------------------------------------------------------------
def fig_lines_eps_GPE():
    rows, eps_GPE, eps_pkg, sigma = load_design_map()
    fig, axes = plt.subplots(1, len(sigma), figsize=(4.0 * len(sigma), 3.6),
                             sharey=True, constrained_layout=True)
    if len(sigma) == 1:
        axes = [axes]
    cmap = plt.get_cmap("viridis")
    for ax, s in zip(axes, sigma):
        for k, ep in enumerate(eps_pkg):
            ys = [next(r["Pi_c"] for r in rows
                       if abs(r["sigma_stack_MPa"] - s) < 1e-9
                       and abs(r["eps_GPE"] - eg) < 1e-9
                       and abs(r["eps_pkg"] - ep) < 1e-9) for eg in eps_GPE]
            ax.plot(eps_GPE, ys, "-s", color=cmap(k / max(1, len(eps_pkg) - 1)),
                    label=f"$\\varepsilon_{{\\rm pkg}}={ep:.2f}$", markersize=4)
        ax.axhline(1.0, color="black", lw=0.8, ls="--", alpha=0.6)
        ax.set_xlabel("$\\varepsilon_{\\rm GPE}$")
        ax.set_title(f"$\\sigma_{{\\rm stack}}={s:.2f}$ MPa", fontsize=10)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("$\\Pi_c$")
    axes[-1].legend(fontsize=8, loc="upper right")
    fig.suptitle("Increasing $\\varepsilon_{\\rm GPE}$ slightly reduces $\\Pi_c$ "
                 "(separator unloading)", fontsize=10)
    out = FIG_DIR / "Pi_c_vs_eps_GPE.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


# -----------------------------------------------------------------------------
# Figure 5: per-interface pc bar chart, several design points
# -----------------------------------------------------------------------------
def fig_interface_bars():
    rows, _, _, _ = load_design_map()

    def pick(eg, ep, sg):
        return next(r for r in rows
                    if abs(r["eps_GPE"] - eg) < 1e-9
                    and abs(r["eps_pkg"] - ep) < 1e-9
                    and abs(r["sigma_stack_MPa"] - sg) < 1e-9
                    and abs(r.get("E_sep_MPa", DEFAULT_E_SEP) - DEFAULT_E_SEP) < 1e-9)

    selections = [
        ("baseline\n(0.05/0.05/0)", pick(0.05, 0.05, 0.0)),
        ("low pkg, high GPE\n(0.12/0.02/0)", pick(0.12, 0.02, 0.0)),
        ("high pkg\n(0.05/0.10/0)", pick(0.05, 0.10, 0.0)),
        ("high pkg + 1 MPa\n(0.05/0.10/1.0)", pick(0.05, 0.10, 1.0)),
    ]
    pp_keys = ["pc_cu_anode", "pc_anode_sep", "pc_sep_cathode", "pc_cathode_al"]
    pp_labels = ["cu|anode", "anode|sep", "sep|cathode", "cathode|al"]

    fig, ax = plt.subplots(figsize=(8.0, 4.0))
    n_iface = len(pp_keys)
    n_design = len(selections)
    x = np.arange(n_iface)
    bar_w = 0.8 / n_design
    cmap = plt.get_cmap("tab10")
    for k, (label, rec) in enumerate(selections):
        ys = [rec[k_] for k_ in pp_keys]
        ax.bar(x + (k - n_design / 2 + 0.5) * bar_w, ys, bar_w,
               color=cmap(k), label=label, edgecolor="black", linewidth=0.4)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(pp_labels)
    ax.set_ylabel("Contact pressure $p_c$ [MPa]")
    ax.set_title("Per-interface contact pressure (negative = tension / debonding risk)",
                 fontsize=11)
    ax.legend(fontsize=8, loc="upper right", ncols=2)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    out = FIG_DIR / "interface_pc_bars.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


# -----------------------------------------------------------------------------
# Figure 6: cycling trajectory
# -----------------------------------------------------------------------------
def fig_cycling():
    if not CYCLE_CSV.exists():
        sys.stderr.write(f"missing {CYCLE_CSV}; skipping cycling figure.\n"
                         "  Re-run with: ../../eel-opt -i pressure_free_stack.i "
                         "cycle_amplitude=1.0 t_cycle=1.0 end_time=3.0 "
                         "Executioner/dt=0.05 Outputs/file_base=rst/cycling_demo\n")
        return
    rows = []
    with open(CYCLE_CSV, newline="") as fh:
        rows = list(csv.DictReader(fh))
    t = np.array([float(r["time"]) for r in rows])
    pi_c = np.array([float(r["Pi_c"]) for r in rows])
    pc_avg = np.array([float(r["pc_avg"]) for r in rows])
    soc = 0.5 * (1 - np.cos(2 * np.pi * np.maximum(t - 1.0, 0.0) / 1.0))
    soc[t < 1.0] = 0.0

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.0, 4.6), sharex=True,
                                   gridspec_kw={"height_ratios": [1, 1.7]},
                                   constrained_layout=True)
    ax1.plot(t, soc, "-", color="purple", label="state of charge $s$")
    ax1.set_ylabel("$s$")
    ax1.axvspan(0, 1, alpha=0.15, color="orange")
    ax1.text(0.5, 0.5, "cure ramp", ha="center", va="center", fontsize=9, color="darkorange")
    ax1.set_ylim(-0.05, 1.05)
    ax1.legend(fontsize=8, loc="upper right")

    ax2.plot(t, pi_c, "-", color="navy", label="$\\Pi_c$")
    ax2.set_xlabel("pseudo-time $t$")
    ax2.set_ylabel("$\\Pi_c$", color="navy")
    ax2.tick_params(axis="y", labelcolor="navy")
    ax2.axvspan(0, 1, alpha=0.15, color="orange")

    ax2b = ax2.twinx()
    ax2b.plot(t, pc_avg, "--", color="tab:red", label="$\\langle p_c \\rangle$")
    ax2b.set_ylabel("$\\langle p_c \\rangle$ [MPa]", color="tab:red")
    ax2b.tick_params(axis="y", labelcolor="tab:red")
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2b.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc="upper left")
    ax2.grid(alpha=0.3)

    ax1.set_title("Lithiation breathing modulates contact pressure during cycling",
                  fontsize=11)
    out = FIG_DIR / "cycling_trajectory.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


# -----------------------------------------------------------------------------
# Figure 7: GPE-modulus lever -- Pi_c vs E_sep at fixed shrink combos
# -----------------------------------------------------------------------------
def fig_E_sep_lever():
    rows, _, _, _ = load_design_map_full()
    if "E_sep_MPa" not in rows[0]:
        print("design_map.csv has no E_sep_MPa column; skipping E_sep lever figure")
        return
    E_seps = sorted({r["E_sep_MPa"] for r in rows})
    eps_pkgs = sorted({r["eps_pkg"] for r in rows})
    eps_GPE_fixed = min({r["eps_GPE"] for r in rows})  # most-favorable corner
    sigma_panels = sorted({r["sigma_stack_MPa"] for r in rows})

    fig, axes = plt.subplots(1, len(sigma_panels), figsize=(4.0 * len(sigma_panels), 3.6),
                             sharey=True, constrained_layout=True)
    if len(sigma_panels) == 1:
        axes = [axes]
    cmap = plt.get_cmap("plasma")
    for ax, sg in zip(axes, sigma_panels):
        for k, ep in enumerate(eps_pkgs):
            ys = []
            for E in E_seps:
                rec = next(r for r in rows
                           if abs(r["E_sep_MPa"] - E) < 1e-9
                           and abs(r["eps_GPE"] - eps_GPE_fixed) < 1e-9
                           and abs(r["eps_pkg"] - ep) < 1e-9
                           and abs(r["sigma_stack_MPa"] - sg) < 1e-9)
                ys.append(rec["Pi_c"])
            ax.plot(E_seps, ys, "-o", color=cmap(k / max(1, len(eps_pkgs) - 1)),
                    label=f"$\\varepsilon_{{\\rm pkg}}={ep:.2f}$", markersize=4)
        ax.axhline(1.0, color="black", lw=0.8, ls="--", alpha=0.6)
        ax.text(E_seps[0] * 1.15, 1.02, "$\\Pi_c = 1$ target",
                fontsize=8, color="black", alpha=0.7)
        ax.set_xscale("log")
        ax.set_xticks(E_seps)
        ax.get_xaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
        ax.set_xlabel("$E_{\\rm sep}$ [MPa]")
        ax.set_title(f"$\\sigma_{{\\rm stack}}={sg:.2f}$ MPa", fontsize=10)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("$\\Pi_c$")
    axes[-1].legend(fontsize=8, loc="lower right")
    fig.suptitle(f"Stiffening the GPE: diminishing returns above $\\sim$25 MPa "
                 f"(at $\\varepsilon_{{\\rm GPE}}={eps_GPE_fixed:.2f}$)",
                 fontsize=10)
    out = FIG_DIR / "Pi_c_vs_E_sep.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


# -----------------------------------------------------------------------------
# Figure 8: proposal-style compact subpanels (sigma_stack = 0 only)
# -----------------------------------------------------------------------------
def fig_proposal_panels():
    # (b) Single-panel heatmap at sigma_stack = 0, E_sep = default.
    rows, eps_GPE, eps_pkg, sigma = load_design_map()
    M = pivot(rows, "Pi_c", eps_GPE, eps_pkg, 0.0)
    vmin = min(r["Pi_c"] for r in rows)
    vmax = max(r["Pi_c"] for r in rows)
    def edges(xs):
        a = np.asarray(xs, dtype=float)
        d = np.diff(a)
        e = np.concatenate(([a[0] - d[0] / 2], a[:-1] + d / 2, [a[-1] + d[-1] / 2]))
        return float(e[0]), float(e[-1])
    x_lo, x_hi = edges(eps_pkg)
    y_lo, y_hi = edges(eps_GPE)
    fig, ax = plt.subplots(figsize=(3.6, 3.0), constrained_layout=True)
    im = ax.imshow(M, origin="lower", aspect="auto",
                   extent=[x_lo, x_hi, y_lo, y_hi],
                   vmin=vmin, vmax=vmax, cmap="viridis")
    try:
        X, Y = np.meshgrid(eps_pkg, eps_GPE)
        cs = ax.contour(X, Y, M, levels=[1.0], colors="white",
                        linewidths=1.4, linestyles="-")
        ax.clabel(cs, inline=True, fontsize=7, fmt="$\\Pi_c=1$")
    except Exception:
        pass
    for i, eg in enumerate(eps_GPE):
        for j, ep in enumerate(eps_pkg):
            color = "white" if M[i, j] < (vmin + vmax) * 0.55 else "black"
            ax.text(ep, eg, f"{M[i,j]:.2f}",
                    ha="center", va="center", fontsize=7, color=color)
    ax.set_xlabel("$\\varepsilon_{\\rm pkg}$")
    ax.set_ylabel("$\\varepsilon_{\\rm GPE}$")
    ax.set_xticks(eps_pkg)
    ax.set_yticks(eps_GPE)
    ax.set_title("$\\Pi_c$ at $\\sigma_{\\rm stack}=0$", fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.85, label="$\\Pi_c$")
    out = FIG_DIR / "prop_heatmap_sigma0.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)

    # (c) Pi_c vs E_sep at sigma_stack = 0.
    rows_full, _, _, _ = load_design_map_full()
    if "E_sep_MPa" not in rows_full[0]:
        return
    E_seps = sorted({r["E_sep_MPa"] for r in rows_full})
    eps_pkgs = sorted({r["eps_pkg"] for r in rows_full})
    eps_GPE_fixed = min({r["eps_GPE"] for r in rows_full})
    fig, ax = plt.subplots(figsize=(3.6, 3.0), constrained_layout=True)
    cmap = plt.get_cmap("plasma")
    for k, ep in enumerate(eps_pkgs):
        ys = [next(r["Pi_c"] for r in rows_full
                   if abs(r["E_sep_MPa"] - E) < 1e-9
                   and abs(r["eps_GPE"] - eps_GPE_fixed) < 1e-9
                   and abs(r["eps_pkg"] - ep) < 1e-9
                   and abs(r["sigma_stack_MPa"]) < 1e-9) for E in E_seps]
        ax.plot(E_seps, ys, "-o", color=cmap(k / max(1, len(eps_pkgs) - 1)),
                label=f"$\\varepsilon_{{\\rm pkg}}={ep:.2f}$", markersize=4)
    ax.axhline(1.0, color="black", lw=0.8, ls="--", alpha=0.6)
    ax.set_xscale("log")
    ax.set_xticks(E_seps)
    ax.get_xaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel("$E_{\\rm sep}$ [MPa]")
    ax.set_ylabel("$\\Pi_c$")
    ax.set_title("Stiffness lever at $\\sigma_{\\rm stack}=0$", fontsize=10)
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(alpha=0.3)
    out = FIG_DIR / "prop_E_sep_sigma0.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


# -----------------------------------------------------------------------------
# Temperature-grid design map (Figures 3-5 analog at multiple T_op)
# -----------------------------------------------------------------------------
def load_T_grid():
    if not T_GRID_CSV.exists():
        sys.stderr.write(f"missing {T_GRID_CSV}; run `python3 sweep.py T` first.\n")
        return None, None, None, None
    rows = []
    with open(T_GRID_CSV, newline="") as fh:
        for r in csv.DictReader(fh):
            rows.append({k: float(v) for k, v in r.items()})
    T_ops = sorted({r["T_op_C"] for r in rows})
    eps_GPE = sorted({r["eps_GPE"] for r in rows})
    eps_pkg = sorted({r["eps_pkg"] for r in rows})
    return rows, eps_GPE, eps_pkg, T_ops


def pivot_T(rows, key, eps_GPE, eps_pkg, T_value):
    M = np.full((len(eps_GPE), len(eps_pkg)), np.nan)
    for r in rows:
        if abs(r["T_op_C"] - T_value) > 1e-9:
            continue
        i = eps_GPE.index(r["eps_GPE"])
        j = eps_pkg.index(r["eps_pkg"])
        M[i, j] = r[key]
    return M


def fig_T_heatmaps():
    """Π_c heatmaps over (eps_GPE, eps_pkg) at sigma_stack = 0, panels = T_op."""
    rows, eps_GPE, eps_pkg, T_ops = load_T_grid()
    if rows is None:
        return
    fig, axes = plt.subplots(1, len(T_ops), figsize=(4.2 * len(T_ops), 3.8),
                             sharey=True, constrained_layout=True)
    if len(T_ops) == 1:
        axes = [axes]
    vmin = min(r["Pi_c"] for r in rows)
    vmax = max(r["Pi_c"] for r in rows)

    def edges(xs):
        a = np.asarray(xs, dtype=float)
        d = np.diff(a)
        e = np.concatenate(([a[0] - d[0] / 2], a[:-1] + d / 2, [a[-1] + d[-1] / 2]))
        return float(e[0]), float(e[-1])
    x_lo, x_hi = edges(eps_pkg)
    y_lo, y_hi = edges(eps_GPE)
    for ax, T in zip(axes, T_ops):
        M = pivot_T(rows, "Pi_c", eps_GPE, eps_pkg, T)
        im = ax.imshow(M, origin="lower", aspect="auto",
                       extent=[x_lo, x_hi, y_lo, y_hi],
                       vmin=vmin, vmax=vmax, cmap="viridis")
        try:
            X, Y = np.meshgrid(eps_pkg, eps_GPE)
            cs = ax.contour(X, Y, M, levels=[1.0], colors="white",
                            linewidths=1.6, linestyles="-")
            ax.clabel(cs, inline=True, fontsize=7, fmt="$\\Pi_c=1$")
        except Exception:
            pass
        for i, eg in enumerate(eps_GPE):
            for j, ep in enumerate(eps_pkg):
                color = "white" if M[i, j] < (vmin + vmax) * 0.55 else "black"
                ax.text(ep, eg, f"{M[i,j]:.2f}",
                        ha="center", va="center", fontsize=7, color=color)
        ax.set_xlabel("$\\varepsilon_{\\rm pkg}$")
        ax.set_title(f"$T_{{\\rm op}}={T:+.0f}^\\circ$C", fontsize=10)
        ax.set_xticks(eps_pkg)
        ax.set_yticks(eps_GPE)
    axes[0].set_ylabel("$\\varepsilon_{\\rm GPE}$")
    fig.colorbar(im, ax=axes, shrink=0.85, label="$\\Pi_c$")
    fig.suptitle("Contact stability $\\Pi_c$ across the DARPA operating window "
                 "(at $\\sigma_{\\rm stack}=0$)", fontsize=11)
    out = FIG_DIR / "T_design_map_heatmap.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def fig_T_lines_eps_pkg():
    """Π_c vs eps_pkg, lines = eps_GPE, panels = T_op."""
    rows, eps_GPE, eps_pkg, T_ops = load_T_grid()
    if rows is None:
        return
    fig, axes = plt.subplots(1, len(T_ops), figsize=(4.0 * len(T_ops), 3.6),
                             sharey=True, constrained_layout=True)
    if len(T_ops) == 1:
        axes = [axes]
    cmap = plt.get_cmap("plasma")
    for ax, T in zip(axes, T_ops):
        for k, eg in enumerate(eps_GPE):
            ys = [next(r["Pi_c"] for r in rows
                       if abs(r["T_op_C"] - T) < 1e-9
                       and abs(r["eps_GPE"] - eg) < 1e-9
                       and abs(r["eps_pkg"] - ep) < 1e-9) for ep in eps_pkg]
            ax.plot(eps_pkg, ys, "-o", color=cmap(k / max(1, len(eps_GPE) - 1)),
                    label=f"$\\varepsilon_{{\\rm GPE}}={eg:.2f}$", markersize=4)
        ax.axhline(1.0, color="black", lw=0.8, ls="--", alpha=0.6)
        ax.set_xlabel("$\\varepsilon_{\\rm pkg}$")
        ax.set_title(f"$T_{{\\rm op}}={T:+.0f}^\\circ$C", fontsize=10)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("$\\Pi_c$")
    axes[-1].legend(fontsize=8, loc="upper left")
    fig.suptitle("$\\Pi_c$ vs.\\ package shrink at three operating temperatures",
                 fontsize=10)
    out = FIG_DIR / "T_Pi_c_vs_eps_pkg.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def fig_T_lines_eps_GPE():
    """Π_c vs eps_GPE, lines = eps_pkg, panels = T_op."""
    rows, eps_GPE, eps_pkg, T_ops = load_T_grid()
    if rows is None:
        return
    fig, axes = plt.subplots(1, len(T_ops), figsize=(4.0 * len(T_ops), 3.6),
                             sharey=True, constrained_layout=True)
    if len(T_ops) == 1:
        axes = [axes]
    cmap = plt.get_cmap("viridis")
    for ax, T in zip(axes, T_ops):
        for k, ep in enumerate(eps_pkg):
            ys = [next(r["Pi_c"] for r in rows
                       if abs(r["T_op_C"] - T) < 1e-9
                       and abs(r["eps_GPE"] - eg) < 1e-9
                       and abs(r["eps_pkg"] - ep) < 1e-9) for eg in eps_GPE]
            ax.plot(eps_GPE, ys, "-s", color=cmap(k / max(1, len(eps_pkg) - 1)),
                    label=f"$\\varepsilon_{{\\rm pkg}}={ep:.2f}$", markersize=4)
        ax.axhline(1.0, color="black", lw=0.8, ls="--", alpha=0.6)
        ax.set_xlabel("$\\varepsilon_{\\rm GPE}$")
        ax.set_title(f"$T_{{\\rm op}}={T:+.0f}^\\circ$C", fontsize=10)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("$\\Pi_c$")
    axes[-1].legend(fontsize=8, loc="upper right")
    fig.suptitle("$\\Pi_c$ vs.\\ GPE shrink at three operating temperatures",
                 fontsize=10)
    out = FIG_DIR / "T_Pi_c_vs_eps_GPE.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


def main() -> int:
    fig_geometry()
    fig_heatmaps()
    fig_lines_eps_pkg()
    fig_lines_eps_GPE()
    fig_interface_bars()
    fig_cycling()
    fig_E_sep_lever()
    fig_proposal_panels()
    fig_T_heatmaps()
    fig_T_lines_eps_pkg()
    fig_T_lines_eps_GPE()
    return 0


if __name__ == "__main__":
    sys.exit(main())
