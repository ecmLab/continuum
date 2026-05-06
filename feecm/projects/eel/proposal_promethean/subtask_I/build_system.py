"""Build a PTFE-on-graphite system for a LAMMPS cleavage simulation.

Outputs:
  system.data       — LAMMPS data file (atom_style full)
  system_info.json  — box dimensions, atom counts, area for post-processing
"""
import json
import math
from pathlib import Path

import numpy as np

# ---------------------- system parameters ----------------------
# Graphite slab (frozen substrate)
A_GR    = 2.46          # graphite lattice constant (Å)
D_CC_GR = 1.42          # in-layer C-C distance
D_LAYER = 3.35          # interlayer spacing
N_X     = 8             # orthorhombic unit cells along x
N_Y     = 5             # along y
N_LAYER = 3

# PTFE chains: CF3-(CF2)_(N_C-2)-CF3
N_C     = 20            # carbons per chain
BOND_CC = 1.529
BOND_CF = 1.332
ANG_CCC = math.radians(112.7)
ANG_FCF = math.radians(109.5)

N_CHAIN_X = 1           # only one chain per row (chain spans most of Lx)
N_CHAIN_Y = 4           # 4 parallel chains along y
GAP_GR_PTFE = 3.5       # initial vertical gap (Å) above top of graphite

# ---------------------- atom / topology types ----------------------
T_CGR = 1   # graphite C
T_CPF = 2   # PTFE C
T_F   = 3   # F

BT_CC = 1   # bond C-C
BT_CF = 2   # bond C-F

AT_CCC = 1  # angle C-C-C
AT_CCF = 2  # angle C-C-F
AT_FCF = 3  # angle F-C-F

DT_CCCC = 1  # dihedral C-C-C-C (backbone only, F-involved skipped)

MASS = {T_CGR: 12.011, T_CPF: 12.011, T_F: 18.998}


# ---------------------- builders ----------------------
def build_graphite(n_x, n_y, n_layers, z_top):
    """Return (atoms, Lx, Ly).  atoms is a list of (type, x, y, z)."""
    a = A_GR
    by = a * math.sqrt(3.0)
    basis = [(0.0, 0.0),
             (0.0, D_CC_GR),
             (a / 2.0, by / 2.0),
             (a / 2.0, by / 2.0 + D_CC_GR)]
    atoms = []
    for k in range(n_layers):
        z = z_top - (n_layers - 1 - k) * D_LAYER
        for i in range(n_x):
            for j in range(n_y):
                for (bx, byy) in basis:
                    atoms.append((T_CGR,
                                  bx + i * a,
                                  byy + j * by,
                                  z))
    return atoms, n_x * a, n_y * by


def build_ptfe_chain(x_start, y_offset, z_start):
    """Build a planar zig-zag PTFE chain along +x.

    Returns (atoms, bonds, angles, dihedrals) where each element is local-indexed.
    Local indices start at 0.  Atom tuples are (type, x, y, z).
    """
    # zig-zag spacing: backbone C alternates between two z-values.
    # solve for dx, dz given C-C bond length and C-C-C angle.
    cos_a = math.cos(ANG_CCC)
    # cos(theta) = (dx^2 - dz^2)/(dx^2 + dz^2);  dx^2+dz^2 = BOND_CC^2
    dx = BOND_CC * math.sqrt((1.0 + cos_a) / 2.0)
    dz = BOND_CC * math.sqrt((1.0 - cos_a) / 2.0)

    fy = BOND_CF * math.sin(ANG_FCF / 2.0)
    fz = BOND_CF * math.cos(ANG_FCF / 2.0)

    atoms = []
    backbone = []  # indices of backbone C atoms in local list
    for i in range(N_C):
        x = x_start + i * dx
        y = y_offset
        z = z_start + (i % 2) * dz
        backbone.append(len(atoms))
        atoms.append((T_CPF, x, y, z))

        # F's: for interior CF2, two F's perpendicular; for terminal CF3, three F's tetrahedral.
        if i == 0 or i == N_C - 1:
            # terminal: 3 F's around -vec_in (axis pointing away from chain)
            if i == 0:
                vec_in = np.array([dx, 0.0, dz - 0.0])  # to next C  (z next is dz)
            else:
                # last C: previous C is at i-1.  z of last alternates: if N_C even last is z_start+dz; else z_start.
                z_prev = z_start + ((i - 1) % 2) * dz
                vec_in = np.array([-dx, 0.0, z_prev - z])
            vec_in /= np.linalg.norm(vec_in)
            axis = -vec_in
            # build basis perpendicular to axis
            tmp = np.array([0.0, 1.0, 0.0]) if abs(axis[1]) < 0.9 else np.array([1.0, 0.0, 0.0])
            e1 = np.cross(axis, tmp); e1 /= np.linalg.norm(e1)
            e2 = np.cross(axis, e1)
            cos_b = math.cos(math.pi - ANG_FCF)  # angle between -axis and F bond ~= 109.5° from C-C bond
            # equivalently: F bond at ANG_FCF from -axis (the inward direction)... let's just use 109.5° from vec_in.
            # F bond vector = cos(theta)*axis + sin(theta)*(perpendicular)
            # We want angle between F bond and vec_in == ANG_FCF.
            # cos(ANG_FCF) = dot(F_dir, vec_in) = -cos(theta_along_axis)
            theta = math.pi - ANG_FCF  # angle from axis
            for k in range(3):
                phi = k * 2.0 * math.pi / 3.0
                F_dir = (math.cos(theta) * axis
                         + math.sin(theta) * (math.cos(phi) * e1 + math.sin(phi) * e2))
                fpos = np.array([x, y, z]) + BOND_CF * F_dir
                atoms.append((T_F, fpos[0], fpos[1], fpos[2]))
        else:
            # interior CF2: 2 F's perpendicular to backbone, opposite side from C-C-C bisector
            # bisector of (vec_to_prev, vec_to_next) points toward neighbors; F's on opposite side
            # for low-z C (i even): neighbors are at z+dz; F's at z-fz, ±fy in y
            # for high-z C (i odd): F's at z+fz
            sign = -1.0 if (i % 2 == 0) else +1.0
            atoms.append((T_F, x, y + fy, z + sign * fz))
            atoms.append((T_F, x, y - fy, z + sign * fz))

    # ---------- topology ----------
    n_atoms = len(atoms)
    backbone_set = set(backbone)
    f_indices_per_C = {bi: [] for bi in backbone}
    for ai in range(n_atoms):
        if ai in backbone_set:
            continue
        # find nearest backbone C
        ax, ay, az = atoms[ai][1:]
        dists = [((atoms[b][1] - ax) ** 2 + (atoms[b][2] - ay) ** 2 + (atoms[b][3] - az) ** 2, b)
                 for b in backbone]
        dists.sort()
        f_indices_per_C[dists[0][1]].append(ai)

    bonds = []
    for k in range(N_C - 1):
        bonds.append((BT_CC, backbone[k], backbone[k + 1]))
    for bi, fs in f_indices_per_C.items():
        for f in fs:
            bonds.append((BT_CF, bi, f))

    angles = []
    # C-C-C
    for k in range(1, N_C - 1):
        angles.append((AT_CCC, backbone[k - 1], backbone[k], backbone[k + 1]))
    # C-C-F and F-C-F
    for k, bi in enumerate(backbone):
        c_neighbors = []
        if k > 0:
            c_neighbors.append(backbone[k - 1])
        if k < N_C - 1:
            c_neighbors.append(backbone[k + 1])
        fs = f_indices_per_C[bi]
        for cn in c_neighbors:
            for f in fs:
                angles.append((AT_CCF, cn, bi, f))
        for ii in range(len(fs)):
            for jj in range(ii + 1, len(fs)):
                angles.append((AT_FCF, fs[ii], bi, fs[jj]))

    dihedrals = []
    # backbone C-C-C-C only
    for k in range(N_C - 3):
        dihedrals.append((DT_CCCC,
                          backbone[k], backbone[k + 1],
                          backbone[k + 2], backbone[k + 3]))

    return atoms, bonds, angles, dihedrals


# ---------------------- assemble & write ----------------------
def main():
    out_dir = Path(__file__).parent
    # graphite slab top at z = 0
    gr_atoms, Lx, Ly = build_graphite(N_X, N_Y, N_LAYER, z_top=0.0)

    # tile chains across the slab
    chain_x_span = (N_C - 1) * BOND_CC * math.sqrt((1.0 + math.cos(ANG_CCC)) / 2.0)
    # x-start centers each chain in the box
    xs = [(Lx - chain_x_span) / 2.0 for _ in range(N_CHAIN_X)] if N_CHAIN_X >= 1 else []
    # y-spacings: evenly spread N_CHAIN_Y chains
    ys = [Ly * (k + 0.5) / N_CHAIN_Y for k in range(N_CHAIN_Y)]

    chains = []
    for x0 in xs:
        for y0 in ys:
            chains.append(build_ptfe_chain(x_start=x0,
                                           y_offset=y0,
                                           z_start=GAP_GR_PTFE))

    # stack atoms / topology with offsets
    all_atoms = []          # (mol_id, type, x, y, z)
    all_bonds = []
    all_angles = []
    all_dihedrals = []

    # graphite is mol-id 1 (frozen)
    for (t, x, y, z) in gr_atoms:
        all_atoms.append((1, t, x, y, z))

    for ci, (atoms, bonds, angles, dihedrals) in enumerate(chains):
        offset = len(all_atoms)
        mol_id = 2 + ci  # each chain its own mol id
        for (t, x, y, z) in atoms:
            all_atoms.append((mol_id, t, x, y, z))
        for (bt, a1, a2) in bonds:
            all_bonds.append((bt, a1 + offset, a2 + offset))
        for (at, a1, a2, a3) in angles:
            all_angles.append((at, a1 + offset, a2 + offset, a3 + offset))
        for (dt, a1, a2, a3, a4) in dihedrals:
            all_dihedrals.append((dt, a1 + offset, a2 + offset, a3 + offset, a4 + offset))

    # box
    z_lo = -(N_LAYER - 1) * D_LAYER - 5.0
    z_hi = 50.0  # plenty of vacuum above for pulling

    # write LAMMPS data
    n_atoms = len(all_atoms)
    n_bonds = len(all_bonds)
    n_angles = len(all_angles)
    n_dihedrals = len(all_dihedrals)

    with open(out_dir / "system.data", "w") as f:
        f.write("LAMMPS data file: PTFE on graphite cleavage demo\n\n")
        f.write(f"{n_atoms} atoms\n")
        f.write(f"{n_bonds} bonds\n")
        f.write(f"{n_angles} angles\n")
        f.write(f"{n_dihedrals} dihedrals\n\n")
        f.write("3 atom types\n")
        f.write("2 bond types\n")
        f.write("3 angle types\n")
        f.write("1 dihedral types\n\n")
        f.write(f"0.0 {Lx:.6f} xlo xhi\n")
        f.write(f"0.0 {Ly:.6f} ylo yhi\n")
        f.write(f"{z_lo:.6f} {z_hi:.6f} zlo zhi\n\n")
        f.write("Masses\n\n")
        for t in (T_CGR, T_CPF, T_F):
            f.write(f"{t} {MASS[t]:.4f}\n")
        f.write("\nAtoms # full\n\n")
        for i, (mol, t, x, y, z) in enumerate(all_atoms, 1):
            f.write(f"{i} {mol} {t} 0.0 {x:.6f} {y:.6f} {z:.6f}\n")
        if n_bonds:
            f.write("\nBonds\n\n")
            for i, (bt, a1, a2) in enumerate(all_bonds, 1):
                f.write(f"{i} {bt} {a1+1} {a2+1}\n")
        if n_angles:
            f.write("\nAngles\n\n")
            for i, (at, a1, a2, a3) in enumerate(all_angles, 1):
                f.write(f"{i} {at} {a1+1} {a2+1} {a3+1}\n")
        if n_dihedrals:
            f.write("\nDihedrals\n\n")
            for i, (dt, a1, a2, a3, a4) in enumerate(all_dihedrals, 1):
                f.write(f"{i} {dt} {a1+1} {a2+1} {a3+1} {a4+1}\n")

    info = {
        "Lx_A": Lx,
        "Ly_A": Ly,
        "z_lo": z_lo,
        "z_hi": z_hi,
        "n_atoms": n_atoms,
        "n_graphite": len(gr_atoms),
        "n_ptfe_atoms": n_atoms - len(gr_atoms),
        "n_chains": len(chains),
        "n_carbons_per_chain": N_C,
        "interface_area_A2": Lx * Ly,
    }
    with open(out_dir / "system_info.json", "w") as f:
        json.dump(info, f, indent=2)

    print(f"Wrote system.data: {n_atoms} atoms ({len(gr_atoms)} graphite + "
          f"{n_atoms - len(gr_atoms)} PTFE), Lx={Lx:.2f} Å, Ly={Ly:.2f} Å")
    print(f"Interface area: {Lx*Ly:.1f} Å²")


if __name__ == "__main__":
    main()
