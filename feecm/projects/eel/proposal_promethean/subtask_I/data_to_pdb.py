"""Convert a LAMMPS 'full'-style data file to a PDB file for VMD.

Usage: python data_to_pdb.py [input.data] [output.pdb]
       (defaults: system.data -> system.pdb)
"""
import re
import sys
from pathlib import Path

# atom-type -> (element, atom-name (cols 13-16), residue name)
TYPE_INFO = {
    1: ("C", "Cgr", "GRA"),   # graphite carbon
    2: ("C", "Cpt", "PTF"),   # PTFE carbon
    3: ("F", "F",   "PTF"),   # PTFE fluorine
}

SECTION_KEYWORDS = {
    "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers",
    "Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs",
    "Dihedral Coeffs", "Improper Coeffs",
}


def parse_data(path):
    """Return (sections-dict, box-dict).  sections values are lists of raw lines."""
    text = Path(path).read_text().splitlines()
    sections = {}
    box = {}
    cur = None
    for raw in text:
        line = raw.split("#")[0].rstrip()
        stripped = line.strip()
        if not stripped:
            continue
        m = re.match(r"^([\d\.\-eE\+]+)\s+([\d\.\-eE\+]+)\s+(xlo xhi|ylo yhi|zlo zhi)$",
                     stripped)
        if m:
            box[m.group(3)] = (float(m.group(1)), float(m.group(2)))
            continue
        first = stripped.split()[0]
        if first in SECTION_KEYWORDS:
            cur = first
            sections.setdefault(cur, [])
            continue
        if cur:
            sections[cur].append(stripped)
    return sections, box


def read_atoms(lines):
    """LAMMPS 'full' atom_style: id mol type q x y z [...]."""
    atoms = []
    for l in lines:
        toks = l.split()
        if len(toks) < 7:
            continue
        atoms.append({
            "id":   int(toks[0]),
            "mol":  int(toks[1]),
            "type": int(toks[2]),
            "q":    float(toks[3]),
            "x":    float(toks[4]),
            "y":    float(toks[5]),
            "z":    float(toks[6]),
        })
    atoms.sort(key=lambda a: a["id"])
    return atoms


def read_bonds(lines):
    bonds = []
    for l in lines:
        toks = l.split()
        if len(toks) < 4:
            continue
        bonds.append((int(toks[2]), int(toks[3])))
    return bonds


def pdb_atom_line(serial, name, resname, chain, resseq, x, y, z, elem):
    """Strict PDB column layout for HETATM records."""
    # Atom name padding rule: 1-char element gets a leading space (cols 13-14).
    if len(elem) == 1 and len(name) <= 3:
        name_field = f" {name:<3s}"
    else:
        name_field = f"{name:<4s}"
    return (f"HETATM{serial:>5d} {name_field} {resname:>3s} {chain:1s}{resseq:>4d}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{0.00:>6.2f}"
            f"          {elem:>2s}")


def write_pdb(atoms, bonds, box, out_path):
    # adjacency (1-indexed)
    adj = {}
    for (a, b) in bonds:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)

    # one chain per molecule id (cycle A..Z if needed)
    chains_used = sorted({a["mol"] for a in atoms})
    chain_map = {}
    for i, mol in enumerate(chains_used):
        c = chr(ord("A") + i) if i < 26 else chr(ord("a") + (i - 26)) if i < 52 else " "
        chain_map[mol] = c

    # one residue per molecule
    resseq_map = {mol: (i + 1) % 10000 for i, mol in enumerate(chains_used)}

    out = []
    out.append("REMARK   PTFE-on-graphite system; generated from LAMMPS data file")
    out.append("REMARK   chain A = graphite slab; chains B,C,... = PTFE oligomers")
    if all(k in box for k in ("xlo xhi", "ylo yhi", "zlo zhi")):
        a = box["xlo xhi"][1] - box["xlo xhi"][0]
        b = box["ylo yhi"][1] - box["ylo yhi"][0]
        c = box["zlo zhi"][1] - box["zlo zhi"][0]
        out.append(f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}"
                   f"  90.00  90.00  90.00 P 1           1")

    prev_mol = None
    for atom in atoms:
        if prev_mol is not None and atom["mol"] != prev_mol:
            out.append("TER")
        elem, name, resname = TYPE_INFO.get(atom["type"], ("X", "X", "UNK"))
        out.append(pdb_atom_line(
            serial=atom["id"] % 100000,
            name=name,
            resname=resname,
            chain=chain_map[atom["mol"]],
            resseq=resseq_map[atom["mol"]],
            x=atom["x"], y=atom["y"], z=atom["z"],
            elem=elem,
        ))
        prev_mol = atom["mol"]
    out.append("TER")

    for atom_id in sorted(adj.keys()):
        partners = sorted(adj[atom_id])
        for i in range(0, len(partners), 4):
            chunk = partners[i:i + 4]
            line = f"CONECT{atom_id % 100000:>5d}"
            for p in chunk:
                line += f"{p % 100000:>5d}"
            out.append(line)

    out.append("END")
    Path(out_path).write_text("\n".join(out) + "\n")

    n_chains = len(chains_used)
    print(f"Wrote {out_path}: {len(atoms)} atoms, {len(bonds)} bonds, "
          f"{n_chains} chains/residues.")


def main():
    here = Path(__file__).parent
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else here / "system.data"
    dst = Path(sys.argv[2]) if len(sys.argv) > 2 else src.with_suffix(".pdb")
    sections, box = parse_data(src)
    atoms = read_atoms(sections.get("Atoms", []))
    bonds = read_bonds(sections.get("Bonds", []))
    if not atoms:
        sys.exit(f"No atoms parsed from {src}")
    write_pdb(atoms, bonds, box, dst)


if __name__ == "__main__":
    main()
