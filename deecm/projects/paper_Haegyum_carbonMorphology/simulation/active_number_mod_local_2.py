# =============================================================================
# active_number.py  —  Active-particle analysis for LIGGGHTS DEM electrodes
# =============================================================================
#
# PURPOSE
#   Reads a LIGGGHTS dump file and identifies which type-1 "DRX" particles
#   (e.g. active-material grains) are electrochemically accessible — i.e.
#   connected to the bottom current-collector via a percolating type-2
#   conductive network (e.g. carbon black / binder).
#
# INPUTS
#   <input.dump>   LIGGGHTS dump file (single or multi-timestep; last used).
#                  Required columns: id  type  x  y  z  radius  mass
#                  Particle types:
#                    1 = DRX  (active material)
#                    2 = C    (conductive network)
#                    3 = Wall (one particle; defines system thickness)
#
# OUTPUTS
#   <output.json>          Summary: total / active DRX count, active volume
#                          fraction, system thickness, and per-particle
#                          details (z, radius) for every active DRX particle.
#   <output.data>          LAMMPS-format data file containing ONLY the active
#                          type-1 and type-2 particles (diameter + density
#                          derived from radius and mass), ready for re-import
#                          into a follow-up simulation.
#
# ALGORITHM
#   1. Parse  — extract particle positions, radii, and masses from the last
#               timestep; record simulation box bounds.
#   2. Contact — build a cKDTree over all particle centres; connect any two
#               particles whose centre-to-centre distance ≤ r_i + r_j.
#   3. Seed    — mark every type-2 particle touching the bottom wall (z ≤ r)
#               as a network seed.
#   4. BFS     — breadth-first search propagates through the type-2 contact
#               network from the seeds, flagging all reachable (percolating)
#               conductive particles as "visited".
#   5. Active  — a type-1 (DRX) particle is active if it sits on the bottom
#               wall OR touches at least one visited type-2 particle.
#   6. Export  — write the JSON summary and a filtered LAMMPS .data file
#               containing only the active particles.
#
# USAGE
#   python active_number.py [input.dump] [output.json]
#   (defaults to post/rst_mesh_cnt10_mr6.dump / post/active_cnt10_mr6.json)
# =============================================================================


import numpy as np
import json
import sys
from collections import deque
from scipy.spatial import cKDTree

# ── Default paths (override via command line) ─────────────────────────
DEFAULT_INPUT  = "post/rst_mesh_cnt10_mr6.dump"
DEFAULT_OUTPUT = "post/active_cnt10_mr6.json"

input_file  = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_INPUT
output_file = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_OUTPUT

# Derive .data filename from JSON output
data_output = output_file.replace('.json', '.data')

print(f"Input:  {input_file}")
print(f"Output: {output_file}")
print(f"Data:   {data_output}")

# ── Parse LIGGGHTS dump (last timestep) ──────────────────────────────
with open(input_file, 'r') as fh:
    all_lines = fh.readlines()

timestep_starts = [i for i, l in enumerate(all_lines) if l.startswith("ITEM: TIMESTEP")]
if not timestep_starts:
    sys.exit("ERROR: No 'ITEM: TIMESTEP' found — is this a LIGGGHTS dump file?")

last_header = timestep_starts[-1]

# Parse box bounds (lines +5, +6, +7 relative to ITEM: TIMESTEP)
box_x = [float(v) for v in all_lines[last_header + 5].split()]
box_y = [float(v) for v in all_lines[last_header + 6].split()]
box_z = [float(v) for v in all_lines[last_header + 7].split()]
print(f"Box X: {box_x[0]:.4f} to {box_x[1]:.4f}")
print(f"Box Y: {box_y[0]:.4f} to {box_y[1]:.4f}")
print(f"Box Z: {box_z[0]:.4f} to {box_z[1]:.4f}")

# Parse column names from ITEM: ATOMS header
atoms_header_line = all_lines[last_header + 8].strip()
if not atoms_header_line.startswith("ITEM: ATOMS"):
    sys.exit("ERROR: Expected 'ITEM: ATOMS ...' at line " + str(last_header + 8))
col_names = atoms_header_line.split()[2:]
print(f"Dump columns: {' '.join(col_names)}")

required = ['id', 'type', 'x', 'y', 'z', 'radius', 'mass']
col_idx = {}
for name in required:
    if name not in col_names:
        sys.exit(f"ERROR: Required column '{name}' not found in dump header.")
    col_idx[name] = col_names.index(name)

lines = all_lines[last_header + 9:]
print(f"Using last timestep (header at line {last_header})")

DRX_lst = []   # type 1
C_lst = []    # type 2

SYS_thk = 0.0
found_wall = False

for line in lines:
    line = line.rstrip()
    if not line or line.startswith("ITEM:"):
        break
    vals = [float(x) for x in line.split()]
    tp   = vals[col_idx['type']]
    x    = vals[col_idx['x']]
    y    = vals[col_idx['y']]
    z    = vals[col_idx['z']]
    r    = vals[col_idx['radius']]
    mass = vals[col_idx['mass']]
    pid  = vals[col_idx['id']]
    # row: id, type, x, y, z, radius, mass
    row = [pid, tp, x, y, z, r, mass]
    if tp == 1:
        DRX_lst.append(row)
    elif tp == 2:
        C_lst.append(row)
    elif tp == 3 and not found_wall:
        SYS_thk = z - r
        found_wall = True

num_DRX = len(DRX_lst)
num_C  = len(C_lst)

print(f"Type 1 (DRX) particles: {num_DRX}")
print(f"Type 2 (C)  particles: {num_C}")
print(f"System thickness:      {SYS_thk:.4f}")

# Rows: [id, type, x, y, z, radius, mass]. Sort by id, keep x y z radius mass.
DRX_arr = np.asarray(DRX_lst)
DRX_arr = DRX_arr[DRX_arr[:, 0].argsort()][:, 2:]   # cols: x y z radius mass

C_arr = np.asarray(C_lst) if num_C > 0 else np.empty((0, 5))
if num_C > 0:
    C_arr = C_arr[C_arr[:, 0].argsort()][:, 2:]

# Combined: [0, num_DRX) = type 1;  [num_DRX, num_DRX+num_C) = type 2
combined = np.concatenate((DRX_arr, C_arr), axis=0) if num_C > 0 else DRX_arr.copy()
num_particle = combined.shape[0]

coords = combined[:, :3]   # x y z
radii  = combined[:, 3]    # radius
masses = combined[:, 4]    # mass
r_max  = radii.max()

r_min  = radii.min()
r_mean = radii.mean()
print(f"Radius  min: {r_min:.4f}  max: {r_max:.4f}  mean: {r_mean:.4f}")

# ── Build KD-tree & adjacency list ───────────────────────────────────
dlta = 0.0

print("Building KD-tree ...")
tree = cKDTree(coords)

print("Finding contacts (per-particle KD-tree queries) ...")
adj = [[] for _ in range(num_particle)]
bottom_seeds = []
num_contacts = 0

batch = max(num_particle // 20, 1)
for i in range(num_particle):
    if i % batch == 0:
        print(f"  {i}/{num_particle} ({100*i//num_particle}%)")

    ri  = radii[i]
    ci  = coords[i]

    # Only type 2 particles seed from z = 0
    if ci[2] <= ri + dlta and i >= num_DRX:
        bottom_seeds.append(i)

    search_r = ri + r_max + dlta
    candidates = tree.query_ball_point(ci, search_r)

    for k in candidates:
        if k <= i:
            continue
        rk = radii[k]
        dist = np.linalg.norm(ci - coords[k])
        if dist <= ri + rk + dlta:
            adj[i].append(k)
            adj[k].append(i)
            num_contacts += 1

print(f"Contacts: {num_contacts}  |  Bottom seeds: {len(bottom_seeds)}")

# ── BFS from bottom (type 2 network only) ────────────────────────────
print("Running BFS from z = 0 (type 2 network only) ...")
visited = np.zeros(num_particle, dtype=bool)
queue = deque()

for s in bottom_seeds:
    if not visited[s]:
        visited[s] = True
        queue.append(s)

while queue:
    u = queue.popleft()
    for v in adj[u]:
        if not visited[v] and v >= num_DRX:  # only traverse type 2
            visited[v] = True
            queue.append(v)

total_visited_C = visited[num_DRX:].sum()
print(f"BFS reached {total_visited_C} / {num_C} type-2 particles")

# ── Count active type-1 particles ────────────────────────────────────
num_active_DRX = 0
vol_DRX        = 0.0
vol_active_DRX = 0.0
active_DRX_mask = np.zeros(num_DRX, dtype=bool)

active_info = {}
for i in range(num_DRX):
    ri = radii[i]
    vol_DRX += ri**3
    on_bottom = coords[i, 2] <= ri + dlta
    touches_network = any(visited[v] for v in adj[i] if v >= num_DRX)
    if on_bottom or touches_network:
        vol_active_DRX += ri**3
        active_DRX_mask[i] = True
        active_info[num_active_DRX] = {
            'particle_id': i,
            'z': float(coords[i, 2]),
            'radius': float(ri)
        }
        num_active_DRX += 1

frac = vol_active_DRX / vol_DRX if vol_DRX > 0 else 0.0

print("=" * 50)
print(f"Total type-1 (DRX):    {num_DRX}")
print(f"Active type-1 (DRX):   {num_active_DRX}")
print(f"Active type-2 (C):    {total_visited_C}")
print(f"Active volume frac:   {frac:.6f}")
print("=" * 50)

# ── Save JSON ────────────────────────────────────────────────────────
result = {
    'num_DRX':            num_DRX,
    'num_active_DRX':     num_active_DRX,
    'thickness':         SYS_thk,
    'active_DRX_percent': frac,
    'active_particles':  active_info
}

with open(output_file, 'w') as fp:
    json.dump(result, fp, indent=2)
print(f"JSON saved to {output_file}")

# ── Save LAMMPS .data file (active particles only) ───────────────────
# Collect active type-1 and active type-2 particles
active_rows = []  # each: (type, diameter, density, x, y, z)

for i in range(num_DRX):
    if not active_DRX_mask[i]:
        continue
    ri = radii[i]
    mi = masses[i]
    vol = (4.0 / 3.0) * np.pi * ri**3
    density = mi / vol if vol > 0 else 0.0
    active_rows.append((1, 2.0 * ri, density,
                         coords[i, 0], coords[i, 1], coords[i, 2]))

for i in range(num_DRX, num_particle):
    if not visited[i]:
        continue
    ri = radii[i]
    mi = masses[i]
    vol = (4.0 / 3.0) * np.pi * ri**3
    density = mi / vol if vol > 0 else 0.0
    active_rows.append((2, 2.0 * ri, density,
                         coords[i, 0], coords[i, 1], coords[i, 2]))

num_active_total = len(active_rows)

with open(data_output, 'w') as f:
    f.write("LAMMPS data file\n\n")
    f.write(f"{num_active_total}  atoms\n")
    f.write(f"   2 atom types\n\n")
    f.write(f"  {box_x[0]:.19f}  {box_x[1]:.19f} xlo xhi\n")
    f.write(f"  {box_y[0]:.19f}  {box_y[1]:.19f} ylo yhi\n")
    f.write(f"  {box_z[0]:.19f}  {box_z[1]:.19f} zlo zhi\n\n")
    f.write("Atoms # sphere\n\n")
    for idx, (tp, diam, dens, x, y, z) in enumerate(active_rows, start=1):
        f.write(f"{idx:>8d}  {tp:>3d}  {diam:.5f}  {dens:.5f}   "
                f"{x:.15f}   {y:.15f}   {z:.15f}\n")

print(f"Data file saved to {data_output}  ({num_active_total} active atoms)")