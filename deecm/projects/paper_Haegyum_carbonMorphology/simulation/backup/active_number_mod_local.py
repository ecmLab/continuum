import numpy as np
import json
import sys
from collections import deque
from scipy.spatial import cKDTree

# ── Default paths (override via command line) ─────────────────────────
DEFAULT_INPUT  = "post/cnt_LD10_mr1_pck500.dump"
DEFAULT_OUTPUT = "post/active_cnt_LD10_mr1_pck500.json"

input_file  = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_INPUT
output_file = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_OUTPUT

print(f"Input:  {input_file}")
print(f"Output: {output_file}")

# ── Parse LIGGGHTS dump ──────────────────────────────────────────────
# Dump files can contain multiple timesteps, each with a 9-line header.
# We read only the LAST timestep in the file.
# Dump columns: id type x y z radius fx fy fz

with open(input_file, 'r') as fh:
    all_lines = fh.readlines()

timestep_starts = [i for i, l in enumerate(all_lines) if l.startswith("ITEM: TIMESTEP")]
if not timestep_starts:
    sys.exit("ERROR: No 'ITEM: TIMESTEP' found — is this a LIGGGHTS dump file?")

last_header = timestep_starts[-1]

# Parse column names from "ITEM: ATOMS id type x y z ..." header
atoms_header_line = all_lines[last_header + 8].strip()
if not atoms_header_line.startswith("ITEM: ATOMS"):
    sys.exit("ERROR: Expected 'ITEM: ATOMS ...' at line " + str(last_header + 8))
col_names = atoms_header_line.split()[2:]  # strip "ITEM:" and "ATOMS"
print(f"Dump columns: {' '.join(col_names)}")

# Map required column names to indices
required = ['id', 'type', 'x', 'y', 'z', 'radius']
col_idx = {}
for name in required:
    if name not in col_names:
        sys.exit(f"ERROR: Required column '{name}' not found in dump header.")
    col_idx[name] = col_names.index(name)

lines = all_lines[last_header + 9:]
print(f"Using last timestep (header at line {last_header})")

Ag_lst = []   # type 1
C_lst = []    # type 2

SYS_thk = 0.0
found_wall = False

for line in lines:
    line = line.rstrip()
    if not line or line.startswith("ITEM:"):
        break
    vals = [float(x) for x in line.split()]
    tp = vals[col_idx['type']]
    x  = vals[col_idx['x']]
    y  = vals[col_idx['y']]
    z  = vals[col_idx['z']]
    r  = vals[col_idx['radius']]
    pid = vals[col_idx['id']]
    row = [pid, tp, x, y, z, r]
    if tp == 1:
        Ag_lst.append(row)
    elif tp == 2:
        C_lst.append(row)
    elif tp == 3 and not found_wall:
        SYS_thk = z - r
        found_wall = True

num_Ag = len(Ag_lst)
num_C  = len(C_lst)

print(f"Type 1 (Ag) particles: {num_Ag}")
print(f"Type 2 (C)  particles: {num_C}")
print(f"System thickness:      {SYS_thk:.4f}")

# Rows are [id, type, x, y, z, radius]. Sort by id, keep x y z radius.
Ag_arr = np.asarray(Ag_lst)
Ag_arr = Ag_arr[Ag_arr[:, 0].argsort()][:, 2:]   # columns: x y z radius

C_arr = np.asarray(C_lst) if num_C > 0 else np.empty((0, 4))
if num_C > 0:
    C_arr = C_arr[C_arr[:, 0].argsort()][:, 2:]

# Combined: [0, num_Ag) = type 1;  [num_Ag, num_Ag+num_C) = type 2
combined = np.concatenate((Ag_arr, C_arr), axis=0) if num_C > 0 else Ag_arr.copy()
num_particle = combined.shape[0]

coords = combined[:, :3]   # (N, 3)
radii  = combined[:, 3]    # (N,)

r_min = radii.min()
r_max = radii.max()
r_mean = radii.mean()
print(f"Radius  min: {r_min:.4f}  max: {r_max:.4f}  mean: {r_mean:.4f}")

# ── Build KD-tree for neighbour search ────────────────────────────────
dlta = 0.0  # overlap tolerance

print("Building KD-tree ...")
tree = cKDTree(coords)

# ── Build adjacency list & identify bottom-touching particles ────────
# Per-particle query: search radius = r_i + r_max so small particles
# search small volumes. Process in batches for progress reporting.
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
    z_i = ci[2]

    if z_i <= ri + dlta and i >= num_Ag:  # only type 2 seeds
        bottom_seeds.append(i)

    # Only search as far as this particle can possibly touch another
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

# ── BFS from bottom interface (type 2 only) ──────────────────────────
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
        if not visited[v] and v >= num_Ag:  # only traverse type 2
            visited[v] = True
            queue.append(v)

total_visited = visited.sum()
print(f"BFS reached {total_visited} / {num_particle} particles")

# ── Count active type-1 (Ag) particles ───────────────────────────────
# Type 1 is active if it contacts any type-2 particle in the connected network
num_active_Ag = 0
vol_Ag        = 0.0
vol_active_Ag = 0.0

active_info = {}
for i in range(num_Ag):
    ri = radii[i]
    vol_Ag += ri**3
    # Check if this type-1 particle touches any visited type-2 neighbor
    touches_network = any(visited[v] for v in adj[i] if v >= num_Ag)
    if touches_network:
        vol_active_Ag += ri**3
        active_info[num_active_Ag] = {
            'particle_id': i,
            'z': float(coords[i, 2]),
            'radius': float(ri)
        }
        num_active_Ag += 1

frac = vol_active_Ag / vol_Ag if vol_Ag > 0 else 0.0

print("=" * 50)
print(f"Total type-1 (Ag):    {num_Ag}")
print(f"Active type-1 (Ag):   {num_active_Ag}")
print(f"Active volume frac:   {frac:.6f}")
print("=" * 50)

result = {
    'num_Ag':            num_Ag,
    'num_active_Ag':     num_active_Ag,
    'thickness':         SYS_thk,
    'active_Ag_percent': frac,
    'active_particles':  active_info
}

with open(output_file, 'w') as fp:
    json.dump(result, fp, indent=2)

print(f"Results saved to {output_file}")