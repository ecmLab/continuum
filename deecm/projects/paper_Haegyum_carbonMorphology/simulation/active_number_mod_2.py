import numpy as np
import json
import sys
from collections import deque
from scipy.spatial import cKDTree

if len(sys.argv) != 3:
    sys.exit("Usage: python active_number_mod.py <input.dump> <output.json>")

input_file  = sys.argv[1]
output_file = sys.argv[2]

# ── Parse LIGGGHTS dump (last timestep) ──────────────────────────────
# Dump may contain multiple timesteps; we use the last one.
# Column names are read from the ITEM: ATOMS header line.

with open(input_file, 'r') as fh:
    all_lines = fh.readlines()

timestep_starts = [i for i, l in enumerate(all_lines) if l.startswith("ITEM: TIMESTEP")]
if not timestep_starts:
    sys.exit("ERROR: No 'ITEM: TIMESTEP' found — is this a LIGGGHTS dump file?")

last_header = timestep_starts[-1]

# Parse column names
atoms_header_line = all_lines[last_header + 8].strip()
if not atoms_header_line.startswith("ITEM: ATOMS"):
    sys.exit("ERROR: Expected 'ITEM: ATOMS ...' at line " + str(last_header + 8))
col_names = atoms_header_line.split()[2:]

required = ['id', 'type', 'x', 'y', 'z', 'radius']
col_idx = {}
for name in required:
    if name not in col_names:
        sys.exit(f"ERROR: Required column '{name}' not found in dump header.")
    col_idx[name] = col_names.index(name)

lines = all_lines[last_header + 9:]

DRX_lst = []   # type 1
C_lst = []    # type 2

SYS_thk = 0.0
found_wall = False

for line in lines:
    line = line.rstrip()
    if not line or line.startswith("ITEM:"):
        break
    vals = [float(x) for x in line.split()]
    tp  = vals[col_idx['type']]
    x   = vals[col_idx['x']]
    y   = vals[col_idx['y']]
    z   = vals[col_idx['z']]
    r   = vals[col_idx['radius']]
    pid = vals[col_idx['id']]
    row = [pid, tp, x, y, z, r]
    if tp == 1:
        DRX_lst.append(row)
    elif tp == 2:
        C_lst.append(row)
    elif tp == 3 and not found_wall:
        SYS_thk = z - r
        found_wall = True

num_DRX = len(DRX_lst)
num_C  = len(C_lst)

# Rows are [id, type, x, y, z, radius]. Sort by id, keep x y z radius.
DRX_arr = np.asarray(DRX_lst)
DRX_arr = DRX_arr[DRX_arr[:, 0].argsort()][:, 2:]

C_arr = np.asarray(C_lst) if num_C > 0 else np.empty((0, 4))
if num_C > 0:
    C_arr = C_arr[C_arr[:, 0].argsort()][:, 2:]

# Combined: [0, num_DRX) = type 1;  [num_DRX, num_DRX+num_C) = type 2
combined = np.concatenate((DRX_arr, C_arr), axis=0) if num_C > 0 else DRX_arr.copy()
num_particle = combined.shape[0]

coords = combined[:, :3]
radii  = combined[:, 3]
r_max  = radii.max()

# ── Build KD-tree & adjacency list ───────────────────────────────────
dlta = 0.0

tree = cKDTree(coords)

adj = [[] for _ in range(num_particle)]
bottom_seeds = []
num_contacts = 0

for i in range(num_particle):
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

# ── BFS from bottom (type 2 network only) ────────────────────────────
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

# ── Count active type-1 particles ────────────────────────────────────
num_active_DRX = 0
vol_DRX        = 0.0
vol_active_DRX = 0.0

active_info = {}
for i in range(num_DRX):
    ri = radii[i]
    vol_DRX += ri**3
    on_bottom = coords[i, 2] <= ri + dlta
    touches_network = any(visited[v] for v in adj[i] if v >= num_DRX)
    if on_bottom or touches_network:
        vol_active_DRX += ri**3
        active_info[num_active_DRX] = {
            'particle_id': i,
            'z': float(coords[i, 2]),
            'radius': float(ri)
        }
        num_active_DRX += 1

frac = vol_active_DRX / vol_DRX if vol_DRX > 0 else 0.0

# stdout: num_DRX  num_active  fraction
print(num_DRX, num_active_DRX, frac)

result = {
    'num_DRX':            num_DRX,
    'num_active_DRX':     num_active_DRX,
    'thickness':         SYS_thk,
    'active_DRX_percent': frac,
    'active_particles':  active_info
}

with open(output_file, 'w') as fp:
    json.dump(result, fp, indent=2)