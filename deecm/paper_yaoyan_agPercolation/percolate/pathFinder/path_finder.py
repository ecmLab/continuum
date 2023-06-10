import json
from queue import PriorityQueue
import sys

args = sys.argv
input = args[1]
output = args[2]

# Load JSON file and variables
with open(input, 'rb') as fp:
    data = json.load(fp)

adj_data = data['adj']
num_particle = data['num_Ag']
vertical_dist_data = data['vertical_dist']
particle_rads_data = data['particle_rads']
particle_type_data = data['type']
SYS_thk = data['thickness']

adj = {}
for key in adj_data.keys():
    bool = False
    value = {}
    for i in adj_data[key].keys():
        value[int(i)] = adj_data[key][i]
    adj[int(key)] = value

vertical_dist = {int(key): value for key, value in vertical_dist_data.items()}
particle_rads = {int(key): value for key, value in particle_rads_data.items()}
particle_type = {int(key): value for key, value in particle_type_data.items()}

# Helper function
def heuristic(x):
    if x < 0 or x >= num_particle:
        raise IndexError('Index out of bound.')
    return vertical_dist[x]

def path_finder(s):
    if particle_type[s] != 'Ag':
        raise IndexError('Only consider Ag particles!')

    dist_so_far = [float('inf')] * num_particle
    parent = [None] * num_particle

    pq = PriorityQueue()

    dist_so_far[s] = 0
    source = (heuristic(s), s)
    parent[s] = s

    pq.put(source)
    while not pq.empty():
        _, temp = pq.get()
        if particle_type[temp] == 'Target':
            break
        for v in adj[temp].keys():
            if v == parent[temp]:
                continue
            dist = dist_so_far[temp] + adj[temp][v][0]
            if dist < dist_so_far[v]:
                parent[v] = temp
                dist_so_far[v] = dist
                pq.put((dist + heuristic(v), v))

    if not parent[-1]:
        return [], float('inf')

    shortest_path = [num_particle - 1]
    path_length = 0
    particle_on_path = num_particle - 1
    while True:
        v = parent[particle_on_path]
        if v == particle_on_path:
            break
        shortest_path.append(v)
        path_length += adj[v][particle_on_path][0]
        particle_on_path = v
    return list(reversed(shortest_path)), path_length
###


shortest_paths = {}
num_active_Ag = 0
vol_Ag = 0.
vol_active_Ag = 0.
for s in range(num_particle-1):
    path, length = path_finder(s)
    vol_Ag += particle_rads[s]**3
    if path:
        vol_active_Ag += particle_rads[s]**3
        shortest_paths[num_active_Ag] = [length, vertical_dist[s],  particle_rads[s], path]
        num_active_Ag += 1

data = {'num_Ag': num_particle, 'num_active_Ag': num_active_Ag, 'thickness': SYS_thk, 'active_Ag_percent': vol_active_Ag/vol_Ag, 'path': shortest_paths}
with open(output, 'w') as fp:
    json.dump(data, fp)
