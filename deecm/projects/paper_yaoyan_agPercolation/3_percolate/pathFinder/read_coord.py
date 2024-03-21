import numpy as np
import json
import itertools
import math
import sys

args = sys.argv
input = args[1]
output = args[2]

# Load txt for coordinates and radius of Ag particles.
file = open(input, 'r')
lines = file.readlines()
Ag_lst = []
for line in lines:
    line = line.rstrip()
    particle_data = [float(x) for x in line.split()]
    coord = particle_data[0:]
    Ag_lst.append(coord)

Ag_tmp1  = np.asarray(Ag_lst)
Ag_tmp2  = Ag_tmp1[Ag_tmp1[:,0].argsort()]
Ag_array = Ag_tmp2[:,2:]

### Verify size distribution. ###
# Ag_radius = np.asarray([e[3] for e in Ag_lst])
# print('Ag radius mean: ', np.mean(Ag_radius))
# print('Ag radius std: ', np.std(Ag_radius))

dlta = 0.4
num_particle = Ag_array.shape[0]
max_x = max(Ag_array[:, 0])
max_y = max(Ag_array[:, 1])
max_z = max(Ag_array[:, 2])
max_coord = max(max_x, max_y, max_z)
max_radius = max(Ag_array[:, 3])
grid_size = 2.01 * max_radius
grid_number = int(math.ceil(max_coord / grid_size))
max_index = pow(grid_number, 3) - 1
SYS_thk = max_z + max_radius

# Helper function
def grid_index_lst(coord):
    return [int(e / grid_size) for e in coord]

def grid_index(lst):
    """
    Converts a list of index to a 1D index
    :param lst: 3D grid_box index.
    :return: grid_box index
    """
    return lst[0] + grid_number * lst[1] + grid_number * grid_number * lst[2]


def neighbor_grid_index(lst):
    neighbor_list = []
    x_index = [lst[0] - 1, lst[0], lst[0] + 1]
    y_index = [lst[1] - 1, lst[1], lst[1] + 1]
    z_index = [lst[2] - 1, lst[2], lst[2] + 1]
    x_index = [e for e in x_index if 0 <= e < grid_number]
    y_index = [e for e in y_index if 0 <= e < grid_number]
    z_index = [e for e in z_index if 0 <= e < grid_number]
    index_lst = [x_index, y_index, z_index]
    for comb in itertools.product(*index_lst):
        neighbor_list.append(grid_index(comb))
    return neighbor_list


def distance(a, b):
    return np.linalg.norm(a-b)


def type_of_particle(i):
    if i == num_particle:
        return 'Target'
    return 'Ag'

# Build adj map for each Ag particle.
# E.g., {0: {1: [5.2, -0.1, 'Ag'], 1: {0: [5.2, -0.1, 'Target']}}
# where 0, 1, 2 are labels of particles; 5.2 and 3.3 are inter-particle
#  distance for particles 0 and 1, and particles 0 and 2, respectively;
# -0.1 and -0.2 are overlap of particles 0 and 1, and particles 0 and 2,
# respectively;

# Assign particles to grid boxes.
index_map = {i: [] for i in range(max_index + 1)}
for i in range(num_particle):
    coord_i = Ag_array[i, :3]
    index_i = grid_index(grid_index_lst(coord_i))
    index_map[index_i].append(i)

count = sum(len(e) for e in index_map.values())
# Build adjacent lists.
adj = {i: {} for i in range(num_particle + 1)}  # dummy node as target.
vertical_dist = dict()
particle_rads = dict()
particle_type = dict()

vertical_dist[num_particle] = 0  # For dummy node.
particle_type[num_particle] = 'Target'
particle_rads[num_particle] = 0.

mxGap = 0.
for i in range(num_particle):
    particle_type[i] = type_of_particle(i)
    radius_i = Ag_array[i, 3]
    coord_i = Ag_array[i, :3]
    z_coord = Ag_array[i, 2]
    vertical_dist[i] = z_coord
    particle_rads[i] = radius_i
    if z_coord <= radius_i+dlta:
        adj[i][num_particle] = [z_coord, 0, 'Target']

    grid_index_i = grid_index_lst(coord_i)
    neighbor_box_i = neighbor_grid_index(grid_index_i)
    for j in neighbor_box_i:
        for k in index_map[j]:
            if k == i:
                continue
            radius_k = Ag_array[k, 3]
            coord_k = Ag_array[k, :3]
            dist = distance(coord_i, coord_k)
            overlap = dist - (radius_i + radius_k)
            if overlap <= dlta:
                adj[i][k] = [dist, overlap, type_of_particle(k)]
                if overlap < mxGap:
                   mxGap = overlap

print('Maximal Gap: ', mxGap)

data = {'num_Ag': num_particle + 1, 'vertical_dist': vertical_dist,
        'particle_rads': particle_rads, 'type': particle_type,
        'adj': adj,  'thickness': SYS_thk
}
with open(output, 'w') as fp:
    json.dump(data, fp)
