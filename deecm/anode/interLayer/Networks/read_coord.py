import numpy as np
import json
import itertools
import math
import sys

args = sys.argv
input = args[1]
#input = 'size1640000'
output = args[2]

# Load txt for coordinates and radius of C and Ag particles.
file = open(input, 'r')

lines = file.readlines()
Ag_lst = []
C_lst = []
Top_lst= []

bool = True
SYS_thk = 0.
for line in lines:
    line = line.rstrip()
    particle_data = [float(x) for x in line.split()]
    type_index = particle_data[1]
    coord = particle_data[0:]
    if type_index == 2:
        C_lst.append(coord)
    elif type_index == 1:
        Ag_lst.append(coord)
    elif type_index == 3:
        Top_lst.append(coord)
    elif bool and type_index == 3:
        SYS_thk = particle_data[4] - particle_data[5]
        bool = False

SYS_thk = particle_data[4] - particle_data[5]

num_Ag = len(Ag_lst)
num_C = len(C_lst)
num_Top = len(Top_lst)

Ag_tmp1  = np.asarray(Ag_lst)
Ag_tmp2  = Ag_tmp1[Ag_tmp1[:,0].argsort()]
Ag_array = Ag_tmp2[:,2:]
C_tmp1  = np.asarray(C_lst)
C_tmp2  = C_tmp1[C_tmp1[:,0].argsort()]
C_array = C_tmp2[:,2:]
Top_tmp1  = np.asarray(Top_lst)
Top_tmp2  = Top_tmp1[Top_tmp1[:,0].argsort()]
Top_array = Top_tmp2[:,2:]

dlta = 0.0
combined_array = np.concatenate((Ag_array, C_array, Top_array), axis=0)
num_particle = combined_array.shape[0]
max_x = max(combined_array[:, 0])
max_y = max(combined_array[:, 1])
max_z = max(combined_array[:, 2])
max_coord = max(max_x, max_y, max_z)
max_radius = max(combined_array[:, 9])
grid_size = 2.01 * max_radius
grid_number = int(math.ceil(max_coord / grid_size))
max_index = pow(grid_number, 3) - 1

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
    if i < num_Ag:
        return 'Ag'
    elif i == num_particle:
        return 'Target'
    elif i >= num_Ag and i < num_Ag+num_C:
        return 'C'
    return 'Top'

# Build adj map for each NMC and LPS particles.
# E.g., {0: {1: [5.2, -0.1, 'NMC'], 2: [3.3, -0.2, 'LPS']},
#       1: {0: [5.2, -0.1, 'NMC']}}
# where 0, 1, 2 are labels of particles; 5.2 and 3.3 are inter-particle
#  distance for particles 0 and 1, and particles 0 and 2, respectively;
# -0.1 and -0.2 are overlap of particles 0 and 1, and particles 0 and 2,
# respectively;


# Assign particles to grid boxes.
index_map = {i: [] for i in range(max_index + 1)}
for i in range(num_particle):
    coord_i = combined_array[i, :3]
    index_i = grid_index(grid_index_lst(coord_i))
    index_map[index_i].append(i)

count = sum(len(e) for e in index_map.values())
# Build adjacent lists.
adj = {i: {} for i in range(num_particle + 1)}  # dummy node as target.
vertical_dist = dict()
particle_rads = dict()
particle_type = dict()

vertical_dist[num_particle] = SYS_thk  # For dummy node.
particle_type[num_particle] = 'Target'
particle_rads[num_particle] = 0.

mxGap = 0.

#for i in range(20):
for i in range(num_particle):
    particle_type[i] = type_of_particle(i)
    radius_i = combined_array[i, 9]
    coord_i = combined_array[i, :3]
    z_coord = combined_array[i, 2]
    vertical_dist[i] = z_coord
    particle_rads[i] = radius_i
    if z_coord <= radius_i+dlta:
        adj[i][num_particle] = [z_coord, z_coord-radius_i, 'Target']

    grid_index_i = grid_index_lst(coord_i)
    neighbor_box_i = neighbor_grid_index(grid_index_i)
    for j in neighbor_box_i:
        for k in index_map[j]:
            if k == i:
                continue
            radius_k = combined_array[k, 9]
            coord_k = combined_array[k, :3]
            dist = distance(coord_i, coord_k)
            overlap = dist - (radius_i + radius_k)
            if overlap <= dlta:
                adj[i][k] = [dist, overlap, type_of_particle(k)]
                if overlap < mxGap:
                   mxGap = overlap

print('Maximal Gap: ', mxGap)

inb = 0.
for i in range(num_Top):
    l=i+num_Ag+num_C
    nbi = len(adj[l])
    if nbi < 1:
         inb = inb + 1
#         print("Particle", i,combined_array[i,:3], "is not well contact, with neighbor list:", adj[i])
print('Number of not well contact Top particles ', inb)

data = {'num_Top': num_Top, 'num_Ag': num_Ag,
        'num_C': num_C,
        'num_particle': num_particle + 1,
        'vertical_dist': vertical_dist,
        'particle_rads': particle_rads,
        'type': particle_type,
        'adj': adj,
        'thickness': SYS_thk
}

with open(output, 'w') as fp:
    json.dump(data, fp)

