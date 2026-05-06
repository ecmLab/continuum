# Cubit python file for create pore with curve 
# Curve is sinusoidal with a spline surface 
# pore dimensions are 2.0 um length and 0.5 um average width
import sys 
sys.path.append("/opt/Coreform-Cubit-2020.2/bin")
import cubit
import math

cubit.init(['cubit'])
#####################################################################################################
cubit.cmd('reset')
import math
import numpy as np
# Model parameters 
ly_pore  = 2.5      # height of pore 
lx_pore = 0.5       # Avg width of pore
w_block = 0.1       # width of sinusoidal block 
l_block = 1.5       # height of sinusoidal block 
n_block = 50        # number of points in block 
y_block_start = 0.11 # start of curve in pore wall
x_block_start = 0.5 # start of nominal pore wall
ly_metal = 2.5
lx_metal = 0.5
# ------- Pore Block -----------------
yprev = 0
# Splits of the interval for pores
ysplit = [0, y_block_start, y_block_start + l_block, ly_pore]
npart = len(ysplit) -1
yy = [y_block_start + l_block*i/n_block for i in range(0,n_block+1,1)]
yy2 = [y for y in ysplit if y < y_block_start or y > y_block_start + l_block]
yy = yy + yy2
yy.sort()
x1 = 0
y1 = 0
xprev = x_block_start
## --- Split at dx if no other partitiions exist 

# Solid Electrolyte 
## --- Split at dx if no other partitiions exist 

v_pore_left_curve = []
v_pore_left_curve_list = []
v_pore_left_curve_ends_list = []
v_pore_left_curve_ends = []

for i in range(npart):
    v_pore_left_curve.clear()
    v_pore_left_curve_ends.clear()
    if i > 0:
        vv = cubit.vertex(cubit.get_entities('vertex')[-1])
        print('Vertex {} at {} added to list'.format(vv.id(), vv.coordinates()))
        v_pore_left_curve.append(vv)
        newyy = [ax for ax in yy if ax > ysplit[i] and ax <= ysplit[i+1]]
    else:
        newyy = [ax for ax in yy if ax >=ysplit[i] and ax <= ysplit[i+1]]
    for y in newyy:
        if y >= ysplit[1] and y<=ysplit[2]:
            x = x_block_start - w_block*(math.sin(math.pi*(y-y_block_start)/l_block/2))
            xprev = x 
        else: 
            x = xprev
        v_pore_left_curve.append(cubit.create_vertex(x,y,0))

    if len(v_pore_left_curve) > 1:
        v_pore_left_curve_ends.append(v_pore_left_curve[0])
        v_pore_left_curve_ends.append(v_pore_left_curve[-1])
    else:
        v_pore_left_curve_ends.append(v_pore_left_curve_ends[-1][-1])
        v_pore_left_curve_ends.append(v_pore_left_curve[-1])
        
    v_pore_left_curve_list.append(v_pore_left_curve.copy())    
    # v_pore_left_curve_ends_list.append(v_pore_left_curve_ends.copy())

# Spline for bottom surface of pore
# loop through all vertex list and create curves 
c_pore_left_list = []
c_pore_curve_left = None
for v in v_pore_left_curve_list:
    if len(v) > 2:
        str1 = "create curve spline location vertex " + str(v[0].id()) + " to " + str(v[-1].id())
        print(str1)
        cubit.cmd(str1)
        c_pore_curve_bot = cubit.curve(cubit.get_entities('curve')[-1])
    else: 
        c_pore_curve_bot = cubit.create_curve(v[0], v[-1])
    print('Creating line from {} at {} to {} at {}'.format(v[0].id(), v[0].coordinates(), v[-1].id(), v[-1].coordinates()))
    c_pore_left_list.append(c_pore_curve_bot)
# --- Recreate vertex and curve lists 
c_pore_left_list.clear()
v_pore_left_curve_list.clear()
v_pore_left_curve.clear()
curves = cubit.get_entities('curve')
for c in curves:
    c_pore_left_list.append(cubit.curve(c))
    v_pore_left_curve_list.append(cubit.curve(c).vertices())
pore_curves = c_pore_left_list.copy()

# Rest of Pore wall
v_pore_top_right = cubit.create_vertex(x_block_start + lx_pore, ly_pore, 0)
v_pore_bot_right = cubit.create_vertex(x_block_start + lx_pore, 0, 0)
pore_curves.append(cubit.create_curve(v_pore_left_curve_list[0][0], v_pore_bot_right))
pore_curves.append(cubit.create_curve(v_pore_bot_right, v_pore_top_right))
pore_curves.append(cubit.create_curve(v_pore_top_right, v_pore_left_curve_list[-1][-1]))
pore_surface = cubit.create_surface(pore_curves)

cubit.cmd('create surface rectangle width 1 height 5 zplane ')
ceramic_surface_id = cubit.get_entities('surface')[-1]
ceramic_surface = cubit.surface(ceramic_surface_id)
cubit.cmd('move Surface 2  x 0.5 y -2.5 include_merged ')
cubit.cmd('create surface rectangle width 0.5 height 0.1 zplane ')
cubit.cmd('move Surface 3  x 0.25 y 0.05 include_merged ')


# ------- metal Block -----------------
ly_metal = 2.5
lx_metal = 0.5

yprev = 0
# Splits of the interval for metals
ysplit = [0.1, y_block_start, y_block_start + l_block, ly_metal]
npart = len(ysplit) -1
yy = [y_block_start + l_block*i/n_block for i in range(0,n_block+1,1)]
yy2 = [y for y in ysplit if y < y_block_start or y > y_block_start + l_block]
yy = yy + yy2
yy.sort()
x1 = 0
y1 = 0
xprev = x_block_start
## --- Split at dx if no other partitiions exist 

# Solid Electrolyte 
## --- Split at dx if no other partitiions exist 

v_metal_right_curve = []
v_metal_right_curve_list = []
v_metal_right_curve_ends_list = []
v_metal_right_curve_ends = []

for i in range(npart):
    v_metal_right_curve.clear()
    v_metal_right_curve_ends.clear()
    if i > 0:
        vv = cubit.vertex(cubit.get_entities('vertex')[-1])
        print('Vertex {} at {} added to list'.format(vv.id(), vv.coordinates()))
        v_metal_right_curve.append(vv)
        newyy = [ax for ax in yy if ax > ysplit[i] and ax <= ysplit[i+1]]
    else:
        newyy = [ax for ax in yy if ax >=ysplit[i] and ax <= ysplit[i+1]]
    for y in newyy:
        if y >= ysplit[1] and y<=ysplit[2]:
            x = x_block_start - w_block*(math.sin(math.pi*(y-y_block_start)/l_block/2))
            xprev = x 
        else: 
            x = xprev
        v_metal_right_curve.append(cubit.create_vertex(x,y,0))

    if len(v_metal_right_curve) > 1:
        v_metal_right_curve_ends.append(v_metal_right_curve[0])
        v_metal_right_curve_ends.append(v_metal_right_curve[-1])
    else:
        v_metal_right_curve_ends.append(v_metal_right_curve_ends[-1][-1])
        v_metal_right_curve_ends.append(v_metal_right_curve[-1])
        
    v_metal_right_curve_list.append(v_metal_right_curve.copy())    
    # v_metal_right_curve_ends_list.append(v_metal_right_curve_ends.copy())

# Spline for bottom surface of metal
# loop through all vertex list and create curves 
c_metal_right_list = []
c_metal_curve_right = None
for v in v_metal_right_curve_list:
    if len(v) > 2:
        str1 = "create curve spline location vertex " + str(v[0].id()) + " to " + str(v[-1].id())
        print(str1)
        cubit.cmd(str1)
        c_metal_curve_bot = cubit.curve(cubit.get_entities('curve')[-1])
    else: 
        c_metal_curve_bot = cubit.create_curve(v[0], v[-1])
    print('Creating line from {} at {} to {} at {}'.format(v[0].id(), v[0].coordinates(), v[-1].id(), v[-1].coordinates()))
    c_metal_right_list.append(c_metal_curve_bot)
# # --- Recreate vertex and curve lists 
# c_metal_right_list.clear()
# v_metal_right_curve_list.clear()
# v_metal_right_curve.clear()
# curves = cubit.get_entities('curve')
# for c in curves:
#     c_metal_right_list.append(cubit.curve(c))
#     v_metal_right_curve_list.append(cubit.curve(c).vertices())
metal_curves = c_metal_right_list.copy()

# Rest of metal wall
v_metal_top_left = cubit.create_vertex(x_block_start - lx_metal, ly_metal, 0)
v_metal_bot_left = cubit.create_vertex(x_block_start - lx_metal, 0.1, 0)
metal_curves.append(cubit.create_curve(v_metal_right_curve_list[0][0], v_metal_bot_left))
metal_curves.append(cubit.create_curve(v_metal_bot_left, v_metal_top_left))
metal_curves.append(cubit.create_curve(v_metal_top_left, v_metal_right_curve_list[-1][-1]))
metal_surface = cubit.create_surface(metal_curves)


