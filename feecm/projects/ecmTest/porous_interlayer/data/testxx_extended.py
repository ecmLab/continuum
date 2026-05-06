# Python test for cubit curved surface model
import sys 
sys.path.append("/opt/Coreform-Cubit-2020.2/bin")
import cubit
import math
import numpy as np
cubit.init(['cubit'])
#####################################################################################################

ly   = 50       # height of the model, in unit um
dw   = 0.2       # width of the defect located at the top middle, in unit um
dh   = 1.0       # length of the defect, in unit um
dInt = 0.1      # Initial thickenss of the interlayer, in unit um
nSE  = 101      # Number of discretization points of the cosine shape defect for SE
nSE2 = 50
nLi  = 101      # Number of discretization points of the cosine shape defect for Li
nInt = 101      # Number of discretization points of the cosine shape defect for interlyer
lx   = 50       # width of the model, in unit um
# Interlayer-1 bottom
interlayer_curves = []

yprev = 0
# Splits of the interval for pores
# xsplit = [0, dw/4, dw/4 + 0.1, dw/4+0.2, 3*dw/4+0.2, 3*dw/4+0.3, 5*dw/4+0.3, 5*dw/4+0.4, 7*dw/4+0.4, lx]
# xspace = [[dw, dw+2]]
xsplit = [0, dw, lx/2]
npart = len(xsplit) -1
xx = [dw*i/nSE for i in range(0,nSE+1,1)]
xx2 = [x for x in xsplit if x > dw]
xx = xx + xx2
xx.sort()
x1 = 0
y1 = 0
yprev = 0
## --- Split at dx if no other partitiions exist 

# Solid Electrolyte 
## --- Split at dx if no other partitiions exist 

v_setop_curve = []
v_setop_curve_list = []
v_setop_curve_ends_list = []
v_setop_curve_ends = []

for i in range(npart):
    v_setop_curve.clear()
    v_setop_curve_ends.clear()
    if i > 0:
        vv = cubit.vertex(cubit.get_entities('vertex')[-1])
        print('Vertex {} at {} added to list'.format(vv.id(), vv.coordinates()))
        v_setop_curve.append(vv)
        newxx = [ax for ax in xx if ax > xsplit[i] and ax <= xsplit[i+1]]
    else:
        newxx = [ax for ax in xx if ax >=xsplit[i] and ax <= xsplit[i+1]]
    for x in newxx:
        if x <= dw:
            y = -dh/2.0*(1.0 + math.cos(math.pi*x/dw))
            yprev = y 
        else: 
            y = yprev
        v_setop_curve.append(cubit.create_vertex(x,y,0))

    if len(v_setop_curve) > 1:
        v_setop_curve_ends.append(v_setop_curve[0])
        v_setop_curve_ends.append(v_setop_curve[-1])
    else:
        v_setop_curve_ends.append(v_setop_curve_ends[-1][-1])
        v_setop_curve_ends.append(v_setop_curve[-1])
        
    v_setop_curve_list.append(v_setop_curve.copy())    
    # v_setop_curve_ends_list.append(v_setop_curve_ends.copy())

# Spline for bottom surface of se
# loop through all vertex list and create curves 
c_se_top_list = []
c_se_curve_top = None
for v in v_setop_curve_list:
    if len(v) > 2:
        str1 = "create curve spline location vertex " + str(v[0].id()) + " to " + str(v[-1].id())
        print(str1)
        cubit.cmd(str1)
        c_se_curve_bot = cubit.curve(cubit.get_entities('curve')[-1])
    else: 
        c_se_curve_bot = cubit.create_curve(v[0], v[-1])
    print('Creating line from {} at {} to {} at {}'.format(v[0].id(), v[0].coordinates(), v[-1].id(), v[-1].coordinates()))
    c_se_top_list.append(c_se_curve_bot)
# --- Recreate vertex and curve lists 
c_se_top_list.clear()
v_setop_curve_list.clear()
v_setop_curve.clear()
curves = cubit.get_entities('curve')
for c in curves:
    c_se_top_list.append(cubit.curve(c))
    v_setop_curve_list.append(cubit.curve(c).vertices())
se_curves = c_se_top_list.copy()
# Rest of Solid Electrolyte 
v_se_bot_left = cubit.create_vertex(v_setop_curve_list[0][0].coordinates()[0], -ly/2.0, 0)
v_se_bot_right = cubit.create_vertex(lx/2, -ly/2.0, 0)
se_curves.append(cubit.create_curve(v_setop_curve_list[0][0], v_se_bot_left))
se_curves.append(cubit.create_curve(v_se_bot_left, v_se_bot_right))
se_curves.append(cubit.create_curve(v_se_bot_right, v_setop_curve_list[-1][-1]))
se_surface = cubit.create_surface(se_curves)
start_id = len(cubit.get_entities('curve'))

# # Interlayer Bottom
interlayer_curves = []
str1 = 'create curve offset curve '
for c in c_se_top_list:
    str1 += ' ' + str(c.id())
distance = 0.0
str1 += ' distance ' + str(distance) + ' normal'
cubit.cmd(str1)
# Now get vertex id's of all of these 2 curves 
# We will end of up 2 vertex id's for each set of curves 
l = len(c_se_top_list)
c_interlayer_bottom_list = list(cubit.get_entities('curve')[start_id:])
# c_interlayer_bottom_list.reverse()
v_interlayer_bottom_vertices = []
for c in c_interlayer_bottom_list:
    v_interlayer_bottom_vertices.append(list(cubit.curve(c).vertices()))
start_id += len(c_interlayer_bottom_list)

interlayer_curves.extend([cubit.curve(c) for c in c_interlayer_bottom_list])

# InterLayer top 
str1 = 'create curve offset curve '
for c in c_se_top_list:
    str1 += ' ' + str(c.id())
distance = dInt
str1 += ' distance ' + str(distance) + ' normal'
cubit.cmd(str1)
c_interlayer_top_list = list(cubit.get_entities('curve')[start_id:])
interlayer_curves.extend([cubit.curve(c) for c in c_interlayer_top_list])
# c_interlayer_top_list.reverse()
v_interlayer_top_vertices = []
for c in c_interlayer_top_list:
    v_interlayer_top_vertices.append(list(cubit.curve(c).vertices()))
start_id += len(c_interlayer_top_list)


if (v_interlayer_top_vertices[0][0].coordinates()[0] < 0):
    str1 = 'move vertex '
    str1 += str(v_interlayer_top_vertices[0][0].id())
    str1 += '0 ' + ' '.join(map(str, v_interlayer_top_vertices[0][0].coordinates()[1:]))
    str1 += ' include_merged'
# Now crate the interface layer
c1 = cubit.create_curve(v_interlayer_top_vertices[0][0], v_interlayer_bottom_vertices[0][0])
c2 = cubit.create_curve(v_interlayer_top_vertices[-1][-1], v_interlayer_bottom_vertices[-1][-1])
interlayer_curves.append(c1)
interlayer_curves.append(c2)
interlayer_surface = cubit.create_surface(interlayer_curves)
start_id += 2

#Metal_top
metal_curves = []
str1 = 'create curve offset curve '
for c in c_se_top_list:
    str1 += ' ' + str(c.id())
distance = dInt
str1 += ' distance ' + str(distance) + ' normal'
cubit.cmd(str1)
c_metal_bottom_list = list(cubit.get_entities('curve')[start_id:])
metal_curves.extend([cubit.curve(c) for c in c_metal_bottom_list])
# c_metal_bottom_list.reverse()
v_metal_bottom_vertices = []
for c in c_metal_bottom_list:
    v_metal_bottom_vertices.append(list(cubit.curve(c).vertices()))
start_id += len(c_metal_bottom_list)


if (v_metal_bottom_vertices[0][0].coordinates()[0] < 0):
    str1 = 'move vertex '
    str1 += str(v_metal_bottom_vertices[0][0].id())
    str1 += '0 ' + ' '.join(map(str, v_metal_bottom_vertices[0][0].coordinates()[1:]))
    str1 += ' include_merged'
v_metalbot_curve_bot = v_metal_bottom_vertices[0][0]
v_metalbot_curve_top = v_metal_bottom_vertices[-1][-1]
v_metal_top_right = cubit.create_vertex(lx/2, ly/4, 0)
v_metal_top_left = cubit.create_vertex(0, ly/4, 0)

c_metal_left = cubit.create_curve(v_metalbot_curve_bot, v_metal_top_left)
c_metal_top = cubit.create_curve(v_metal_top_right, v_metal_top_left)
c_metal_right = cubit.create_curve(v_metalbot_curve_top, v_metal_top_right)
start_id += 3

metal_curves.extend([c_metal_left, c_metal_top, c_metal_right])
metal_surf = cubit.create_surface(metal_curves)

# Now the quartz
v_quartz_top_left = cubit.create_vertex(0, ly/2, 0)
v_quartz_top_right = cubit.create_vertex(lx/2, ly/2, 0)
c_quartz_left = cubit.create_curve(v_metal_top_left, v_quartz_top_left)
c_quartz_top = cubit.create_curve(v_quartz_top_left, v_quartz_top_right)
c_quartz_right = cubit.create_curve(v_quartz_top_right, v_metal_top_right)
c_quartz_bottom = cubit.create_curve(v_metal_top_left, v_metal_top_right)
quartz_curves = [c_quartz_bottom, c_quartz_left, c_quartz_top, c_quartz_right]
quartz_surf = cubit.create_surface(quartz_curves)

# now merge the curves of the two surfaces 
str1 = 'merge curve ' + str(c_metal_top.id()) + ' ' + str(c_quartz_bottom.id())
cubit.cmd(str1)
# --- now create mesh 
surfaces = [se_surface, interlayer_surface, metal_surf, quartz_surf]
# for s in surfaces:
#     # str1 = 'mesh surface '
#     # str2 = 'surface ' 
#     # str1 += str(s.id()) + ' '
#     # str2 += str(s.id()) + ' '
#     # str2 += 'size auto factor 1'
#     # cubit.cmd(str2)
#     # cubit.cmd(str1)
#     str1 = 'surface ' + str(s.id())
#     str1 += 'sizing function type skeleton min_size auto max_size auto max_gradient 1.5 min_num_layers_2d 6 min_num_layers_1d 6 scale 3 time_accuracy_level 2 min_depth 2 max_depth 7 facet_extract_ang 10 max_span_ang_surf 45 max_span_ang_curve 45'
#     cubit.cmd(str1)

#surface 1 sizing function type skeleton min_size auto max_size auto max_gradient 1.5 min_num_layers_2d 6 min_num_layers_1d 6 scale 3 time_accuracy_level 2 min_depth 2 max_depth 7 facet_extract_ang 10 max_span_ang_surf 45 max_span_ang_curve 45 
#refine curve 1 2 numsplit 2 bias 2 depth 2 smooth
#refine curve 12 13 numsplit 2 bias 2 depth 2 smooth
#refine surface 2 numsplit 3 bias 1.0 depth 1 smooth
#refine surface 1 numsplit 2 bias 2 depth 1 smooth