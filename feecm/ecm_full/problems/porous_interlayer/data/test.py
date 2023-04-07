# Python test for cubit curved surface model
import sys 
sys.path.append("/opt/Coreform-Cubit-2020.2/bin")
import cubit
cubit.init(['cubit'])
#####################################################################################################
cubit.cmd('reset')
import math
import numpy as np
ly   = 50       # height of the model, in unit um
dw   = 2.0       # width of the defect located at the top middle, in unit um
dh   = 2.0       # length of the defect, in unit um
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

# # Interlayer Bottom

v_intbot_curve = []
v_intbot_curve_list = []
# v_intbot_curve_ends_list = []
# v_intbot_curve_ends = []

for i in range(npart):
    v_intbot_curve.clear()
    # v_intbot_curve_ends.clear()
    if i > 0:
        vv = cubit.vertex(cubit.get_entities('vertex')[-1])
        print('Vertex {} at {} added to list'.format(vv.id(), vv.coordinates()))
        v_intbot_curve.append(vv)
        newxx = [ax for ax in xx if ax > xsplit[i] and ax <= xsplit[i+1]]
    else:
        newxx = [ax for ax in xx if ax >=xsplit[i] and ax <= xsplit[i+1]]
    for x in newxx:
        if x <= dw:
            y = -dh/2.0*(1.0 + math.cos(math.pi*x/dw))
            yprev = y 
        else: 
            y = yprev
        v_intbot_curve.append(cubit.create_vertex(x,y,0))

    # if len(v_intbot_curve) > 1:
    #     v_intbot_curve_ends.append(v_intbot_curve[0])
    #     v_intbot_curve_ends.append(v_intbot_curve[-1])
    # else:
    #     v_intbot_curve_ends.append(v_intbot_curve_ends[-1][-1])
    #     v_intbot_curve_ends.append(v_intbot_curve[-1])
        
    v_intbot_curve_list.append(v_intbot_curve.copy())    
    # v_intbot_curve_ends_list.append(v_intbot_curve_ends.copy())

# Spline for bottom surface of interlayer
# loop through all vertex list and create curves 
c_interlayer_bot_list = []
c_interlayer_curve_bot = None
for v in v_intbot_curve_list:
    if len(v) > 2:
        str1 = "create curve spline location vertex " + str(v[0].id()) + " to " + str(v[-1].id())
        print(str1)
        cubit.cmd(str1)
        c_interlayer_curve_bot = cubit.curve(cubit.get_entities('curve')[-1])
    else: 
        c_interlayer_curve_bot = cubit.create_curve(v[0], v[-1])
    print('Creating line from {} at {} to {} at {}'.format(v[0].id(), v[0].coordinates(), v[-1].id(), v[-1].coordinates()))
    c_interlayer_bot_list.append(c_interlayer_curve_bot)
# --- Recreate vertex and curve lists 
c_interlayer_bot_list.clear()
v_intbot_curve_list.clear()
v_intbot_curve.clear()
curves = cubit.get_entities('curve')
for c in curves:
    v_intbot_curve.clear()
    c_interlayer_bot_list.append(cubit.curve(c))
    v_intbot_curve_list.append(cubit.curve(c).vertices())

# InterLayer top 

v_inttop_curve = []
k = 0
# now add point at x = 0 since it is not guaranteed
v_inttop_curve.append(cubit.create_vertex(0.0,v_intbot_curve_top.coordinates()[1] + dInt,0))
newx = 0 
newy = 0
yprev = 0
for x in xx:
    if x <= dw:
        y = -dh/2.0 * (1.0 + math.cos(math.pi * x/dw))
        dx = -dw
        dy = -dh/2.0 * math.sin(math.pi * x/dw)
        den = math.sqrt(dx*dx + dy*dy)
        newx = x + dInt/den * dy
        newy = y - dInt/den * dx
        if newx > 0:
            v_inttop_curve.append(cubit.create_vertex(newx,newy,0))
            yprev = newy
    else:
        v_inttop_curve.append(cubit.create_vertex(x,yprev,0))


v_inttop_curve_top = v_inttop_curve[0]
v_inttop_curve_bot = v_inttop_curve[-1]
# Spline for top surface of interlayer
str1 = "create curve spline location vertex " + str(v_inttop_curve_top.id()) + " to " + str(v_inttop_curve_bot.id())
print(str1)
cubit.cmd(str1)
c_interlayer_curve_top = cubit.curve(cubit.get_entities('curve')[-1])

