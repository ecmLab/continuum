# Python test for cubit curved surface model
import sys 
sys.path.append("/opt/Coreform-Cubit-2020.2/bin")
import cubit
import math
import numpy as np
lx   = 50       # width of the model, in unit um
ly   = 50       # height of the model, in unit um
dw   = 2.0       # width of the defect located at the top middle, in unit um
dh   = 4.0       # length of the defect, in unit um
dInt = dw/10      # Initial thickenss of the interlayer, in unit um
nSE  = 51      # Number of discretization points of the cosine shape defect for SE
nLi  = 101      # Number of discretization points of the cosine shape defect for Li
nInt = 101      # Number of discretization points of the cosine shape defect for interlyer

nSE_mesh = 51 
nLi_mesh = 81

pore_width = dw/2.0
pore_wall = 0.1 

# Create SE geometry


m0SE  = lx/10    # mesh characteristic length for bottom points of SE
m1SE  = dInt/4   # mesh characteristic length for interface points of SE
m0Li  = lx/10    # mesh characteristic length for top points of Li
m1Li  = dInt   # mesh characteristic length for bottom points of Li
m0Int = dInt/2   # mesh characteristic length for top points of interlayer
m1Int = dInt/4   # mesh characteristic length for interface points of interlayer

cubit.init(['cubit'])
v_se_curve_top = None
v_se_curve_bot = None
v_se_curve = []
for i in range(0,nSE+1,1):
    x = dw - dw*i/(nSE)
    y = -dh/2.0*(1.0 + math.cos(math.pi*x/dw))
    v_se_curve.append(cubit.create_vertex(x,y,0))
v_se_curve_top = v_se_curve[0]
v_se_curve_bot = v_se_curve[-1]
v_se_bot_left = cubit.create_vertex(v_se_curve_bot.coordinates()[0], -ly/2.0, 0)
v_se_bot_right = cubit.create_vertex(lx/2, -ly/2.0, 0)
v_se_top_right = cubit.create_vertex(lx/2, v_se_curve_top.coordinates()[1],0)

str1 = "create curve spline location vertex " + str(v_se_curve_top.id()) + " to " + str(v_se_curve_bot.id())
cubit.cmd(str1)
# Curve id se_top_curve
c_se = []
c_se_curve_top = cubit.curve(cubit.get_entities('curve')[-1])
c_se.append(c_se_curve_top)
c_se_left = cubit.create_curve(v_se_curve_bot, v_se_bot_left)
c_se.append(c_se_left)
c_se_bot = cubit.create_curve(v_se_bot_left, v_se_bot_right)
c_se.append(c_se_bot)
c_se_right = cubit.create_curve(v_se_bot_right, v_se_top_right)
c_se.append(c_se_right)
c_se_top = cubit.create_curve(v_se_top_right, v_se_curve_top)
c_se.append(c_se_top)
se_surf = cubit.create_surface(c_se)

#####################################################################################################

# Interlayer-1 bottom
interlayer_curves = []
v_intbot_curve = []
for i in range(0,nSE+1,1):
    x = dw - dw*i/(nSE)
    y = -dh/2.0*(1.0 + math.cos(math.pi*x/dw))
    if x < pore_width:
        v_intbot_curve.append(cubit.create_vertex(x,y,0))
v_intbot_curve_top = v_intbot_curve[0]
v_intbot_curve_bot = v_intbot_curve[-1]

# Spline for bottom surface of interlayer
str1 = "create curve spline location vertex " + str(v_intbot_curve_top.id()) + " to " + str(v_intbot_curve_bot.id())
print(str1)
cubit.cmd(str1)
c_interlayer_curve_bot = cubit.curve(cubit.get_entities('curve')[-1])

# Interlayer top
v_inttop_curve = []
k = 0
for i in range(0,nSE+1,1):
    t = i / (nSE + 1)
    x = dw * (1.0 - t)
    y = -dh/2.0 * (1.0 + math.cos(math.pi * x/dw))
    dx = -dw
    dy = -dh/2.0 * math.sin(math.pi * x/dw)
    den = math.sqrt(dx*dx + dy*dy)
    newx = x + dInt/den * dy
    newy = y - dInt/den * dx
    if newx > 0 and newx <= pore_width:
        v_inttop_curve.append(cubit.create_vertex(newx,newy,0))
        k = k + 1
# now add point at x = 0 since it is not guaranteed
v_inttop_curve.append(cubit.create_vertex(0.0,v_intbot_curve_bot.coordinates()[1] + dInt,0))

v_inttop_curve_top = v_inttop_curve[0]
v_inttop_curve_bot = v_inttop_curve[-1]

# Spline for top surface of interlayer
str1 = "create curve spline location vertex " + str(v_inttop_curve_top.id()) + " to " + str(v_inttop_curve_bot.id())
print(str1)
cubit.cmd(str1)
c_interlayer_curve_top = cubit.curve(cubit.get_entities('curve')[-1])

c_interlayer_left = cubit.create_curve(v_inttop_curve_bot, v_intbot_curve_bot)
c_interlayer_right = cubit.create_curve(v_inttop_curve_top, v_intbot_curve_top)

# create loop for interlayer
c_interlayer = [c_interlayer_curve_bot, c_interlayer_left, c_interlayer_curve_top, c_interlayer_right]
surf_interlayer = cubit.create_surface(c_interlayer)
v_inttop_curve_top = c_interlayer_curve_top.vertices()[0]
v_inttop_curve_bot = c_interlayer_curve_top.vertices()[-1]


# Create Li metal surface 

v_metalbot_curve = []
k = 0
for i in range(0,nSE+1,1):
    t = i / (nSE + 1)
    x = dw * (1.0 - t)
    y = -dh/2.0 * (1.0 + math.cos(math.pi * x/dw))
    dx = -dw
    dy = -dh/2.0 * math.sin(math.pi * x/dw)
    den = math.sqrt(dx*dx + dy*dy)
    newx = x + dInt/den * dy
    newy = y - dInt/den * dx
    if newx > 0 and newx <= pore_width:
        v_metalbot_curve.append(cubit.create_vertex(newx,newy,0))
        k = k + 1
# now add point at x = 0 since it is not guaranteed
v_metalbot_curve.append(cubit.create_vertex(v_inttop_curve_bot.coordinates()[0], v_inttop_curve_bot.coordinates()[1],0))

v_metalbot_curve_top = v_metalbot_curve[0]
v_metalbot_curve_bot = v_metalbot_curve[-1]

# Spline for top surface of metal
str1 = "create curve spline location vertex " + str(v_metalbot_curve_top.id()) + " to " + str(v_metalbot_curve_bot.id())
print(str1)
cubit.cmd(str1)
c_metal_curve_top = cubit.curve(cubit.get_entities('curve')[-1])

new_pore_width = v_metalbot_curve_top.coordinates()[0]

c_metal = []
v_metal_top_right = cubit.create_vertex(new_pore_width, ly/4, 0)
v_metal_top_left = cubit.create_vertex(0, ly/4, 0)
c_metal_left = cubit.create_curve(v_metalbot_curve_bot, v_metal_top_left)
c_metal_top = cubit.create_curve(v_metal_top_right, v_metal_top_left)
c_metal_right = cubit.create_curve(v_metalbot_curve_top, v_metal_top_right)
c_metal = [c_metal_curve_top, c_metal_left, c_metal_top, c_metal_right]
metal_surf= cubit.create_surface(c_metal)


#####################################################################################################

#####################################################################################################

# Interlayer-2 bottom
interlayer_curves = []
v_intbot_curve = []
for i in range(0,nSE+1,1):
    x = dw - dw*i/(nSE)
    y = -dh/2.0*(1.0 + math.cos(math.pi*x/dw))
    if x >= pore_width + pore_wall:
        v_intbot_curve.append(cubit.create_vertex(x,y,0))
v_intbot_curve_top = v_intbot_curve[0]
v_intbot_curve_bot = v_intbot_curve[-1]

# Spline for bottom surface of interlayer
str1 = "create curve spline location vertex " + str(v_intbot_curve_top.id()) + " to " + str(v_intbot_curve_bot.id())
print(str1)
cubit.cmd(str1)
c_interlayer_curve_bot = cubit.curve(cubit.get_entities('curve')[-1])

# Interlayer top
v_inttop_curve = []
k = 0
for i in range(0,nSE+1,1):
    t = i / (nSE + 1)
    x = dw * (1.0 - t)
    y = -dh/2.0 * (1.0 + math.cos(math.pi * x/dw))
    dx = -dw
    dy = -dh/2.0 * math.sin(math.pi * x/dw)
    den = math.sqrt(dx*dx + dy*dy)
    newx = x + dInt/den * dy
    newy = y - dInt/den * dx
    if newx > pore_width and newx <= pore_width + pore_wall:
        v_inttop_curve.append(cubit.create_vertex(newx,newy,0))
        k = k + 1
# now add point at x = 0 since it is not guaranteed
v_inttop_curve.append(cubit.create_vertex(0.0,v_intbot_curve_bot.coordinates()[1] + dInt,0))

v_inttop_curve_top = v_inttop_curve[0]
v_inttop_curve_bot = v_inttop_curve[-1]

# Spline for top surface of interlayer
str1 = "create curve spline location vertex " + str(v_inttop_curve_top.id()) + " to " + str(v_inttop_curve_bot.id())
print(str1)
cubit.cmd(str1)
c_interlayer_curve_top = cubit.curve(cubit.get_entities('curve')[-1])

c_interlayer_left = cubit.create_curve(v_inttop_curve_bot, v_intbot_curve_bot)
c_interlayer_right = cubit.create_curve(v_inttop_curve_top, v_intbot_curve_top)

# create loop for interlayer
c_interlayer = [c_interlayer_curve_bot, c_interlayer_left, c_interlayer_curve_top, c_interlayer_right]
surf_interlayer = cubit.create_surface(c_interlayer)
v_inttop_curve_top = c_interlayer_curve_top.vertices()[0]
v_inttop_curve_bot = c_interlayer_curve_top.vertices()[-1]


# Create Li metal surface 

v_metalbot_curve = []
k = 0
for i in range(0,nSE+1,1):
    t = i / (nSE + 1)
    x = dw * (1.0 - t)
    y = -dh/2.0 * (1.0 + math.cos(math.pi * x/dw))
    dx = -dw
    dy = -dh/2.0 * math.sin(math.pi * x/dw)
    den = math.sqrt(dx*dx + dy*dy)
    newx = x + dInt/den * dy
    newy = y - dInt/den * dx
    if newx > 0 and newx <= pore_width:
        v_metalbot_curve.append(cubit.create_vertex(newx,newy,0))
        k = k + 1
# now add point at x = 0 since it is not guaranteed
v_metalbot_curve.append(cubit.create_vertex(v_inttop_curve_bot.coordinates()[0], v_inttop_curve_bot.coordinates()[1],0))

v_metalbot_curve_top = v_metalbot_curve[0]
v_metalbot_curve_bot = v_metalbot_curve[-1]

# Spline for top surface of metal
str1 = "create curve spline location vertex " + str(v_metalbot_curve_top.id()) + " to " + str(v_metalbot_curve_bot.id())
print(str1)
cubit.cmd(str1)
c_metal_curve_top = cubit.curve(cubit.get_entities('curve')[-1])

new_pore_width = v_metalbot_curve_top.coordinates()[0]

c_metal = []
v_metal_top_right = cubit.create_vertex(new_pore_width, ly/4, 0)
v_metal_top_left = cubit.create_vertex(0, ly/4, 0)
c_metal_left = cubit.create_curve(v_metalbot_curve_bot, v_metal_top_left)
c_metal_top = cubit.create_curve(v_metal_top_right, v_metal_top_left)
c_metal_right = cubit.create_curve(v_metalbot_curve_top, v_metal_top_right)
c_metal = [c_metal_curve_top, c_metal_left, c_metal_top, c_metal_right]
metal_surf= cubit.create_surface(c_metal)


#####################################################################################################
