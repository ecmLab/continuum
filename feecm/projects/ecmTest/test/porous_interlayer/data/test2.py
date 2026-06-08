# Python test for cubit curved surface model
import sys 
sys.path.append("/opt/Coreform-Cubit-2020.2/bin")
sys.path.append("/usr/lib/python3/dist-packages")
import cubit
import math
import numpy as np
cubit.init(['cubit'])
#####################################################################################################
lx   = 50       # width of the model, in unit um
ly   = 50       # height of the model, in unit um
dw   = 2.0       # width of the defect located at the top middle, in unit um
dh   = 4.0       # length of the defect, in unit um
dInt = dw/10      # Initial thickenss of the interlayer, in unit um
nSE  = 51      # Number of discretization points of the cosine shape defect for SE
nSE2 = 50
nLi  = 101      # Number of discretization points of the cosine shape defect for Li
nInt = 101      # Number of discretization points of the cosine shape defect for interlyer

# Interlayer-1 bottom
interlayer_curves = []
v_intbot_curve = []
yprev = 0
x = 0
y = 0
for i in range(0,nSE+1,1):
    x = dw*i/(nSE)
    y = -dh/2.0*(1.0 + math.cos(math.pi*x/dw))
    v_intbot_curve.append(cubit.create_vertex(x,y,0))

xprev = x
yprev = y

for i in range(0, nSE2+1, 1):
    x = xprev + (lx-dw) * i/nSE2
    y = yprev 
    v_intbot_curve.append(cubit.create_vertex(x,y,0))
    print(x)
print(x, y)
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
# now add point at x = 0 since it is not guaranteed
v_inttop_curve.append(cubit.create_vertex(0.0,v_intbot_curve_top.coordinates()[1] + dInt,0))
newx = 0 
newy = 0
for i in range(0,nSE+1,1):
    t = i / (nSE + 1)
    x = dw * t
    y = -dh/2.0 * (1.0 + math.cos(math.pi * x/dw))
    dx = -dw
    dy = -dh/2.0 * math.sin(math.pi * x/dw)
    den = math.sqrt(dx*dx + dy*dy)
    newx = x + dInt/den * dy
    newy = y - dInt/den * dx
    if newx > 0:
        v_inttop_curve.append(cubit.create_vertex(newx,newy,0))

xprev = newx
yprev = newy
for i in range(0, nSE2+1, 1):
    x = xprev + (lx-xprev) * i/nSE2
    y = yprev 
    v_inttop_curve.append(cubit.create_vertex(x,y,0))
print(x, y)

v_inttop_curve_top = v_inttop_curve[0]
v_inttop_curve_bot = v_inttop_curve[-1]
# Spline for top surface of interlayer
str1 = "create curve spline location vertex " + str(v_inttop_curve_top.id()) + " to " + str(v_inttop_curve_bot.id())
print(str1)
cubit.cmd(str1)
c_interlayer_curve_top = cubit.curve(cubit.get_entities('curve')[-1])

