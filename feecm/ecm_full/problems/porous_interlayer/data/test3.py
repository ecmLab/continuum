# Python test for cubit curved surface model
import sys 
sys.path.append("/opt/Coreform-Cubit-2020.2/bin")
import cubit
import math
import numpy as np
cubit.init(['cubit'])
#####################################################################################################

ly   = 5       # height of the model, in unit um
dw   = 2.0       # width of the defect located at the top middle, in unit um
dh   = 2.0       # length of the defect, in unit um
dInt = 0.1      # Initial thickenss of the interlayer, in unit um
nSE  = 101      # Number of discretization points of the cosine shape defect for SE
nSE2 = 50
nLi  = 101      # Number of discretization points of the cosine shape defect for Li
nInt = 101      # Number of discretization points of the cosine shape defect for interlyer
lx   = 2*dw       # width of the model, in unit um
# Interlayer-1 bottom
interlayer_curves = []

yprev = 0
# Splits of the interval for pores
xsplit = [dw/4, dw/4 + 0.1, dw/4+0.2, 3*dw/4+0.2, 3*dw/4+0.3, 5*dw/4+0.3, 5*dw/4+0.4, 7*dw/4+0.4]
xspace = [[dw, dw+2]]
npart = len(xsplit) -1
xx = [dw*i/nSE for i in range(0,nSE+1,1)]
xx2 = [lx]
xx = xx + xx2
xx.sort()
x1 = 0
y1 = 0
yprev = 0
v_intbot_curve = []
for x in xx:
    if x <= dw:
        y = -dh/2.0*(1.0 + math.cos(math.pi*x/dw))
        yprev = y 
    else: 
        y = yprev
    v_intbot_curve.append(cubit.create_vertex(x,y,0))
str1 = "create curve spline location vertex " + str(v_intbot_curve[0].id()) + " to " + str(v_intbot_curve[-1].id())
print(str1)
cubit.cmd(str1)
c_interlayer_curve_bot = cubit.curve(cubit.get_entities('curve')[-1])
