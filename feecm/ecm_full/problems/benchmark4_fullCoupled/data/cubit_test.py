# Cubit test for simple geometry
import sys
import cubit
import math
import numpy as np
#start cubit - this step is key

cubit.init([''])

dw = 0.1
dh = 0.2
nInt = 50
ipt = np.linspace(0,50, nInt)
for i in ipt:
    x = dw - dw*i/(nInt+1)
    y = -dh/2.0 * (1.0 + np.cos(math.pi * x/dw))
    cubit.create_vertex(x,y,0)
# Vertices 

