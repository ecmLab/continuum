# -*- coding: utf-8 -*-
"""
Created on Tue May  7 09:05:52 2024

@author: fdceme
"""

import json
from queue import PriorityQueue
import sys

#args = sys.argv
#input1 = args[1]
#input2 = args[2]
# output = args[2]
input1 = 'adj_size1720000.json'
input2 = 'result_size1720000.json'
# Load JSON file and variables
with open(input2, 'rb') as fp:
    data = json.load(fp)

with open(input1, 'rb') as fp1:
    data1 = json.load(fp1)


Length_part = data1['particle_rads']
Path=data['path']
particle_type_data = data1['type']
particle_type = {int(key): value for key, value in particle_type_data.items()}
len_part = {int(key): value for key, value in Length_part.items()}

Ag_cond=63
C_cond=0.1276
T_cond=63


Area_part=[]
Len_part=[]
for key in range(0,len(len_part)):
    Area_part.append(3.141592*len_part[key]*len_part[key])
    Len_part.append(len_part[key])

Area_part[-1]=Area_part[-2]
Cond=[]
for s in range(0,len(Len_part)):
    if particle_type[s] == 'Ag':
        Cond.append(Ag_cond)
    elif particle_type[s] == 'C':
        Cond.append(C_cond)
    elif particle_type[s] == 'Top':
        Cond.append(T_cond)
    elif particle_type[s] == 'Target':
        Cond.append(T_cond)
R1=[]

for key in Path.keys():
    Case1=Path[key][3]
    Res=[]
    for i in Case1:
        Res.append(Len_part[i]/(Area_part[i]*Cond[i]))
    R1.append(sum(Res))
    
R2=[]

for j in range(0,len(R1)):
    R2.append(1/R1[j])

Rf=1/sum(R2)

print(Rf)