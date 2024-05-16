# -*- coding: utf-8 -*-
"""
Created on Tue May  7 09:05:52 2024

@author: fdceme
"""

import json
from queue import PriorityQueue
import sys
import math
import numpy as np
from scipy import interpolate

args = sys.argv
input1 = args[1]
input2 = args[2]
#input1 = 'adj_size1720000.json'
#input2 = 'result_size1720000.json'
# Load JSON file and variables
with open(input2, 'rb') as fp:
    data = json.load(fp)

with open(input1, 'rb') as fp1:
    data1 = json.load(fp1)

adj_data = data1['adj']
num_Ag = data1['num_Ag']
num_C = data1['num_C']
num_Top = data1['num_Top']
num_particle = data1['num_particle']
vertical_dist_data = data1['vertical_dist']
particle_rads_data = data1['particle_rads']
particle_type_data = data1['type']
SYS_thk = data1['thickness']
spath=data['path']

vertical_dist = {int(key): value for key, value in vertical_dist_data.items()}
particle_rads = {int(key): value for key, value in particle_rads_data.items()}
particle_type = {int(key): value for key, value in particle_type_data.items()}


# 1. Define global parameters
pi      = 3.1415926
# Faraday constant, in unit s*A/mol
cn_F    = 96485.3329
# Combined constant: R*T/F, in unit Volt (T is temperature at 298K, F is Faraday constant)
RT_F    = 1/38.87677

# Combined constant: F^2/RT, in unit sS/mol (T is temperature at 298K, F is Faraday constant)
F2_RT   = (cn_F)**2/(298*8.314463)
# Molar concentration of Lithium
CC      = 0.028

DAg     = 2e-12
DC      = 1.5E-10
DT      = 2E-12

Ag_cond = F2_RT*CC*DAg/10000
C_cond  = F2_RT*CC*DC/10000
T_cond  = F2_RT*CC*DAg/10000

sgm1=[]
for s in range(0,len(particle_type)):
    if particle_type[s] == 'Ag':
        sgm1.append(Ag_cond)
    elif particle_type[s] == 'C':
        sgm1.append(C_cond)
    elif particle_type[s] == 'Top':
        sgm1.append(T_cond)
    elif particle_type[s] == 'Target':
        sgm1.append(T_cond)
        


# 4. Define parameters for the interaction between NMC and LPS
# The exchange current density of NMC with LPS, in unit mA/cm^2 
# (this data is from NMC in liqud electrolyte, subject to change)
I_exc   = 1.0

# 5. System parameters
# The current density used in the cell, in unit mA/cm^2
I0      = 0.05
# The area of the model, in unit um^2
A0      = 6400

sgm     = 0.3
# relative conductivity of GB to bulkLPS: lammada_1 = 2*Sgm_bulkLPS*Thickness_GB/Sgm_GBLPS, in unit um
lmd_1   = 0.001
sgm     = sgm*1e-7
# The total current that path through the model, in unit A
Imdl    = I0*A0*1e-11
# Constant to compute charge-transfer resistance of NMC in LPS, in unit Ohm/um^2
lmd_2   = RT_F / I_exc * 1e6

Om_path = []
Omv_tot  = 0.

for s in range(0,len(spath)):
    path=spath[str(s)][3]
# Compute charge-transfer resistance, resistance from Li+ diffusion in NCM is neglicted 
    R_c   = particle_rads[path[0]]
    tmp1  = adj_data[str(path[0])]
    tmp2  = tmp1[str(path[1])]
    Dlt_c = abs(tmp2[1])
    if len(path) == 2:
        R_Nc = R_c
    else:
        R_N  = particle_rads[path[1]]
        R_Nc = R_c*R_N/(R_c + R_N)
        
    Om_c  = lmd_2/(pi * 2*R_Nc*Dlt_c)
    
    if len(path) == 2:
        Om_b  = 0.
        Om_g  = 0.
    else:
        R_1    = particle_rads[path[-2]]
        R_2    = particle_rads[path[-3]]
        adj1   = adj_data[str(path[-2])]
        tmp    = adj1[str(path[-1])]
        Dlt_01 = abs(tmp[1])
        tmp    = adj1[str(path[-3])]
        Dlt_12 = abs(tmp[1])            
        R_12   = R_1*R_2/(R_1 + R_2)
        cn1    = math.sqrt(1 - 2*Dlt_01/R_1)
        cn2    = math.sqrt(1 - 2*R_12*Dlt_12/(R_1*R_1))
        Om_g   = 1.0/(2*R_1*Dlt_01)
        Om_b   = math.log((1+cn1)/(1-cn1) * (1+cn2)/(1-cn2))/R_1           

        if len(path) > 3:
#                print(path)
            for ilps in range(1, len(path)-2):
                s_i      = sgm1[path[ilps]]
                R_i      = particle_rads[path[ilps]]
                R_im     = particle_rads[path[ilps-1]]
                R_ip     = particle_rads[path[ilps+1]]
                adji     = adj_data[str(path[ilps])]
                tmp      = adji[str(path[ilps-1])]
                Dlt_im   = abs(tmp[1])
                tmp      = adji[str(path[ilps+1])]
                Dlt_ip   = abs(tmp[1])
                R_imi    = R_im*R_i/(R_im + R_i)
                R_ipi    = R_ip*R_i/(R_ip + R_i)
                cn1      = math.sqrt(1 - 2*R_imi*Dlt_im/(R_i*R_i))
                cn2      = math.sqrt(1 - 2*R_ipi*Dlt_ip/(R_i*R_i))
    #                  print([path[ilps],R_i, R_ip, Dlt_ip])
                Om_g     = Om_g + (1.0/R_i + 1.0/R_ip)/(2*Dlt_ip*s_i)
                Om_b     = Om_b + math.log((1+cn1)/(1-cn1) * (1+cn2)/(1-cn2))/(R_i*s_i)

    Om_g     = lmd_1*Om_g/(2*pi)
    Om_b     = Om_b/(2*pi)
    Om_pathi = Om_g + Om_b + Om_c
    Om_path.append(Om_pathi)
    
    if Om_pathi > 0:
        Omv_tot = Omv_tot + 1.0/Om_pathi

A=1/Omv_tot

A=1/Omv_tot
B=A/64000000
print(B)