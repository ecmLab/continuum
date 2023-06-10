import json
import sys
from Queue import PriorityQueue

iNcm = 1
imr  = 5
iLps = 5

output = 'c1_net_mr'+str(imr)+'_ncm'+str(iNcm)+'_lps'+str(iLps)    
input = '../data/massRatio/mr'+str(imr)+'/result/rst_ncm'+str(iNcm)+'_lps'+str(iLps)+'.json'

# Load JSON file and variables
with open(input, 'rb') as fi:
   data = json.load(fi)

# Extract active NCM
num_NMC = data['num_NMC']
num_active_NMC = data['num_active_NMC']

# Extract connectivity
path = data['path']
with open(output, 'w') as fp:
  fp.write('%s\n' %(num_active_NMC))
  for ipt in range(num_active_NMC):
     path_i = path[str(ipt)]
     tmp = " ".join(str(x) for x in path_i[3])
     fp.write('%s\n' %(tmp))  

