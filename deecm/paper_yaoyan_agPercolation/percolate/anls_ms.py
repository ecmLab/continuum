import json
from Queue import PriorityQueue
import sys

for fname in range(1,4):
    output='c1_mr_size'+str(fname+1)
    with open(output, 'w') as fp:
      for irt in range(35):
         input='../data/massRatio/rt'+str(irt+1)+'/result/result_size'+str(fname+1)+'.json'
# Load JSON file and variables
         with open(input, 'rb') as fi:
            data = json.load(fi)
# Abstract active NCM
         num_NMC = data['num_NMC']
         num_active_NMC = data['num_active_NMC']
         wt_active_NMC = data['active_NMC_percent']
         fp.write('%6d%6d%10.6f\n' %(num_active_NMC,num_NMC,wt_active_NMC))
