import json
from queue import PriorityQueue
import sys

output='mr1.txt'
with open(output, 'w') as fp:
  for irt in range(333):
    input='./result/rst'+str(irt)+'.json'
    with open(input, 'rb') as fi:
      data = json.load(fi)
# Abstract active Ag
    num_Ag = data['num_Ag']
    num_active_Ag = data['num_active_Ag']
    wt_active_Ag  = data['active_Ag_percent']
    fp.write('%6d%6d%10.6f\n' %(num_active_Ag,num_Ag,wt_active_Ag))
