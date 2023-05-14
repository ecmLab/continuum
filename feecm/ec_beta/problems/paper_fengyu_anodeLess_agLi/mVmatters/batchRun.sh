#!/bin/bash

for iSgm in {-2..1}
 do
# Run the simulation
  sed "s/ionic_conductivity = 0.1/ionic_conductivity = 10e${iSgm}/" < dep.i > iRun.i
   ../../../ec_beta-opt -i iRun.i

# Change file name for later data analysis
  mv rst/mdl_anode_current_0001.csv rst/andCrnt_sgm${iSgm}.csv
  mv rst/mdl_anode_potential_0001.csv rst/andPot_sgm${iSgm}.csv
   
done

