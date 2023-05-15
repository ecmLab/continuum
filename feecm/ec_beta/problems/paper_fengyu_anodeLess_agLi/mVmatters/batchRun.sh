#!/bin/bash

#for iPrs in {1..3}
for iExc in {-1..1}
 do
# Run the simulation
#  sed "s/ionic_conductivity = 0.1/ionic_conductivity = 10e${iSgm}/" < dep.i > iRun.i
  sed "s/ex_current= 13/ex_current= 10e${iExc}/" < dep.i > iRun.i
#  sed "s/LiPotElectrode = -0.135/LiPotElectrode = $(echo "-0.135 * $iPrs" | bc -l)/" < dep.i > iRun.i
   ../../../ec_beta-opt -i iRun.i

# Change file name for later data analysis
  mv rst/mdl_anode_current_0001.csv rst/andCrnt_exc${iExc}.csv
  mv rst/mdl_anode_potential_0001.csv rst/andPot_exc${iExc}.csv
   
done

