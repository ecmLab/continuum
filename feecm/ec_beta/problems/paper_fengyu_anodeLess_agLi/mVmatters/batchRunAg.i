#!/bin/bash

for iAg in {0..5}
 do
# Run the simulation
  sed "s/LiPotElectrode = 1.35/LiPotElectrode = ${iAg}/" < depAg.i > iRunAg.i
   ../../../ec_beta-opt -i iRunAg.i

# Change file name for later data analysis
  mv rst/mdlAg_anode_current_0001.csv rst/andCrnt_Ag${iAg}.csv
  mv rst/mdlAg_anode_potential_0001.csv rst/andPot_Ag${iAg}.csv
   
done

