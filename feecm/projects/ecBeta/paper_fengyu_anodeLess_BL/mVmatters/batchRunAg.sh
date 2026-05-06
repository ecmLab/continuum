#!/bin/bash

for iAg in {-3..-1}
 do
# Run the simulation
  sed "s/LiPotRef = -0.35/LiPotRef = ${iAg}/" < depAg.i > iRunAg.i
   ../../../ec_beta-opt -i iRunAg.i

# Change file name for later data analysis
  mv rst/4Ag3defects_anode_current_0001.csv rst/4Ag3defects_andCrnt_Ag${iAg}.csv
  mv rst/4Ag3defects_anode_potential_0001.csv rst/4Ag3defects_andPot_Ag${iAg}.csv
   
done

