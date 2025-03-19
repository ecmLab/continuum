#!/bin/bash

for iAg in {-10..-1}
 do
# Run the simulation
  sed "s/LiPotRef = -1.0/LiPotRef = ${iAg}/" < hetero.i > iRunAg.i
   ../../ec_beta-opt -i iRunAg.i

# Change file name for later data analysis
  mv rst/hetero_SE_potential_0001.csv rst/vlt_sePot_Ag${iAg}.csv
  mv rst/hetero_pore_potential_0001.csv rst/vlt_porePot_Ag${iAg}.csv
  mv rst/hetero_ag_potential_0001.csv rst/vlt_agPot_Ag${iAg}.csv
  mv rst/hetero_CC_potential_0001.csv rst/vlt_ccPot_Ag${iAg}.csv
   
done

