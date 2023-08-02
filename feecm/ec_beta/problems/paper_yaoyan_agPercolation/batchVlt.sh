#!/bin/bash

for iAg in {-10..-1}
 do
# Run the simulation
  sed "s/LiPotRef = -1.0/LiPotRef = ${iAg}/" < agPore.i > iRunAg.i
   ../../ec_beta-opt -i iRunAg.i

# Change file name for later data analysis
  mv rst/agPore_SE_potential_0001.csv rst/vlt_sePot_Ag${iAg}.csv
  mv rst/agPore_pore_potential_0001.csv rst/vlt_porePot_Ag${iAg}.csv
  mv rst/agPore_ag_potential_0001.csv rst/vlt_agPot_Ag${iAg}.csv
  mv rst/agPore_CC_potential_0001.csv rst/vlt_ccPot_Ag${iAg}.csv
  mv rst/agPore_Y1Sec1_potential_0001.csv rst/vlt_z1Pot_Ag${iAg}.csv
  mv rst/agPore_Y1Sec2_potential_0001.csv rst/vlt_z2Pot_Ag${iAg}.csv
  mv rst/agPore_Y1Sec3_potential_0001.csv rst/vlt_z3Pot_Ag${iAg}.csv
  mv rst/agPore_center_potential_0001.csv rst/vlt_z0Pot_Ag${iAg}.csv
   
done

