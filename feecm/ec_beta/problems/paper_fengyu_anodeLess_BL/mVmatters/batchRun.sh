#!/bin/bash

#for iSgm in {-2..1}
#for iExc in {-1..1}
for iVlt in {-5..0}
#for iPrs in {1..3}
 do

# generate the input data
#  sed "s/ionic_conductivity = 1/ionic_conductivity = 10e${iSgm}/"  < dep.i > iRun.i
#  sed "s/ex_current= 13/ex_current= 1.3e${iExc}/" < dep.i > iRun.i
  sed "s/LiPotRef = 0/LiPot = ${iVlt}/" < dep.i > iRun.i
#  sed "s/LiPot = -0.135/LiPot = $(echo "-0.135 * $iPrs" | bc -l)/" < dep.i > iRun.i

# Run the simulation
   ../../../ec_beta-opt -i iRun.i

# Change file name for later data analysis
#  mv rst/mdl_anode_current_0001.csv rst/andCrnt_sgm${iSgm}.csv
#  mv rst/mdl_anode_potential_0001.csv rst/andPot_sgm${iSgm}.csv

#  mv rst/mdl_anode_current_0001.csv rst/andCrnt_exc${iExc}.csv
#  mv rst/mdl_anode_potential_0001.csv rst/andPot_exc${iExc}.csv

  mv rst/mdl_anode_current_0001.csv rst/andCrnt_vlt${iVlt}.csv
  mv rst/mdl_anode_potential_0001.csv rst/andPot_vlt${iVlt}.csv

#  mv rst/mdl_anode_potential_0001.csv rst/andPot_prs${iPrs}.csv
#  mv rst/mdl_anode_potential_0001.csv rst/andPot_prs${iPrs}.csv
   
done

