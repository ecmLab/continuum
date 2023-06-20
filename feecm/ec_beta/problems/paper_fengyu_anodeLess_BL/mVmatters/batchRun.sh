#!/bin/bash

#for iSgm in {-2..1}
for iExc in {-1..1}
#for iVlt in {-5..0}
 do

# generate the input data
#  sed "s/ionic_conductivity = 1/ionic_conductivity = 10e${iSgm}/"  < dep.i > iRun.i
  sed "s/ex_current= 1.3/ex_current= 1.3e${iExc}/" < dep.i > iRun.i
#  sed "s/LiPotRef = 0/LiPot = ${iVlt}/" < dep.i > iRun.i

# Run the simulation
   ../../../ec_beta-opt -i iRun.i

# Change file name for later data analysis
#  mv rst/3defects_anode_current_0001.csv rst/3defects_andCrnt_sgm${iSgm}.csv
#  mv rst/3defects_anode_potential_0001.csv rst/3defects_andPot_sgm${iSgm}.csv

  mv rst/3defects_anode_current_0001.csv rst/3defects_andCrnt_exc${iExc}.csv
  mv rst/3defects_anode_potential_0001.csv rst/3defects_andPot_exc${iExc}.csv

#  mv rst/3defects_anode_current_0001.csv rst/3defects_andCrnt_vlt${iVlt}.csv
#  mv rst/3defects_anode_potential_0001.csv rst/3defects_andPot_vlt${iVlt}.csv

done

