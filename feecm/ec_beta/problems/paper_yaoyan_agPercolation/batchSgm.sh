#!/bin/bash
nLi=0
for iLi in $(seq -1.1 0.2 1.5)
do
# Run the simulation
  iSgmLi=$(bc -l <<< "e(${iLi}*l(10))")
#  echo $iSgmLi
  sed "s/ionic_conductivity = 0.8/ionic_conductivity = ${iSgmLi}/" < homo.i > iRunAg.i
  ((nLi++))

  ../../ec_beta-opt -i iRunAg.i

# Change file name for later data analysis
  mv rst/homo_SE_potential_0001.csv rst/sgm_sePot_Li${nLi}.csv
  mv rst/homo_ag_potential_0001.csv rst/sgm_agPot_Li${nLi}.csv
  mv rst/homo_CC_potential_0001.csv rst/sgm_ccPot_Li${nLi}.csv

done

