#!/bin/bash
nLi=0
nEn=0
for iLi in $(seq -1.1 0.2 1.1)
do
# Run the simulation
  iSgmLi=$(bc -l <<< "e(${iLi}*l(10))")
#  echo $iSgmLi
  sed "s/ionic_conductivity = 0.8/ionic_conductivity = ${iSgmLi}/" < agPore.i > tmp
  ((nLi++))

    for iEn in $(seq 1 0.2 3.7)
    do
      iSgmEn=$(bc -l <<< "e(${iEn}*l(10))")
      sed "s/electronic_conductivity = 172/electronic_conductivity = ${iSgmEn}/" < tmp > iRunAg.i
      ((nEn++))
 
   ../../ec_beta-opt -i iRunAg.i

# Change file name for later data analysis
  mv rst/agPore_SE_potential_0001.csv rst/sgm_sePot_Li${nLi}En${nEn}.csv
  mv rst/agPore_pore_potential_0001.csv rst/sgm_porePot_Li${nLi}En${nEn}.csv
  mv rst/agPore_ag_potential_0001.csv rst/sgm_agPot_Li${nLi}En${nEn}.csv
  mv rst/agPore_CC_potential_0001.csv rst/sgm_ccPot_Li${nLi}En${nEn}.csv
  mv rst/agPore_zeroX1_potential_0001.csv rst/sgm_z1Pot_Li${nLi}En${nEn}.csv
  mv rst/agPore_zeroX2_potential_0001.csv rst/sgm_z2Pot_Li${nLi}En${nEn}.csv
  mv rst/agPore_zeroX3_potential_0001.csv rst/sgm_z3Pot_Li${nLi}En${nEn}.csv
#  echo $iSgmEn

   done

done

