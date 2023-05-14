#!/bin/bash

gmsh="/usr/bin/gmsh"

for iH in {2..20..2}
 do
# Generate model and mesh	 
  sed "s/xLi   = 5/xLi  = ${iH}/" < data/t0.geo > data/t${iH}.geo
  $gmsh data/t${iH}.geo -2 -o data/t${iH}.msh

# Run the simulation
  sed "s/lLi24/t${iH}/" < phase1.i > tmp.i
  sed "s/end_time = 200/end_time = $(echo "10 * $iH" | bc -l)/" < tmp.i > iRun.i
   ../../ecm_electrochem-opt -i iRun.i

# Change file name for later data analysis
#  mv rst/t${iH}_SEAnd_cntLi_0001.csv rst/t${iH}.csv
   
done

