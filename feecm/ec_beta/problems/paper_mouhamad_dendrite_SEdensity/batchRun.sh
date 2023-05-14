#!/bin/bash

gmsh="/usr/bin/gmsh"

for iH in {1..5}
 do

# Generate model and mesh	 
  sed "s/hLi  = 0.2/hLi  = 0.${iH}/" < data/t0.geo > data/t${iH}.geo
  $gmsh data/t${iH}.geo -2 -o data/t${iH}.msh

# Run the simulation
  sed "s/t0/t${iH}/" < microRZ.i > iRun.i
   ../../deposition-opt -i iRun.i

# Change file name for later data analysis
  mv rst/t${iH}_SEAnd_cntLi_0001.csv rst/t${iH}.csv
   
done

