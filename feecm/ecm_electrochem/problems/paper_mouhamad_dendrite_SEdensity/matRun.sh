#!/bin/bash

gmsh="/usr/bin/gmsh"

# Generate model and mesh	 
sed "s/hLi  = 0.2/hLi  = 19.1478/" < data/t1.geo > data/t20.geo
$gmsh data/t20.geo -2 -o data/t20.msh

# Run the simulation
sed "s/t1/t20/" < microRZ.i > iRun.i
../../deposition-opt -i iRun.i

# Change file name for later data analysis
mv rst/t20_SEAnd_cntLi_0001.csv rst/t20.csv
