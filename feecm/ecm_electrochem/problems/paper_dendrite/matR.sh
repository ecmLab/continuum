#!/bin/bash

gmsh="/usr/bin/gmsh"

# Generate model and mesh	 
sed "s/hLi  = 0.2/hLi  = howardtu/" < data/t1.geo > data/tttt.geo
$gmsh data/tttt.geo -2 -o data/tttt.msh

# Run the simulation
sed "s/t1/tttt/" < microRZ.i > iRun.i
../../deposition-opt -i iRun.i

# Change file name for later data analysis
mv rst/tttt_SEAnd_cntLi_0001.csv rst/tttt.csv
