#!/bin/bash

gmsh="/usr/bin/gmsh"

for iH in {2..44}
 do
# Generate model and mesh	 
  sed "s/video_input.e/video_l$(echo "$iH - 1" | bc -l).e/" < video_step1.i > tmp
  sed "s/liLength/$(echo "0.5 * $iH" | bc -l)/" < tmp > tmp1
  sed "s/start_time = 0/start_time = $(echo "5 * ($iH-1)" | bc -l)/" < tmp1 > tmp2
  sed "s/end_time = 5/end_time = $(echo "5 * $iH" | bc -l)/" < tmp2 > tmp3
  sed "s/video_output/video_l$(echo "$iH" | bc -l)/" < tmp3 > iR.i
# Run the simulation
   ../../ecm_electrochem-opt -i iR.i

# Change file name for later data analysis
#  mv rst/t${iH}_SEAnd_cntLi_0001.csv rst/t${iH}.csv
   
done

