#!/bin/bash

for iStrs in 10 8 6
 do
  sed "s/strs    = 1e6/strs    = ${iStrs}e5/" < compression_constTrueStress.i > tmp1
#  sed "s/strs    = 1e6/strs    = ${iStrs}e5/" < tension_constTrueStress.i > tmp1
  sed "s/dtmax = 000/dtmax = $(echo "8.0 / $iStrs" | bc -l)/" < tmp1 > tmp2
  sed "s/end_time = 000/end_time = $(echo "300.0 / $iStrs" | bc -l)/" < tmp2 > tmp3
  sed "s/Stress000/StressC${iStrs}/" < tmp3 > R_Run.i
  ../../electro_chemo_mech2-opt -i R_Run.i
done

rm tmp*
