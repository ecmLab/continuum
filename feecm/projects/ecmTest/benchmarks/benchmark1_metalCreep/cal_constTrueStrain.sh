#!/bin/bash

for iStrnR in 2000 300 30
 do
#  sed "s/strn_rt = 4e-5/strn_rt = ${iStrnR}e-5/" < compression_constTrueStrainRate.i > tmp1
  sed "s/strn_rt = 4e-5/strn_rt = ${iStrnR}e-5/" < tension_constTrueStrainRate.i > tmp1
  sed "s/dtmax = 0.5/dtmax = $(echo "300.0 / $iStrnR" | bc -l)/" < tmp1 > tmp2
  sed "s/end_time = 10000.0/end_time = $(echo "50000.0 / $iStrnR" | bc -l)/" < tmp2 > R_Run.i
  mpirun -np 4 ../../../../ecm_test/ecm-opt -i R_Run.i
done

rm tmp*
