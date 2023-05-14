#!/bin/bash

for iStrnR in 2000 300 30 4
 do
#  sed "s/strn_rt = 2e-2/strn_rt = ${iStrnR}e-5/" < compression_constTrueStrainRate.i > tmp1
  sed "s/strn_rt = 2e-2/strn_rt = ${iStrnR}e-5/" < tension_constTrueStrainRate.i > tmp1
  sed "s/dtmax = 000/dtmax = $(echo "300.0 / $iStrnR" | bc -l)/" < tmp1 > tmp2
  sed "s/end_time = 000/end_time = $(echo "50000.0 / $iStrnR" | bc -l)/" < tmp2 > tmp3
  sed "s/Rate000/RateR${iStrnR}/" < tmp3 > R_Run.i
  ../../electro_chemo_mech2-opt -i R_Run.i
done

rm tmp*
