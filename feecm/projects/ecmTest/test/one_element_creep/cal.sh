#!/bin/bash

for iT in 198 298 398
 do
  sed "s/temperature = 000/temperature = $iT/" < tension_constTrueStrainRate.i > tmp 
  for iStrnR in 2000 300 30 4
   do
    sed "s/value = 000/value = '1.0e-3 * ${iStrnR}e-5 * exp(${iStrnR}e-5*t)'/" < tmp > tmp1
    sed "s/dtmax = 000/dtmax = $(echo "300.0 / $iStrnR" | bc -l)/" < tmp1 > tmp2
    sed "s/end_time = 000/end_time = $(echo "50000.0 / $iStrnR" | bc -l)/" < tmp2 > tmp3
    sed "s/Rate000/RateT${iT}R${iStrnR}/" < tmp3 > TR_Run.i
    ../../electro_chemo_mech2-opt -i TR_RUN.i
  done
done

rm tmp*
