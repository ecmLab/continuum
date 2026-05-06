#!/bin/bash

# For compression test
#for iStrnR in 2000 300 30
# For tension test
for iStrnR in 1
 do
  sed "s/strn_rt = 2e-2/strn_rt = ${iStrnR}e-5/" < compression_constTrueStrainRate.i > tmp1
#  sed "s/strn_rt = 2e-2/strn_rt = ${iStrnR}e-5/" < tension_constTrueStrainRate.i > tmp1
  sed "s/dtmax = 0.5/dtmax = $(echo "3.0 / $iStrnR" | bc -l)/" < tmp1 > tmp2
  sed "s/end_time = 10.0/end_time = $(echo "300.0 / $iStrnR" | bc -l)/" < tmp2 > tmp3
  sed "s/Rate000/RateR${iStrnR}/" < tmp3 > R_Run.i
  ../../ecm_mechanics-opt -i R_Run.i
#  mpirun -np 4 ../../ecm_mechanics-opt -i R_Run.i
done

rm tmp*
