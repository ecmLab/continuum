#!/bin/bash
# Job name:
#SBATCH --job-name=NMC_LPS
#
# Partition:
#SBATCH --partition=tier3
#
# Processors:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=4
#
# Wall clock limit:
#SBATCH --time=02:30:00
#SBATCH --mem=20g

nMdl=1        # number of models studied
iNcm=1        # NMC size studied, iNcm=1 for 5.0um, iNcm=2 for 10um, iNcm=3 for 12um
iLps=2        # LPS size studied, iLps=3 for 3um, iLps=4 for 4um
sxy=50        # area length studied
sz=100         # initial height

for irt in $(seq 1 $nMdl);
do
###Create folder for splitted calculations
#  mkdir mr_ncm${iNcm}_lps${iLps}_${irt}

###Edit the lammps file for each calculation
  sed "s/index000/index ${irt}/" < lmp_mr.in > tmp0.in
  sed "s/iNcm   equal 1/iNcm   equal ${iNcm}/" < tmp0.in > tmp1.in
  sed "s/iLps   equal 3/iLps   equal ${iLps}/" < tmp1.in > tmp2.in
  sed "s/sxy    equal 80/sxy    equal ${sxy}/" < tmp2.in > tmp3.in
  sed "s/sz     equal 83/sz     equal ${sz}/" < tmp3.in > mr.in  

###Edit the submit file to call different number of nodes
#  sed "s/--nodes=1/--nodes=$(($irt))/" < myjob.sh > job0.sh
#  sed "s/-np 20/-np $((100*$irt-80))/" < job0.sh > subjobtest.sh

##Copy files into each directory
  cp mr.in    massratio/mr${irt}/
  cp myjob.sh massratio/mr${irt}/ 

###Submit file for calculation
  cd massratio/mr${irt}/  
  sbatch myjob.sh
  cd ../../
  rm -f mr.in tmp* job*

done
