#!/bin/bash
nMdl=5       # number of models studied

for irt in $(seq 1 $nMdl);
do
###Create folder for splitted calculations

###Edit the lammps file for each calculation
  sed "s/imdl   equal 1/imdl   equal ${irt}/" < lmp_grow.in > grw.in

###Edit the submit file to call different number of nodes
#  sed "s/--nodes=1/--nodes=$(($irt))/" < myjob.sh > job0.sh
#  sed "s/-np 20/-np $((100*$irt-80))/" < job0.sh > subjobtest.sh

##Copy files into each directory
  cp grw.in    massRatio/mr${irt}/
  cp mygrw.sh massRatio/mr${irt}/

###Submit file for calculation
  cd massRatio/mr${irt}/  
  sbatch mygrw.sh
  cd ../../
  rm -f grw.in tmp* job*

done
