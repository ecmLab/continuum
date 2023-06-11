#!/bin/bash
nMdl=5       # number of models studied

for irt in $(seq 1 $nMdl);
do
###Create folder for splitted calculations

###Edit the lammps file for each calculation
  sed "s/imdl   equal 1/imdl   equal ${irt}/" < pck_main.in > pck_run.in

###Edit the submit file to call different number of nodes
#  sed "s/--nodes=1/--nodes=$(($irt))/" < myjob.sh > job0.sh
#  sed "s/-np 20/-np $((100*$irt-80))/" < job0.sh > subjobtest.sh

##Copy files into each directory
  cp pck_run.in    ../data/mr_uniform/mr${irt}/
  cp pck_sub.sh    ../data/mr_uniform/mr${irt}/

###Submit file for calculation
  cd ../data/mr_uniform/mr${irt}/  
  sbatch pck_sub.sh
  cd ../../../1_pack/
  rm -f pck_run.in tmp* job*

done
