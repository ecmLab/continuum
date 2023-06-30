#!/bin/bash
nMdl=4       # number of models studied

for irt in $(seq 1 $nMdl);
do
###Create folder for splitted calculations

###Edit the lammps file for each calculation
  sed "s/imdl   equal 1/imdl   equal ${irt}/" < pck_main_sz.in > pck_run.in

###Edit the submit file to call different number of nodes
#  sed "s/--nodes=1/--nodes=$(($irt))/" < myjob.sh > job0.sh
#  sed "s/-np 20/-np $((100*$irt-80))/" < job0.sh > subjobtest.sh

##Copy files into each directory
  cp pck_run.in       ../data/sz_lognorm_062723/sz${irt}/
  cp pck_sub_sz.sh    ../data/sz_lognorm_062723/sz${irt}/

###Submit file for calculation
  cd ../data/sz_lognorm_062723/sz${irt}/  
  sbatch pck_sub_sz.sh
  cd ../../../1_pack/
  rm -f pck_run.in tmp* job*

done
