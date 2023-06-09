#!/bin/bash
nMdl=4       # number of models studied

for irt in $(seq 1 $nMdl);
do
###Create folder for splitted calculations

###Edit the lammps file for each calculation
  sed "s/imdl   equal 1/imdl   equal ${irt}/" < grw_main.in > grw_run.in

###Edit the submit file to call different number of nodes
#  sed "s/--nodes=1/--nodes=$(($irt))/" < myjob.sh > job0.sh
#  sed "s/-np 20/-np $((100*$irt-80))/" < job0.sh > subjobtest.sh

##Copy files into each directory
#  cp grw_run.in    ../data/mr_uniform/mr${irt}/
#  cp grw_sub.sh    ../data/mr_uniform/mr${irt}/
  cp grw_run.in    ../data/sz_lognorm_062723/sz${irt}/
  cp grw_sub.sh    ../data/sz_lognorm_062723/sz${irt}/

###Submit file for calculation
  cd ../data/sz_lognorm_062723/sz${irt}/  
#  sbatch grw_sub.sh
  cd ../../../2_grow/
  rm -f grw_run.in tmp* job*

done
