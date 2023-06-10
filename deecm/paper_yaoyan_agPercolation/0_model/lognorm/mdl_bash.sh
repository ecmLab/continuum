#!/bin/bash
nMdl=4       # number of models studied

for irt in $(seq 1 $nMdl);
do
###Create folder for splitted calculations
  mkdir -p ../../data/sz_lognorm/sz${irt}/

###Edit the mdl_all.m file for each calculation
  sed "s/isz    = 4/isz  = ${irt}/" < mdl_main.m > mdl_run.m

###Edit the submit file to call different number of nodes
#  sed "s/--nodes=1/--nodes=$(($irt))/" < myjob.sh > job0.sh
#  sed "s/-np 20/-np $((100*$irt-80))/" < job0.sh > subjobtest.sh

##Copy files into each directory
  cp cmp_nmbr.m ../../data/sz_lognorm/sz${irt}/
  cp cmp_prob.m ../../data/sz_lognorm/sz${irt}/
  cp create_sys.m ../../data/sz_lognorm/sz${irt}/
  cp insrtPtc.m ../../data/sz_lognorm/sz${irt}/
  cp particle_info.m ../../data/sz_lognorm/sz${irt}/
  cp particleDistr.m ../../data/sz_lognorm/sz${irt}/
  cp mdl_run.m ../../data/sz_lognorm/sz${irt}/
  cp mdl_sub.sh ../../data/sz_lognorm/sz${irt}/

###Submit file for calculation
  cd ../../data/sz_lognorm/sz${irt}/
  sbatch mdl_sub.sh
  cd ../../../model/lognorm/
  rm -f mdl_run.m

done
