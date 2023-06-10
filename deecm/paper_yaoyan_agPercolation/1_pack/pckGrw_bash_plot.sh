#!/bin/bash
nMdl=5       # number of models studied

for irt in $(seq 2 $nMdl);
do
###Create folder for splitted calculations

###Edit the lammps file for each calculation
  sed "s/imdl   equal 1/imdl   equal ${irt}/" < lmp_packGrow_plot.in > mr.in

###Edit the submit file to call different number of nodes
#  sed "s/--nodes=1/--nodes=$(($irt))/" < myjob.sh > job0.sh
#  sed "s/-np 20/-np $((100*$irt-80))/" < job0.sh > subjobtest.sh

##Copy files into each directory
  cp mr.in    mrPlot/mr${irt}/
  cp myjob.sh mrPlot/mr${irt}/

###Submit file for calculation
  cd mrPlot/mr${irt}/  
  sbatch myjob.sh
  cd ../../
  rm -f mr.in tmp* job*

done
