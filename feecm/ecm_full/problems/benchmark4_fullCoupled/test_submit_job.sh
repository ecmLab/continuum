#!/bin/bash
#BSUB -P "default"
#BSUB -q usercode_16cpu
#BSUB -R "rusage[usercode=16]"
#BSUB -n 16
#BSUB -J cur0.20.05
#BSUB -o .%J.out
#BSUB -e .%J.err
export PETSC_DIR=$PACKAGES_DIR/petsc
### Variables Setting   ####################################################
WORK="/home/sr2/s1.chakravar/moose_tests/temp"
FILE_SER=login05
SCRATCH=/scratch/$LOGNAME/$LSB_JOBID
### INPUT Copy          ####################################################
if [ -e $SCRATCH ]
   then
   rm -rf $SCRATCH
fi

echo "Make scratch dir SCRATCH"
mkdir -p $SCRATCH
chmod 700 $SCRATCH
cd $SCRATCH

echo "Moving to temporarty scratch directory"


cp -auv "$WORK/"* $SCRATCH/
RSYNCERRCODE=$?

COUNT=1
while [ 0 -ne $RSYNCERRCODE ] && [ $COUNT -le 3 ]
do
  sleep 10
  cp -auv "$WORK/"* $SCRATCH/
  let "COUNT+=1"
done

module load openMPI_v3.1
module load python-3.7.3


cp $FILE_SER:"'$WORK/'" $SCRATCH/

/apps/mpi/gcc/RHEL7/openmpi-3.1.0/bin/mpirun -n 16  /home/sr2/s1.chakravar/github/electro_chemo_mech2/electro_chemo_mech2-opt -i full_model.i
### Result File Copy  ####################################################
cp -auv $SCRATCH/* "$WORK/"
RSYNCERRCODE=$?

COUNT=1
while [ 0 -ne $RSYNCERRCODE ] && [ $COUNT -le 3 ]
do
  sleep 10
  cp -auv $SCRATCH/* "$WORK/"
  RSYNCERRCODE=$?
  let "COUNT+=1"
done

### Remove Scratch file  ####################################################
#if [ 0 -eq $RSYNCERRCODE ]
#then
#   cd ../
#   rm -vrf $SCRATCH
#fi

### EXIT Code           ####################################################
#if [ $VASP_ERRCODE -ne 0 ] && [ $RSYNCERRCODE -eq 0 ]
#then
#   exit $VASP_ERRCODE
#fi

