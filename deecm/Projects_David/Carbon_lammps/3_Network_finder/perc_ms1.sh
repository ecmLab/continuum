#!/bin/bash
# Job name:
#SBATCH --job-name=NMC_LPS
#
# Partition:
#SBATCH --partition=tier3
#
# Processors:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=04:30:00
#SBATCH --mem=20g

for cc in {1..5}
do
  cd ../New_loglog/LPS25_0.1_NMC45_0.$cc/massratio/
  for irt in {1..25}
  do
    cd mr$irt/
    mkdir -p result
    cp ../../../../../Percolate/read_coord.py   ./result/
    cp ../../../../../Percolate/path_finder.py  ./result/
    cd ./result

    for fname in {1..1}
    do
      csplit -f ts -s ../rst.lmp /'ITEM: TIMESTEP'/ {5}
      aa=$(ls -dq ts* | wc -l)
      lst=$(($aa-1))
      cp ts0$lst size$fname.lmp
      sed -e '1, 9d' < size$fname.lmp > size$fname

      echo processing: size$cc of mass$irt
      python read_coord.py size$fname adj_size$fname.json
      python path_finder.py adj_size$fname.json result_size$fname.json

      rm ts* size$fname.lmp
    done
  cd ../../
  done
cd ../../../../Percolate/
done
