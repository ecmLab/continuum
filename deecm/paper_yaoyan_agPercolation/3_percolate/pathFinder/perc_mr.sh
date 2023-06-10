#!/bin/bash

#cd ../pack/massRatio/
for irt in {1..1}
do
#  cd mr$irt/
  mkdir -p result
  cp read_coord.py   ./result/
  cp path_finder.py  ./result/
  cd ./result

## Split files
  stps=334
  csplit -f ts -s -n 3 ../rst3.lmp /'ITEM: TIMESTEP'/ {$(($stps-2))}

## Rename files and delete heads  
  for fname in {0..9}
  do
    cp ts00$fname size$fname.tmp
    sed -e '1, 9d' < size$fname.tmp > size$fname
  done
  for fname in {10..99}
  do
    cp ts0$fname size$fname.tmp
    sed -e '1, 9d' < size$fname.tmp > size$fname
  done
  for (( fname=100; fname<$stps; fname++ ))
  do
    cp ts$fname size$fname.tmp
    sed -e '1, 9d' < size$fname.tmp > size$fname
  done

## Main code
  for (( fname=0; fname<$stps; fname++ ))
  do
    echo processing: size$fname of mr$irt
    python read_coord.py size$fname adj$fname.json
    python path_finder.py adj$fname.json rst$fname.json
  done

  rm ts* size* 

# cd ../../
done
