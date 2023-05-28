#!/bin/bash
cd ../pack/massRatio/
for irt in {1..5}
do
  cd mr_ncm1_lps1_$irt/
  mkdir -p result
  cp ../../../percolate/read_coord.py   ./result/
  cp ../../../percolate/path_finder.py  ./result/
  cd ./result

  for fname in {1..1}
  do
    csplit -f ts -s ../pck_size1_$fname.lmp /'ITEM: TIMESTEP'/ {*}
    aa=$(ls -dq ts* | wc -l)
    lst=$(($aa-1))
    cp ts0$lst size$fname.lmp
    sed -e '1, 9d' < size$fname.lmp > size$fname

    echo processing: size$fname of mass$irt
    python read_coord.py size$fname adj_size$fname.json
    python path_finder.py adj_size$fname.json result_size$fname.json

    rm ts* size$fname.lmp
  done

 cd ../../
done
