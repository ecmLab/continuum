#!/bin/bash

for fname in {0..48}
do
    A=1640000
    B=$((A+10000*fname))
    python read_coord_perc.py size$B adj_size_perc$B.json
    python path_finder_perc.py adj_size_perc$B.json result_size_perc$B.json 
done


