#!/bin/bash

filename="values.txt"

> "$filename"

for fname in {0..80}
do
    A=1720000
    B=$((A+1000*fname))
    #python read_coord.py size$B adj_size$B.json
    #python path_finder.py adj_size$B.json result_size$B.json
    python resistance_1.py adj_size$B.json result_size$B.json >> "$filename"    
done


