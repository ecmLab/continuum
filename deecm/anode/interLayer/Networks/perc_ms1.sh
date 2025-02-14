#!/bin/bash

  cd ./postMg
    for fname in {0..43}
    do

    A=1640000
    B=$((A+10000*fname))

      sed -e '1, 9d' < dump$B.lmp > size$B
    done


