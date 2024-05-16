#!/bin/bash

  cd ../post
    for fname in {0..80}
    do

    A=1720000
    B=$((A+1000*fname))

      sed -e '1, 9d' < rstgr_$B.lmp > size$B
    done


