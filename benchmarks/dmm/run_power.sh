#!/bin/bash
source /opt/intel/mkl/bin/mklvars.sh intel64

for N in 4096 8192 16384; do
    echo "Processing file matrix of Dimention $N"
    bench=$(echo "Vector_$N")
    read -p "Start s-tui -c --csv-file power.txt and press enter... " -n1 
    ./mkl_dense_vector $N
    read -p "Stop s-tui -c --csv-file power.txt and press enter... " -n1 
    mv power.txt power/power_$bench.txt
done

