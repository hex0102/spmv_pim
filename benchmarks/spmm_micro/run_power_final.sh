#!/bin/bash
# source /opt/intel/mkl/bin/mklvars.sh intel64
source /home/fengsy/intel/mkl/bin/mklvars.sh intel64
# mkdir -p power_vlsi19/mat
# mkdir -p power_vlsi19/rmat
# mkdir -p power_vlsi19/mat2
# mkdir -p power_vlsi19/mat3
mkdir -p power_micro19
make spmm_single_power

for mat_file in `ls -t ./SPMM/* | sort -V`; do
    echo "Processing file: $mat_file"
    bench=$(echo "$mat_file" | cut -d '/' -f 3 | cut -d 'm' -f 1)
    # folder=$(echo "$mat_file" | cut -d '/' -f 3)
    cp $mat_file matA.txt
    ./inner_mkl_spmm_power
    mv power.txt power_micro19/${bench}txt &
    #read -p "Stop s-tui -c --csv-file power.txt and press enter... " -n1 
    #mv power.txt power/power_$bench.txt
done

