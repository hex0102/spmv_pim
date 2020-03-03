#!/bin/bash
# source /opt/intel/mkl/bin/mklvars.sh intel64
source /home/fengsy/intel/mkl/bin/mklvars.sh intel64
# mkdir -p power_vlsi19/mat
# mkdir -p power_vlsi19/rmat
# mkdir -p power_vlsi19/mat2
# mkdir -p power_vlsi19/mat3
make spmm_single

rm spmm_time.csv
echo "N, r, NNZ, Mult Phase 1, Mult Phase 2, Malloc Time" >> spmm_time.csv

for mat_file in `ls -t ./SPMM/* | sort -V`; do
    echo "Processing file: $mat_file"
    bench=$(echo "$mat_file" | cut -d '/' -f 3 | cut -d 'm' -f 1)
    # folder=$(echo "$mat_file" | cut -d '/' -f 3)
    N=$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '_' -f 2)
    r="0.$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '_' -f 4 | cut -d '.' -f 2)"
    cp $mat_file matA.txt
    echo -n "$N, $r, " >> spmm_time.csv
    echo -n $(./inner_mkl_spmm | tail -3 | head -1 | cut -d ' ' -f 4) >> spmm_time.csv
    echo -n ", " >> spmm_time.csv
    echo $(./inner_mkl_spmm | tail -1) >> spmm_time.csv
    # mv power.txt power_vlsi19/${bench}txt &
    #read -p "Stop s-tui -c --csv-file power.txt and press enter... " -n1 
    #mv power.txt power/power_$bench.txt
done

