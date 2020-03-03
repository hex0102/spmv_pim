#!/bin/bash
source /home/fengsy/intel/mkl/bin/mklvars.sh intel64

if [[ "$1" = "spmm" ]]
then
    for mat_file in `ls -t ./SPMM2/* | sort -V`; do
        echo "Processing file: $mat_file"
        bench=$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '.' -f 1)
        precision=$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '.' -f 2)
        echo -n "$bench"
        echo ".$precision"
        cp $mat_file matA.txt
        read -p "Start s-tui -c --csv-file power.txt and press enter... " -n1 
        ./mkl_spmm
        read -p "Stop s-tui -c --csv-file power.txt and press enter... " -n1 
        mv power.txt power/power_$bench.$precision.txt
    done

elif [[ "$1" = "spmv" ]]
then
    for mat_file in `ls -t ./SPMM2/* | sort -V`; do
        echo "Processing file: $mat_file"
        bench=$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '.' -f 1)
        precision=$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '.' -f 2)
        echo -n "$bench"
        echo ".$precision"
        cp $mat_file matA.txt
        read -p "Start s-tui -c --csv-file power.txt and press enter... " -n1 
        ./mkl_spmv
        read -p "Stop s-tui -c --csv-file power.txt and press enter... " -n1 
        mv power.txt power/power_$bench.$precision.txt
    done

fi

