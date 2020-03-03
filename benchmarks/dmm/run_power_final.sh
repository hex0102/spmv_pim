#!/bin/bash
# source /opt/intel/mkl/bin/mklvars.sh intel64
source /home/fengsy/intel/mkl/bin/mklvars.sh intel64
# mkdir -p power_vlsi19/mat
# mkdir -p power_vlsi19/rmat
# mkdir -p power_vlsi19/mat2
# mkdir -p power_vlsi19/mat3
mkdir -p dmm_power_micro19/
make matrix_power

for N in 8 16 32 64 128 256 512; do
    echo "Processing matrix: $N"
    # folder=$(echo "$mat_file" | cut -d '/' -f 3)
    # cp $mat_file matA.txt
    ./mkl_dense_matrix $N > /dev/null
    mv power.txt dmm_power_micro19/$N.txt &
    #read -p "Stop s-tui -c --csv-file power.txt and press enter... " -n1 
    #mv power.txt power/power_$bench.txt
done

mkdir -p dmv_power_micro19/
make vector_power

for N in 16 32 64 128 256 512 1024 2048 4096 8192; do
    echo "Processing vector: $N"
    # folder=$(echo "$mat_file" | cut -d '/' -f 3)
    # cp $mat_file matA.txt
    ./mkl_dense_vector $N > /dev/null
    mv power.txt dmv_power_micro19/$N.txt &
    #read -p "Stop s-tui -c --csv-file power.txt and press enter... " -n1 
    #mv power.txt power/power_$bench.txt
done