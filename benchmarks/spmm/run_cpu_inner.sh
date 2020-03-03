#!/bin/bash
source /home/fengsy/intel/mkl/bin/mklvars.sh intel64
p=6
rm -f ./run_cpu_inner.log
rm -f matA.txt
max=10
type=$1

if [ "$1" = "spmm_old" ]
then
echo "Time elapsed for Mult Phase1 (s), Mult Phase2 (s), Malloc time (s)" >> run_cpu_inner.log
#    for mat_file in `ls -t ./1M/* | sort -V`; do
#        echo "Processing file: $mat_file"
#        cp $mat_file matA.txt
#        head -3 matA.txt | tail -1
#        echo -n "$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '_' -f 2),,," >> run_cpu_inner.log
#        time1=0
#        time2=0
#        time3=0
#        for i in `seq 1 $max`
#        do
#            IFS=, read var1 var2 var3 var4 <<< $(./inner_mkl_spmm 2>&1 | tail -1)
#            time1=$(echo "scale=7; $time1 + $var1" | bc)
#            time2=$(echo "scale=7; $time2 + $var2" | bc)
#            time3=$(echo "scale=7; $time3 + $var3" | bc)
#        done
#        echo "$(echo "scale=7 ; $time1 / $max" | bc),$(echo "scale=7 ; $time2 / $max" | bc),$(echo "scale=7 ; $time3 / $max" | bc)" >> run_cpu_inner.log
#        #break;
#    done
#
#    for mat_file in `ls -t ./ufl/* | sort -V`; do
#        echo "Processing file: $mat_file"
#        cp $mat_file matA.txt
#        #head -3 matA.txt | tail -1
#        echo -n "$(echo "$mat_file" | cut -d '.' -f 2 | cut -d '/' -f 3 ),,," >> run_cpu_inner.log
#        time1=0
#        time2=0
#        time3=0
#        for i in `seq 1 $max`
#        do
#            IFS=, read var1 var2 var3 <<< $(./inner_mkl_spmm 2>&1 | tail -1)
#            time1=$(echo "scale=7; $time1 + $var1" | bc)
#            time2=$(echo "scale=7; $time2 + $var2" | bc)
#            time3=$(echo "scale=7; $time3 + $var3" | bc)
#        done
#        echo "$(echo "scale=7 ; $time1 / $max" | bc),$(echo "scale=7 ; $time2 / $max" | bc),$(echo "scale=7 ; $time3 / $max" | bc)" >> run_cpu_inner.log
#        #break;
#    done

    for mat_file in `ls -t ./SPMM2/* | sort -V`; do
        echo "Processing file: $mat_file"
        cp $mat_file matA.txt
        #head -3 matA.txt | tail -1
        echo -n "$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '.' -f 1 )." >> run_cpu_inner.log
        echo -n "$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '.' -f 2 ),,," >> run_cpu_inner.log
        time1=0
        time2=0
        time3=0
        for i in `seq 1 $max`
        do
            IFS=, read var1 var2 var3 <<< $(./inner_mkl_spmm 2>&1 | tail -1)
            time1=$(echo "scale=7; $time1 + $var1" | bc)
            time2=$(echo "scale=7; $time2 + $var2" | bc)
            time3=$(echo "scale=7; $time3 + $var3" | bc)
        done
        echo "$(echo "scale=7 ; $time1 / $max" | bc),$(echo "scale=7 ; $time2 / $max" | bc),$(echo "scale=7 ; $time3 / $max" | bc)" >> run_cpu_inner.log
        #break;
    done
elif [ "$1" = "papi" ] 
then
    echo "Number of Phase1 L2_DCM, Phase1 L2_DCA, Phase2 L2_DCM, Phase2 L2_DCA" >> run_cpu_inner.log
    for mat_file in `ls -t ./1M/* | sort -V`; do
        echo "Processing file: $mat_file"
        cp $mat_file matA.txt
        head -3 matA.txt | tail -1
        echo -n "$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '_' -f 2),,," >> run_cpu_inner.log
        num1=0
        num2=0
        num3=0
        num4=0
        for i in `seq 1 $max`
        do
            IFS=, read var1 var2 var3 var4 <<< $(./inner_mkl_papi 2>&1 | tail -1)
            num1=$(echo "scale=0; $num1 + $var1" | bc)
            num2=$(echo "scale=0; $num2 + $var2" | bc)
            num3=$(echo "scale=0; $num3 + $var3" | bc)
            num4=$(echo "scale=0; $num4 + $var4" | bc)
        done
        echo "$(echo "scale=0 ; $num1 / $max" | bc),$(echo "scale=0 ; $num2 / $max" | bc),$(echo "scale=0 ; $num3 / $max" | bc),$(echo "scale=0 ; $num4 / $max" | bc)" >> run_cpu_inner.log
        #break;
    done
elif [ "$1" = "spmv" ] 
then
    echo "" >> run_cpu_inner.log
    for mat_file in `ls -t ./SPMV/* | sort -V`; do
        echo "Processing file: $mat_file"
        cp $mat_file matA.txt
        #head -3 matA.txt | tail -1
        # ./mkl_spmv 1000
        echo -n "$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '_' -f 2 )," >> run_cpu_inner.log
        echo -n "$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '_' -f 4 | cut -d '.' -f 1 )." >> run_cpu_inner.log
        echo -n "$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '_' -f 4 | cut -d '.' -f 2 )," >> run_cpu_inner.log
        echo $(./mkl_spmv 10000000 | tail -1) >> run_cpu_inner.log
        #break;
    done

elif [ "$1" = "spmm" ] 
then
    echo "" >> run_cpu_inner.log
    for mat_file in `ls -t ./SPMM/* | sort -V`; do
        echo "Processing file: $mat_file"
        cp $mat_file matA.txt
        #head -3 matA.txt | tail -1
        # ./mkl_spmv 1000
        echo -n "$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '_' -f 2 )," >> run_cpu_inner.log
        echo -n "$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '_' -f 4 | cut -d '.' -f 1 )." >> run_cpu_inner.log
        echo -n "$(echo "$mat_file" | cut -d '/' -f 3 | cut -d '_' -f 4 | cut -d '.' -f 2 )," >> run_cpu_inner.log
        echo $(./mkl_spmm 1000 | tail -1) >> run_cpu_inner.log
        #break;
    done

fi 