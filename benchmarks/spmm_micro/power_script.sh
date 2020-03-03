

sleep 10
E1=$(cat /sys/class/powercap/intel-rapl/intel-rapl\:0/energy_uj)
T1=$(date +%s%N | cut -b1-13)
#echo $T1
sleep 10
E2=$(cat /sys/class/powercap/intel-rapl/intel-rapl\:0/energy_uj)
T2=$(date +%s%N | cut -b1-13)
#echo $T2
T3=$(echo "scale=2;($T2-$T1)" | bc)
echo $T3
PWR=$(echo "scale=2;($E2-$E1)/($T3)/1000" | bc)
echo $PWR > power.txt
