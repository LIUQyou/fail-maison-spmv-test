#!/bin/bash
input_file=$1
dtb_load_misses=$(cat $input_file | grep "dtlb_load_misses.miss_causes_a_walk " | awk '{print $1}' | sed -e $'s/,//g')
dtb_store_misses=$(cat $input_file | grep "dtlb_store_misses.miss_causes_a_walk " | awk '{print $1}' | sed -e $'s/,//g')
itb_misses=$(cat $input_file | grep "itlb_misses.miss_causes_a_walk " | awk '{print $1}' | sed -e $'s/,//g')
instructions=$(cat $input_file | grep "instructions " | awk '{print $1}' | sed -e $'s/,//g')

echo $dtb_load_misses $dtb_store_misses $itb_misses $instructions
python -c "print('MPKI: ', ($dtb_load_misses+$dtb_store_misses+$itb_misses)/($instructions/1000))"

dtb_load_misses_duration=$(cat $input_file | grep "dtlb_load_misses.walk_duration " | awk '{print $1}' | sed -e $'s/,//g')
dtb_store_misses_duration=$(cat $input_file | grep "dtlb_store_misses.walk_duration " | awk '{print $1}' | sed -e $'s/,//g')
itb_misses_duration=$(cat $input_file | grep "itlb_misses.walk_duration  " | awk '{print $1}' | sed -e $'s/,//g')
cycles=$(cat $input_file | grep "cycles " | awk '{print $1}' | sed -e $'s/,//g')
echo $dtb_load_misses_duration $dtb_store_misses_duration $itb_misses_duration $cycles
python -c "print('TLB miss overhead (%): ', 100*($dtb_load_misses_duration+$dtb_store_misses_duration+$itb_misses_duration)/($cycles))"

llc_load_miss=$(cat $input_file | grep "LLC-load-misses " | awk '{print $1}' | sed -e $'s/,//g')
llc_load=$(cat $input_file | grep "LLC-loads " | awk '{print $1}' | sed -e $'s/,//g')
llc_store_miss=$(cat $input_file | grep "LLC-store-misses " | awk '{print $1}' | sed -e $'s/,//g')
llc_store=$(cat $input_file | grep "LLC-stores " | awk '{print $1}' | sed -e $'s/,//g')
echo $llc_load_miss $llc_load $llc_store_miss $llc_store
python -c "print('LLC miss rate (%): ', 100*($llc_load_miss+$llc_store_miss)/($llc_load+$llc_store+$llc_load_miss+$llc_store_miss))"
