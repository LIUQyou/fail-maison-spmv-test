#!/bin/bash
mode=$1
matrix_path=$2
matrix=$(grep -oP '(?<=spm/).*?(?=/)' <<< "$matrix_path")
perf_re_path='/home/quliu/fait-maison-spmv/perf'
./run_float $matrix_path $mode 10&

# check if text file contains the string "flag_start_spmv"
while true; do
    if grep -q "flag_start_spmv" re.txt; then
        PIDOFperf=$(pgrep -f run_float)
        echo "PIDOFperf: ${PIDOFperf}"
        perf stat -p $PIDOFperf -o ${perf_re_path}/perf_${matrix}_${mode}.txt -e cycles,instructions,LLC-load-misses,LLC-loads,LLC-store-misses,LLC-stores,dtlb_load_misses.miss_causes_a_walk,dtlb_load_misses.stlb_hit,dtlb_store_misses.miss_causes_a_walk,dtlb_store_misses.stlb_hit,itlb_misses.miss_causes_a_walk,itlb_misses.stlb_hit,tlb_flush.dtlb_thread,tlb_flush.stlb_any,itlb.itlb_flush
        break
    fi
    sleep 0.1
done
wait
bash process_re.sh ${perf_re_path}/perf_${matrix}_${mode}.txt