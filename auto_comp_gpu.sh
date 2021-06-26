#!/bin/bash

for ((i = 1; i <= 100; i = i + 1)); do
    echo $((1 + RANDOM % 65536)) >>sample.txt
done

rm -rfv means.txt
touch means.txt

rm -rfv photons_results.txt
touch photons_results.txt

cat sample.txt | while read m; do
    nvcc -O3 -use_fast_math -DPHOTONS=${m} -o tiny_mc_gpu tiny_mc_shared_block.cu  && ./tiny_mc_gpu;
done

echo python; python3 mean_calc.py photons_results.txt "gpu"; rm -rfv photons_results.txt


rm -rfv sample.txt
