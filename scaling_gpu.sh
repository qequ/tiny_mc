#!/bin/bash

for ((i = 1; i <= 64; i = i + 1)); do
    echo $(("1024*i")) >>scaling_sample.txt
done

rm -rfv means.txt
touch means.txt

cat scaling_sample.txt | while read m; do

    for i in {1..30}; do
        nvcc -O3 -use_fast_math -DPHOTONS="${m}" -o tiny_mc_gpu tiny_mc_shared_block.cu  && ./tiny_mc_gpu;

    done

    echo python; python3 mean_calc.py photons_results.txt "${m}"; rm -rfv photons_results.txt


done


python3 plot_means.py "scaling"

rm -rfv scaling_sample.txt

