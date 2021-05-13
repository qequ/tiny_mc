#!/bin/bash

for ((i = 1; i <= 64; i = i + 1)); do
    echo $(("1024*i")) >>scaling_sample.txt
done

rm -rfv means.txt
touch means.txt

cat scaling_sample.txt | while read m; do

    for i in {1..30}; do
        make clean && make vector CPPFLAGS="-DPHOTONS=${m}" && ./tiny_mc_vectorized;


    done

    echo python; python3 mean_calc.py photons_results.txt "${m}"; rm -rfv photons_results.txt


done


python3 plot_means.py "scaling"

rm -rfv scaling_sample.txt

