#!/bin/bash

for ((i = 1; i <= 100; i = i + 1)); do
    echo $((1 + RANDOM % 65536)) >>sample.txt
done

rm -rfv means.txt
touch means.txt

rm -rfv photons_results.txt
touch photons_results.txt

cat sample.txt | while read m; do
    make clean && make omp_vector CPPFLAGS="-DPHOTONS=${m}" && ./tiny_mc_vectorized_omp;
done

echo python; python3 mean_calc.py photons_results.txt "omp_vec"; rm -rfv photons_results.txt


rm -rfv sample.txt
