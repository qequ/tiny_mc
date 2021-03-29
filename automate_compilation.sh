#!/bin/bash

for ((i = 1; i <= 30; i = i + 1)); do
    echo $((1 + RANDOM % 32768)) >>sample.txt
done

cat parameters.txt | while read p; do cat sample.txt | while read m; do
    make clean &&
        make EXTRA_CFLAGS="${p}" CPPFLAGS="-DPHOTONS=${m}" &&
        ./tiny_mc
done; done

rm -rfv sample.txt
