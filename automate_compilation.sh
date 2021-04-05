#!/bin/bash

for ((i = 1; i <= 60; i = i + 1)); do
    echo $((1 + RANDOM % 32768)) >>sample.txt
done

rm -rfv means.txt
touch means.txt

#Lo que hace esto es leer líneas de parameters.txt que son tipos de compilaciones.
#A cada tipo de compilación le hace correr tiny_mc con los 30 números aleatorios generados.
#Los resultados los guarda en un archivo de texto.
#Después calcula la media de los resultados de la métrica con los 30 números aleatorios generados para la compilación.
cat parameters.txt | while read p; do
    cat sample.txt | while read m; do
        make clean && make EXTRA_CFLAGS="${p}" CPPFLAGS="-DPHOTONS=${m}" && ./tiny_mc >>photons_results.txt;

        if [[ $? -ne 0 ]]; then echo $?; exit 1; fi
    done

    echo python; python3 mean_calc.py photons_results.txt "${p}"; rm -rfv photons_results.txt
done

python3 plot_means.py
rm -rfv sample.txt
