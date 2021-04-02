#!/bin/bash

#¿Esto genera 30 números aleatorios y los guarda en un archivo de texto?

for ((i = 1; i <= 30; i = i + 1)); do
    echo $((1 + RANDOM % 32768)) >>sample.txt
done

#No entiendo qué hace esta línea.
rm -rfv means.txt; touch means.txt

#Lo que hace esto es leer líneas de parameters.txt que son tipos de compilaciones. 
#A cada tipo de compilación le hace correr tiny_mc con los 30 números aleatorios generados. 
#Los resultados los guarda en un archivo de texto.
#Mi pregunta es: ¿por qué no poner las extra flags directamente en el archivo parameters.txt?
#Después calcula la media de los resultados de la métrica con los 30 números aleatorios generados para la compilación.
cat parameters.txt | while read p; do
    cat sample.txt | while read m; do
        make clean && make EXTRA_CFLAGS="${p}" CPPFLAGS="-DPHOTONS=${m}" && ./tiny_mc >>photons_results.txt
    done;

    echo python;python3 mean_calc.py photons_results.txt "${p}"; rm -rfv photons_results.txt
done

python3 plot_means.py
#No entiendo qué hace esta línea.
rm -rfv sample.txt
