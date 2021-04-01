#Importamos las librerías necesarias.
import sys
import statistics

#Abrimos el archivo means.txt generado con el modo append (porque le queremos agregar nueva información al final).
#Lo llamamos mf.
with open("means.txt", "a") as mf:

#Abrimos el archivo sys.argv[1] con el modo de sólo lectura. Lo llamamos f.
#¿Qué es el archivo sys.argv[1]?
    with open(sys.argv[1], "r") as f:
#Guardamos en una variable auxiliar llamada data algo.
#¿Ese algo qué sería?
        data = f.read().split("\n")

#Esto guarda en photons_per_second algo.
    photons_per_second = list(map(float, data[:-1]))
#Se escribe en mf la media de los datos.
    mf.write("{}: {}\n".format(
        sys.argv[2], statistics.mean(photons_per_second)))
