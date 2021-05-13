# Importamos las librerías necesarias.
import sys
import statistics

with open("means.txt", "a") as mf:

    with open(sys.argv[1], "r") as f:
        data = f.read().split("\n")

    photons_per_second = list(map(float, data[:-1]))
    mf.write("{}: {}\n".format(
        sys.argv[2], statistics.mean(photons_per_second)))
