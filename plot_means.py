import matplotlib.pyplot as plt
import numpy
import sys

with open("means.txt") as df:
    datafile = df.read().split("\n")


if sys.argv[1] == "scaling":
    data = list(map(lambda s: s.split(":"), datafile))[:-1]
    num_photons = numpy.array(list(map(int, [d[0] for d in data])))
    photons_per_second = numpy.array(
        list(map(lambda s: float(s.lstrip()), [d[1] for d in data])))

    a, b = numpy.polyfit(numpy.log(num_photons), photons_per_second, 1)

    phots_per_sec_extrapol = list(
        map(lambda x: a * numpy.log(x) + b, num_photons))

    plt.scatter(num_photons, photons_per_second,
                label="promedio de fotones/seg", color="red", s=30)

    plt.plot(num_photons, phots_per_sec_extrapol,
             label='extrapolación de fotones/seg', color="blue")

    plt.xlabel("Número de Fotones")
    plt.ylabel("Fotones/segundo")

    plt.title("Scaling de performance vs tamaño del problema")

    plt.legend()
    plt.show()

else:

    data = list(map(lambda s: s.split(":"), datafile))[:-1]
    flags = [d[0] for d in data]
    means = list(map(lambda s: float(s.lstrip()), [d[1] for d in data]))

    left = [x for x in range(len(flags))]

    plt.bar(left, means, tick_label=flags, width=0.8)

    plt.gcf().subplots_adjust(bottom=0.40)
    plt.xticks(rotation=45)
    plt.xlabel("flags")
    plt.ylabel("Photons/second")
    plt.title("Media de fotones por segundo según flags seteadas al compilador")

    plt.show()
