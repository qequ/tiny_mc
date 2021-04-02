import matplotlib.pyplot as plt

with open("means.txt") as df:
    datafile = df.read().split("\n")

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
