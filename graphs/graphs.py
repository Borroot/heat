import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# 0 = seq/major
# 1 = seq/original
# 2 = openmp/major
# 3 = openmp/reduction
# 4 = openmp/schedule
# 5 = mpi/major,1
# 6 = mpi/major,2
# 7 = mpi/major,3
# 8 = mpi/major,4

data = pd.read_csv('output4.csv')
data = data.sort_values(axis=0, by="N")

std = []
means = []
for col in range(1, len(data.columns)):
    std.append([])
    means.append([])

Groups = data.groupby('N')
timeouts = []
N_values = []
for name, g in Groups:
    N_values.append(name)
    for col in range(1, len(g.columns)):
        mean = g[g.columns[col]].mean()
        if mean != -1:
            means[col-1].append(mean)
            std[col-1].append(g[g.columns[col]].std())

ls = ['-', '--', '-.', ':']
for col in range(0, len(data.columns)-1):
    plt.errorbar(N_values[:len(means[col])], means[col], yerr=std[col], label=data.columns[col+1], ls=ls[col % 4])
plt.xscale("log")
plt.legend()
plt.ylabel("seconds")
plt.xlabel("N")
title = "Comparison of openMP scheduling strategies for different values of N"
plt.title(title)
plt.savefig(title + '.png')
