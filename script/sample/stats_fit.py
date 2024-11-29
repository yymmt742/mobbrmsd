import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import colorsets


def linfit(x, y):

    sxx = np.sum(np.power(x, 2))
    sxy = np.sum(x * y)
    sx = np.sum(x)
    mx = np.mean(x)
    sy = np.sum(y)
    my = np.mean(y)
    a = (sxy - sx * my) / (sxx - sx * mx)
    return a, my - a * mx


df = pd.concat([pd.read_pickle(a, compression="gzip") for a in sys.argv[1:]])
x = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
y = [1, 2, 3, 10, 100, 1000]
ea = np.empty([len(x), len(y)])
print(f"$\\alpha$ &    1 &    2 &    3 &   10 &  100 & 1000 \\\\")
for i, alpha in enumerate(x):
    for j, n_apm in enumerate(y):
        mask = (df["n_apm"] == n_apm) & (df["n_sym"] == 1) & (df["alpha"] == alpha)
        n_mol = df["n_mol"][mask].to_numpy()
        n_mean = df["n_mean"][mask].to_numpy()
        a, b = linfit(n_mol, np.log(n_mean))
        ea[i, j] = np.exp(a)

    print(
        f"{alpha:8.1f} & {ea[i,0]:4.1f} & {ea[i,1]:4.1f} & {ea[i,2]:4.1f} & {ea[i,3]:4.1f} & {ea[i,4]:4.1f} & {ea[i,5]:4.1f} \\\\"
    )


fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(ea, cmap="cividis")
for i in range(ea.shape[0]):
    for j in range(ea.shape[1]):
        ax.text(j, i, f"{ea[i,j]:4.1f}", ha="center", va="center")
plt.savefig("a_n.eps")
# plt.show()
plt.clf()
plt.cla()
plt.close()

x = [0.0, 0.6, 0.9]
y = [3, 4, 5, 6]
ea = np.empty([len(x), len(y)])
print(f"$\\alpha$ &    3 &    4 &    5 &    6 \\\\")
for i, alpha in enumerate(x):
    for j, n_mol in enumerate(y):
        mask = (df["n_apm"] == 10) & (df["n_mol"] == n_mol) & (df["alpha"] == alpha)
        n_sym = df["n_sym"][mask].to_numpy()
        n_mean = df["n_mean"][mask].to_numpy()
        a, b = linfit(np.log(n_sym), np.log(n_mean))
        ea[i, j] = a
    print(
        f"{alpha:8.1f} & {ea[i,0]:4.1f} & {ea[i,1]:4.1f} & {ea[i,2]:4.1f} & {ea[i,3]:4.1f} \\\\"
    )

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(ea, cmap="cividis")
for i in range(ea.shape[0]):
    for j in range(ea.shape[1]):
        ax.text(j, i, f"{ea[i,j]:4.1f}", ha="center", va="center")
plt.savefig("a_S.eps")
# plt.show()
plt.clf()
plt.cla()
plt.close()
