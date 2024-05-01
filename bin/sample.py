import numpy as np
import mobbrmsd as mo
import sys
import matplotlib.pyplot as plt


def percentiles(x, q):
    c = 0.0
    for i, xi in enumerate(x):
        c += xi
        if c > q:
            return i


def sample_run(a, b, nhist, nsample):
    hist = np.zeros(nhist)
    for i in range(nsample):
        # x = np.arange(3 * n_apm * n_mol).reshape([-1, n_apm, 3])
        # y = np.arange(3 * n_apm * n_mol).reshape([-1, n_apm, 3])
        x, y = cog.generate_pair(n_apm, n_mol, a, b)
        ret = mrmsd.run(x, y)
        hist[ret.n_eval - 1] += 1.0
    hsum = np.sum(hist)
    p25 = percentiles(hist, hsum / 4)
    p50 = percentiles(hist, hsum / 2)
    p75 = percentiles(hist, 3 * hsum / 4)
    ii1 = np.linspace(1, nhist, num=nhist)
    ii2 = np.power(ii1, 2)
    mean = np.sum(hist * ii1) / hsum
    std = np.sqrt(np.sum(hist * ii2) / hsum - mean * mean)
    return hist, p25, p50, p75, mean, std


n_apm = int(sys.argv[1])
n_mol = int(sys.argv[2])
nsample = int(sys.argv[3])
nhist = int(sys.argv[4])

mrmsd = mo.mobbrmsd()
mrmsd.add_molecule(n_apm, n_mol)
cog = mo.coord_generator()
ax = np.linspace(0.025, 1.0, 40)
by = np.empty((2, 40, 4))

for i, b in enumerate([0.0, 1.0]):
    for j, a in enumerate(ax):
        hist, p25, p50, p75, mean, std = sample_run(a, b, nhist, nsample)
        # plt.plot(hist)
        print(
            f"{a:8.3f} {b:8.3f} {p25:16.9f} {p50:16.9f} {p75:16.9f} {mean:16.9f} {std:16.9f}"
        )
        by[i, j, 0] = p25
        by[i, j, 1] = p50
        by[i, j, 2] = p75
        by[i, j, 3] = mean
    print()
# plt.show()
# plt.clf()
plt.fill_between(ax, by[0, :, 0], by[0, :, 2], alpha=0.2)
plt.fill_between(ax, by[1, :, 0], by[1, :, 2], alpha=0.2)
plt.plot(ax, by[0, :, 3], label="beta=0.0")
plt.plot(ax, by[1, :, 3], label="beta=1.0")
plt.legend()
plt.show()
