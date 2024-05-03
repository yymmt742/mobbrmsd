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
    return 0


def lj(r, r0):
    z = r0 / r
    z6 = np.power(z, 6)
    z12 = np.power(z6, 2)
    return z12 - z6


def energy(x):
    ret = 0.0
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            for k in range(i + 1, x.shape[0]):
                for l in range(x.shape[1]):
                    ret += lj(np.linalg.norm(x[i, j] - x[k, l]), 0.5)
    return 50 * ret


def conflict(x, dtol) -> False:
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            for k in range(i + 1, x.shape[0]):
                for l in range(x.shape[1]):
                    if np.linalg.norm(x[i, j] - x[k, l]) < dtol:
                        return True
    return False


def sample_run(a, b, nhist, nsample, dtol):
    hist = np.zeros(nhist)
    sdtol = dtol * np.sin(0.5 * np.pi * a)
    for i in range(nsample):
        x = cog.generate(n_apm, n_mol, a, b, n_sample=2)
        """
        while True:
            x = cog.generate(n_apm, n_mol, a, b, n_sample=2)
            if conflict(x[0], dtol):
                continue
            if conflict(x[1], dtol):
                continue
            break
        """
        ret = mrmsd.run(x[0].reshape([-1, 3]), x[1].reshape([-1, 3]))
        hist[ret.n_eval - 1] += 1.0
    hsum = np.max([np.sum(hist), 1e-16])
    p25 = percentiles(hist, hsum / 4)
    p50 = percentiles(hist, hsum / 2)
    p75 = percentiles(hist, 3 * hsum / 4)
    ii1 = np.linspace(1, nhist, num=nhist)
    ii2 = np.power(ii1, 2)
    mean = np.sum(hist * ii1) / hsum
    std = np.sqrt(np.max([0.0, np.sum(hist * ii2) / hsum - mean * mean]))
    return hist, p25, p50, p75, mean, std


n_apm = int(sys.argv[1])
n_mol = int(sys.argv[2])
nsample = int(sys.argv[3])
nhist = int(sys.argv[4])
dtol = float(sys.argv[5])

mrmsd = mo.mobbrmsd()
mrmsd.add_molecule(n_apm, n_mol)
cog = mo.coord_generator()
nalpha = 21
nbeta = 3
ax = np.linspace(0.2, 0.8, nalpha)
bx = np.linspace(0.0, 1.0, nbeta)
# ax = np.linspace(0.0125, 0.9875, 79)
by = np.empty((nbeta, nalpha, 4))

for i, b in enumerate(bx):
    for j, a in enumerate(ax):
        hist, p25, p50, p75, mean, std = sample_run(a, b, nhist, nsample, dtol)
        print(
            f"{a:8.5f} {b:8.3f} {p25:16.9f} {p50:16.9f} {p75:16.9f} {mean:16.9f} {std:16.9f}"
        )
        by[i, j, 0] = p25
        by[i, j, 1] = p50
        by[i, j, 2] = p75
        by[i, j, 3] = mean
    print()
# plt.show()
# plt.clf()
for i in range(nbeta):
    plt.fill_between(ax, by[i, :, 0], by[i, :, 2], alpha=0.3)
for i in range(nbeta):
    plt.plot(ax, by[i, :, 3], label=f"beta={bx[i]:f}")
plt.legend()
plt.show()
# plt.savefig(f"m{n_mol:03}n{n_apm:03}.png")
# plt.savefig(f"m{n_mol:03}n{n_apm:03}.eps")
