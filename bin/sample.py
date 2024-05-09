import numpy as np
import mobbrmsd as mo
import sys
import matplotlib.pyplot as plt


def sample_run(a, b, nsample):
    dat = np.empty(nsample)
    for i in range(nsample):
        x = cog.generate(n_apm, n_mol, a, b, n_sample=2)
        ret = mrmsd.run(x[0].reshape([-1, 3]), x[1].reshape([-1, 3]))
        dat[i] = ret.n_eval
    return dat


n_apm = int(sys.argv[1])
n_mol = int(sys.argv[2])
nsample = int(sys.argv[3])

mrmsd = mo.mobbrmsd(molecules={"n_apm": n_apm, "n_mol": n_mol})
cog = mo.coord_generator()
nalpha = 5
nbeta = 1
# ax = np.linspace(0.2, 0.8, nalpha)
ax = np.linspace(0.0, 1.0, nalpha)
bx = [0.0]
# bx = np.linspace(0.0, 1.0, nbeta)
by = np.empty((nbeta, nalpha, 4))
hs = []

for i, b in enumerate(bx):
    hh = []
    for j, a in enumerate(ax):
        dat = sample_run(a, b, nsample)
        print(f"{a:8.5f} {b:8.3f} {np.mean(dat):16.9f}")
        hh.append(dat)
    hs.append(hh)
    print()

for i in range(nbeta):
    ind = np.arange(nalpha) + (0.5 * i / nbeta)
    quartile1, medians, quartile3 = np.percentile(hs[i], [25, 50, 75], axis=1)
    plt.violinplot(hs[i], positions=ind, showextrema=False)
    plt.scatter(ind, medians, marker="o", color="white", s=30, zorder=3)
    plt.vlines(ind, quartile1, quartile3, color="b", linestyle="-", lw=5)
    plt.vlines(
        ind,
        [np.min(h) for h in hs[i]],
        [np.max(h) for h in hs[i]],
        color="b",
        linestyle="-",
        lw=1,
    )
plt.show()
# plt.savefig(f"m{n_mol:03}n{n_apm:03}.png")
# plt.savefig(f"m{n_mol:03}n{n_apm:03}.eps")
