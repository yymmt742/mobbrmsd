import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import colorsets
from pathlib import Path

cg = colorsets.uppl
cm = colorsets.uorg
cl = colorsets.ublk


def logymax(v):
    return 10 ** np.floor(np.log10(v) + 1)


plt.rcParams["font.size"] = 8

df = pd.concat([pd.read_pickle(a, compression="gzip") for a in sys.argv[1:]])
fig = plt.figure()
ax = fig.add_subplot(111)
n_mols = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
for n_mol in n_mols:
    for key, c, mec, mfc, delta, marker, ms, mss, mew, elinewidth in zip(
        ["time_grmsd", "time_mobbrmsd", "time_lprmsd"],
        [cg, cm, cl],
        [None, cm, cl],
        [colorsets.uwht, colorsets.uwht, None],
        [0.15, -0.2, 0.0],
        ["s", "o", "_"],
        [3, 3, 4],
        [2, 3, 4],
        [1.0, 1.0, 1.0],
        [2, 4, 1],
    ):
        t = df[key][(df["n_mol"] == n_mol)]
        x = t[~(t.isnull())]
        med = np.percentile(x, 50)
        n_points = x.shape[0]
        if n_points > 20:
            err = [[med - np.percentile(x, 25)], [np.percentile(x, 75) - med]]
            ax.errorbar(
                n_mol + delta,
                med,
                yerr=err,
                color=c,
                linestyle="",
                capsize=0,
                marker=marker,
                mfc=mfc,
                mec=mec,
                mew=mew,
                elinewidth=elinewidth,
                ms=ms,
            )

        else:
            ax.scatter(
                [n_mol + delta] * n_points,
                x,
                facecolors=c,
                edgecolors="none",
                s=mss,
                lw=0.2,
            )
            ax.scatter(
                [n_mol + delta],
                [med],
                edgecolors=c,
                facecolors=mfc,
                marker=marker,
                s=4 * mss,
            )

ax.set_xticks(n_mols)
ax.set_yscale("log")
ax.set_ylim([0.00001, 100000])
ax.set_yticks([0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000])
ax.set_box_aspect(1.414)
# ax.set_box_aspect(0.707)
ax.tick_params(which="minor", direction="out", axis="y", length=2, width=0.5)
plt.savefig(f"vs_grmsd_time.png")
plt.savefig(f"vs_grmsd_time.eps")
plt.clf()
plt.close()

exit()

for n_mol in n_mols:
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    t = df["mobbrmsd"][(df["n_mol"] == n_mol)]
    xm = t[~(t.isnull())]
    t = df["grmsd"][(df["n_mol"] == n_mol)]
    xg = t[~(t.isnull())]
    ax1.scatter(
        xm,
        xg,
        facecolors="none",
        edgecolors="r",
        lw=0.1,
    )
    t = df["time_mobbrmsd"][(df["n_mol"] == n_mol)]
    tm = t[~(t.isnull())]
    t = df["time_grmsd"][(df["n_mol"] == n_mol)]
    tg = t[~(t.isnull())]
    ax2.scatter(
        tm,
        tg,
        facecolors="none",
        # facecolors="r",
        edgecolors="r",
        lw=0.1,
    )
    ax1.set_xlim([0, 2.5])
    ax1.set_ylim([0, 2.5])
    ax1.set_xlabel("mobbrmsd")
    ax1.set_ylabel("X-GRMSD")
    ax1.set_box_aspect(1)
    xlim = logymax(
        np.max(
            [
                np.max(tm),
                np.max(tg),
            ]
        )
    )
    xmin = (
        logymax(
            np.min(
                [
                    np.min(tm),
                    np.min(tg),
                ]
            )
        )
        * 0.1
    )
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim([xmin, xlim])
    ax2.set_ylim([xmin, xlim])
    ax2.set_xlabel("mobbrmsd")
    ax2.set_ylabel("X-GRMSD")
    ax2.set_box_aspect(1)
    plt.tight_layout()
    plt.savefig(f"vs_grmsd_{n_mol}.png")
    plt.savefig(f"vs_grmsd_{n_mol}.eps")
    plt.clf()
    plt.close()
