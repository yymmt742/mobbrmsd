import sys
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
from pathlib import Path


def hex_to_RGB(hex_str):
    """#FFFFFF -> [255,255,255]"""
    # Pass 16 to the integer function for change of base
    return [int(hex_str[i : i + 2], 16) for i in range(1, 6, 2)]


def get_color_gradient(c1, c2, n):
    """
    Given two hex colors, returns a color gradient
    with n colors.
    """
    assert n > 1
    c1_rgb = np.array(hex_to_RGB(c1)) / 255
    c2_rgb = np.array(hex_to_RGB(c2)) / 255
    mix_pcts = [x / (n - 1) for x in range(n)]
    rgb_colors = [((1 - mix) * c1_rgb + (mix * c2_rgb)) for mix in mix_pcts]
    return [
        "#" + "".join([format(int(round(val * 255)), "02x") for val in item])
        for item in rgb_colors
    ]


uwht = "#FFFFFF"
ugrb = "#C8C8C8"
ugrd = "#7F878F"
ublk = "#000000"

ured = "#FF2800"
ublu = "#0041FF"
ugrn = "#35A16B"
uylw = "#FAF500"
usbu = "#66CCFF"
upnk = "#FF99A0"
uorg = "#FF9900"
uppl = "#9A0079"
ubrw = "#663300"

sred = "#FFD1D1"
sblu = "#B4EBFA"
sgrn = "#87E7B0"
sygr = "#CBF266"
sylw = "#FFFF99"
sorg = "#EDC58F"
sppl = "#C7B2DE"

mycmap1 = colors.LinearSegmentedColormap.from_list(
    "mycmap1", [ublk, ublu, sblu, sylw, uwht]
)


def logymax(v):
    return 10 ** np.floor(np.log10(v) + 1)


def plot_corr_dummy_map(ax, x, y, lim, bw=None):
    kde = sp.stats.gaussian_kde(
        np.log10(np.vstack([x, y])),
        bw_method=bw,
    )
    xs = np.linspace(np.log10(lim[0]), np.log10(lim[1]), 100)
    ys = np.linspace(np.log10(lim[0]), np.log10(lim[1]), 100)
    xx, yy = np.meshgrid(xs, ys)
    zz = np.vstack([xx.flatten(), yy.flatten()])
    z = kde.pdf(zz).reshape([ys.size, xs.size])
    level = np.max(z) * np.array([0.01, 0.1])
    cs = ax.contour(
        xs,
        ys,
        z,
        levels=level,
        colors=[sorg, ured],
    )


df = pd.concat([pd.read_pickle(a, compression="gzip") for a in sys.argv[1:]])

plt.rcParams["font.size"] = 8
fig = plt.figure()
ax = fig.add_subplot(111)

systems = np.unique(df["system"])[[0, 5, 2, 3, 6, 8, 7, 1, 4]]

for key, c, mec, mfc, delta, marker, ms, mss, mew, elinewidth in zip(
    ["time_rdkit", "time_mobbrmsd", "time_lprmsd_pg"],
    [uppl, uorg, ublk],
    [uppl, uorg, ublk],
    [uwht, uwht, uwht],
    [0.0, 0.0, 0.0],
    ["s", "o", "_"],
    [6, 6, 8],
    [2, 1, 4],
    [1.0, 1.0, 0.8],
    [2, 4, 2],
):
    for k, s in enumerate(systems):
        t = df[key][(df["system"] == s)]
        x = t[~(t.isnull())]
        if len(x) < 1:
            continue
        med = np.percentile(x, 50)
        n_points = x.shape[0]
        if n_points > 9:
            err = [[med - np.percentile(x, 25)], [np.percentile(x, 75) - med]]
            ax.errorbar(
                [k + delta],
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
                [k + delta] * n_points,
                x,
                facecolors=c,
                edgecolors="none",
                s=mss,
                lw=0.2,
            )
            ax.scatter(
                [k + delta],
                [med],
                edgecolors=c,
                facecolors=mfc,
                marker=marker,
                s=4 * mss,
            )

ax.set_xticks(range(len(systems)), systems, rotation=90)
ax.set_yscale("log")
ax.set_ylim([0.000001, 1000])
ax.set_yticks([0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000])
ax.set_ylabel("Wall-clock time (sec/run)")
ax.set_box_aspect(1.414)
ax.tick_params(which="minor", direction="out", axis="y", length=2, width=0.5)
plt.tight_layout()
plt.savefig(f"conformer_of_hydrogens_summary.png")
plt.clf()
plt.close()


for method, label in zip(
    ["rdkit", "lprmsd", "lprmsd_pg"], ["GetBestRMS", "LP-RMSD", "LP-RMSD-PG"]
):
    # for method, label in zip(["rdkit", "lprmsd"], ["GetBestRMS", "LP-RMSD"]):
    for system in systems:
        fig = plt.figure()
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        t = df["mobbrmsd"][(df["system"] == system)]
        xm = t[~(t.isnull())]
        t = df[method][(df["system"] == system)]
        xg = t[~(t.isnull())]
        if len(xg) < 1:
            continue

        n_points = xm.shape[0]
        lim = 4.0

        ax1.plot([0, lim], [0, lim], c=ugrd)
        ax1.scatter(
            xg[:100],
            xm[:100],
            facecolors="none",
            edgecolors="r",
            lw=0.1,
        )
        ax1.set_xlim([0, lim])
        ax1.set_ylim([0, lim])
        ax1.set_xlabel(label)
        ax1.set_ylabel("mobbrmsd")
        ax1.set_xticks([0, 2, 4])
        ax1.set_yticks([0, 2, 4])
        ax1.set_box_aspect(1)

        t = df["time_mobbrmsd"][(df["system"] == system)]
        tm = t[~(t.isnull())]
        t = df["time_" + method][(df["system"] == system)]
        tg = t[~(t.isnull())]

        xmin = 0.000001
        xlim = 1000.0
        if n_points > 100:
            plot_corr_dummy_map(
                ax2,
                tg,
                tm,
                [xmin, xlim],
            )
            ax2.set_xticks([-6, -3, 0, 3])
            ax2.set_yticks([-6, -3, 0, 3])
        else:
            ax2.scatter(
                tm,
                tg,
                facecolors="none",
                edgecolors="r",
                lw=0.1,
            )
            ax2.set_xscale("log")
            ax2.set_yscale("log")
            ax2.set_xlim([xmin, xlim])
            ax2.set_ylim([xmin, xlim])

            ax2.set_xticks([0.000001, 0.001, 1.0, 1000.0])
            ax2.set_yticks([0.000001, 0.001, 1.0, 1000.0])

        ax2.set_xlabel(label)
        ax2.set_ylabel("mobbrmsd")
        ax2.set_box_aspect(1)
        plt.savefig(f"{method}_{system}.png")
        plt.clf()
        plt.close()
