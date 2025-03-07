import sys
import numpy as np
import scipy as sp
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


plt.rcParams["font.size"] = 8

df = pd.concat([pd.read_pickle(a, compression="gzip") for a in sys.argv[1:]])

fig = plt.figure()
ax = fig.add_subplot(111)
n_mols = np.unique(df["n_mol"])

ax.plot([5, 5], [0.000001, 1000000], ls=":", lw=0.5, c=ugrb)
ax.plot([10, 10], [0.000001, 1000000], ls=":", lw=0.5, c=ugrb)
ax.plot([15, 15], [0.000001, 1000000], ls=":", lw=0.5, c=ugrb)
ax.plot([0, max(n_mols) + 1], [0.001, 0.001], ls=":", lw=0.5, c=ugrb)
ax.plot([0, max(n_mols) + 1], [1.0, 1.0], ls=":", lw=0.5, c=ugrb)
ax.plot([0, max(n_mols) + 1], [1000.0, 1000.0], ls=":", lw=0.5, c=ugrb)

for system, c, mec, mfc, delta, marker, ms, mss, mew, elinewidth, pline in zip(
    ["grmsd", "frmsd", "mobbrmsd", "rdkit", "lprmsd"],
    [uppl, ublu, uorg, ubrw, ublk],
    [None, None, None, ubrw, ublk],
    [uwht, uwht, uwht, ubrw, ublk],
    [0.0, 0.0, 0.0, 0.0, 0.0],
    ["s", "v", "o", "x", "_"],
    [4, 6, 6, 4, 6],
    [2, 4, 4, 2, 3],
    [1.0, 1.0, 1.0, 1.0, 0.8],
    [1, 2, 3, 1, 1],
    [True, True, True, True, False],
):
    if pline:
        nm = []
        meds = []
        for n_mol in n_mols:
            t = df[f"time_{system}"][(df["n_mol"] == n_mol)]
            if len(t) == 0:
                continue
            x = t[~(t.isnull())]
            try:
                meds += [np.percentile(x, 50)]
            except:
                continue
            nm += [n_mol]
        ax.plot(nm, meds, c=c, ls="-", lw=0.5)

    for n_mol in n_mols:
        t = df[f"time_{system}"][(df["n_mol"] == n_mol)]
        if len(t) == 0:
            continue
        x = t[~(t.isnull())]
        try:
            med = np.percentile(x, 50)
        except:
            continue
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
                label=system if n_mol == 2 else None,
            )

        else:
            ax.scatter(
                [n_mol + delta] * n_points,
                x,
                facecolors=c,
                edgecolors="none",
                s=2,
                lw=0.2,
            )
            ax.scatter(
                [n_mol + delta],
                [med],
                edgecolors=c,
                facecolors=mfc,
                marker=marker,
                s=8 * mss,
                label=system if n_mol == 2 else None,
            )

ax.set_xlim([min(n_mols) - 0.5, max(n_mols) + 0.5])
ax.set_xticks(n_mols)
ax.set_yscale("log")
ax.set_ylim([0.000001, 1000000])
ax.set_yticks(
    [
        0.000001,
        0.00001,
        0.0001,
        0.001,
        0.01,
        0.1,
        1,
        10,
        100,
        1000,
        10000,
        100000,
        1000000,
    ]
)
ax.legend()
ax.set_box_aspect(1.414)
ax.tick_params(which="minor", direction="out", axis="y", length=2, width=0.5)
plt.savefig(f"random_coordinate_time.png")
plt.clf()
plt.close()
