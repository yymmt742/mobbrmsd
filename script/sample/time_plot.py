import sys
import math
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import colorsets

cmap = colors.LinearSegmentedColormap.from_list(
    "mymap",
    [colorsets.ublk, colorsets.ublu, colorsets.sblu, colorsets.sylw, colorsets.uwht],
)


def xbound(x):
    xd = 10 ** math.floor(np.log10(x))
    return (math.floor(x / xd) + 1) * xd


for a in sys.argv[1:]:
    stem = Path(a).stem
    dat = np.load(a)
    print(stem)
    xl = np.mean(dat[:, 1]) - np.std(dat[:, 1])
    xu = np.mean(dat[:, 1]) + np.std(dat[:, 1])
    yl = np.mean(dat[:, 2]) - np.std(dat[:, 2])
    yu = np.mean(dat[:, 2]) + np.std(dat[:, 2])
    bins = (np.linspace(xl, xu, 101), np.linspace(yl, yu, 101))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    h, xedges, yedges, img = ax.hist2d(
        dat[:, 1],
        dat[:, 2],
        bins=bins,
        cmap=cmap,
    )
    cmax = xbound(np.max(h))
    ax.set_title(stem)
    ax.set_xlabel("Neval")
    ax.set_ylabel("time")
    img.set_clim(0, cmax)
    fig.colorbar(img, ax=ax)
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.savefig(stem + ".eps")
    plt.savefig(stem + ".png")
    plt.clf()
    plt.cla()
    plt.close()

    continue

    xm = xbound(np.mean(dat[:, 0]) + np.std(dat[:, 0]) * 2)
    ym = xbound(np.mean(dat[:, 1]) + np.std(dat[:, 1]) * 2)
    bins = (np.linspace(0.0, xm, 101), np.linspace(0, ym, min(ym, 100) + 1))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    h, xedges, yedges, img = ax.hist2d(
        dat[:, 0],
        dat[:, 1],
        bins=bins,
        cmap=cmap,
    )
    cmax = xbound(np.max(h))
    ax.set_title(stem)
    ax.set_xlabel("moRMSD")
    ax.set_ylabel("Neval")
    img.set_clim(0, cmax)
    fig.colorbar(img, ax=ax)
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.savefig(stem + "_hmap.eps")
    plt.savefig(stem + "_hmap.png")
    plt.clf()
    plt.cla()
    plt.close()
