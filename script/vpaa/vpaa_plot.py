#! /usr/bin/env python3
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import seaborn as sns
import pandas as pd
import scipy as sp

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

cmap = colors.LinearSegmentedColormap.from_list(
    "mymap",
    [ublk, ublu, sblu, sylw, uwht],
)
cmap.set_bad(ublk)


def logymax(v):
    return 10 ** (np.floor(np.log10(v)) + 1)


def hex_to_RGB(hex_str):
    return [int(hex_str[i : i + 2], 16) for i in range(1, 6, 2)]


def get_color_gradient(c1, c2, n):
    assert n > 1
    c1_rgb = np.array(hex_to_RGB(c1)) / 255
    c2_rgb = np.array(hex_to_RGB(c2)) / 255
    mix_pcts = [x / (n - 1) for x in range(n)]
    rgb_colors = [((1 - mix) * c1_rgb + (mix * c2_rgb)) for mix in mix_pcts]
    return [
        "#" + "".join([format(int(round(val * 255)), "02x") for val in item])
        for item in rgb_colors
    ]


def nstats(m):
    return np.exp(np.log(2.0) * m + np.sum([np.log(i) for i in range(m, 0, -1)]))


def nnodes(m):
    ret = 0
    c = 1
    for i in range(m, 0, -1):
        c *= i * 2
        ret += c
    return ret


def plot_violins(ax, data, y):
    ind = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10])
    mask = (
        (data["coff"] == 1000.0)
        & (data["uoff"] == 1000.0)
        & ((data["sys"] == "vpaa1") | (data["sys"] == "vpaa2"))
    )
    hue = data[mask]["sys"]
    quarts = np.array(
        [
            [
                np.percentile(
                    (data[mask & (data["n_mol"] == n_mol) & (data["sys"] == s)])[y],
                    [50, 25, 75],
                )
                for n_mol in ind
            ]
            for s in ["vpaa1", "vpaa2"]
        ]
    )
    xdata = data[mask]["n_mol"]
    ydata = data[mask][y]
    sns.violinplot(
        x=xdata,
        y=ydata,
        palette={"vpaa1": ublu, "vpaa2": uorg},
        saturation=1.0,
        width=1.4,
        hue=hue,
        density_norm="width",
        split=True,
        orient="v",
        gap=0.45,
        log_scale=True,
        inner=None,
        bw_adjust=0.1,
        fill=True,
        ax=ax,
    )
    med0 = quarts[0, :, 0]
    err0 = np.array([med0 - quarts[0, :, 1], quarts[0, :, 2] - med0])
    med1 = quarts[1, :, 0]
    err1 = np.array([med1 - quarts[1, :, 1], quarts[1, :, 2] - med1])
    lb = (ind + 1) * ind
    brute = np.array([nstats(m) for m in ind])
    ub = np.array([nnodes(m) for m in ind])
    errb = [brute - lb, ub - brute]
    ax.errorbar(
        ind - 2.0,
        brute,
        color="gray",
        yerr=errb,
        marker="_",
        elinewidth=1,
        ms=8,
        linestyle="",
        capsize=4,
    )
    ax.errorbar(
        ind - 2.07,
        med0,
        color=ublu,
        yerr=err0,
        marker=".",
        mfc="white",
        mew=0,
        elinewidth=4,
        ms=5,
        linestyle="",
        capsize=0,
    )
    ax.errorbar(
        ind - 1.93,
        med1,
        color=uorg,
        yerr=err1,
        marker=".",
        mfc="white",
        mew=0,
        elinewidth=4,
        ms=5,
        linestyle="",
        capsize=0,
    )
    xmin = 0
    xmax = 8
    ymin = logymax(np.min(ydata)) / 10
    ymax = logymax(brute[-1])
    ax.set_xlim(
        (
            xmin - 0.7,
            xmax + 0.7,
        )
    )
    ax.set_ylim(
        (
            ymin,
            ymax,
        )
    )
    ax.set_box_aspect(0.707)
    ax.legend()


def plot_corr_dummy_map(ax, data, n_mol, bw=0.05, xlim=1.2, ylim=1.2):
    mask = (data["n_mol"] == n_mol) & (data["sys"] == "dummy")
    kde = sp.stats.gaussian_kde(
        np.vstack([(data[mask])["rmsd"], (data[mask])["n_eval"]]),
        bw_method=bw,
    )
    x = np.linspace(0.0, np.max((data[mask])["rmsd"]) * xlim, 100)
    y = np.linspace(0.0, np.max((data[mask])["n_eval"]) * ylim + 10, 100)
    xx, yy = np.meshgrid(x, y)
    zz = np.vstack([xx.flatten(), yy.flatten()])
    z = kde.pdf(zz).reshape([y.size, x.size])
    level = np.max(z) * np.array([0.01, 0.1])
    cs = ax.contour(
        x,
        y,
        z,
        levels=level,
        colors="#555555",
    )
    plt.clabel(
        cs,
        inline_spacing=2.0,
        fmt={
            level[0]: "1%",
            level[1]: "10%",
        },
    )


def plot_corr_map(ax, data, n_mol, coff, bw=0.05, xmax=None, ymax=None):

    print(n_mol, coff)
    ax.set_box_aspect(1.414)
    mask = (data["n_mol"] == n_mol) & (data["coff"] == coff) & (data["uoff"] == 1000.0)
    mask1 = mask & (data["sys"] == "vpaa1")
    mask2 = mask & (data["sys"] == "vpaa2")
    kde1 = sp.stats.gaussian_kde(
        np.vstack([(data[mask1])["rmsd"], (data[mask1])["n_eval"]]),
        bw_method=bw,
    )
    kde2 = sp.stats.gaussian_kde(
        np.vstack([(data[mask2])["rmsd"], (data[mask2])["n_eval"]]),
        bw_method=bw,
    )
    if xmax is None:
        xmax = np.max(data["rmsd"])
    if ymax is None:
        ymax = logymax(np.max(data["n_eval"]))
    x = np.linspace(0.0, xmax, 100)
    y = np.linspace(0.0, ymax, 100)
    xx, yy = np.meshgrid(x, y)
    zz = np.vstack([xx.flatten(), yy.flatten()])
    z1 = kde1.pdf(zz).reshape([y.size, x.size])
    z2 = kde2.pdf(zz).reshape([y.size, x.size])
    level1 = np.max(z1) * np.array([0.01, 0.1])
    level2 = np.max(z2) * np.array([0.01, 0.1])
    cs1 = ax.contour(
        x,
        y,
        z1,
        levels=level1,
        linestyles=["dotted", "solid"],
        colors=ublu,
    )
    cs2 = ax.contour(
        x,
        y,
        z2,
        levels=level2,
        linestyles=["dotted", "solid"],
        colors=uorg,
    )
    """
    plt.clabel(
        cs1,
        inline_spacing=2.0,
        fmt={
            level1[0]: "1%",
            level1[1]: "10%",
        },
    )
    plt.clabel(
        cs2,
        inline_spacing=2.0,
        fmt={
            level2[0]: "1%",
            level2[1]: "10%",
        },
    )
    """


def plot_corr_lines(ax, data, n_mol, coff, ymin=None, ymax=None):
    perf = 100 - np.arange(50) * 2
    print(n_mol, coff)
    mask = (data["n_mol"] == n_mol) & (data["coff"] == coff) & (data["uoff"] == 1000.0)
    hue = ["vpaa1", "vpaa2"]
    rmsds = np.array(
        [
            [
                np.percentile(
                    (data[mask & (data["sys"] == s)])["rmsd"],
                    p,
                )
                for s in hue
            ]
            for p in perf
        ]
    )
    nevals = np.array(
        [
            [
                np.percentile(
                    (data[mask & (data["sys"] == s) & (data["rmsd"] < conds)])[
                        "n_eval"
                    ],
                    [50, 25, 75],
                )
                for s, conds in zip(hue, cond)
            ]
            for cond in rmsds
        ]
    )
    colors = [ublu, uorg, ugrb]
    for i in range(nevals.shape[1]):
        """
        ax.fill_between(
            perf / 100,
            nevals[:, i, 1],
            nevals[:, i, 2],
            color=colors[i],
            linewidth=0,
            alpha=0.2,
        )
        """
        ax.plot(
            perf / 100,
            nevals[:, i, 1],
            color=colors[i],
            linewidth=1,
            linestyle="dotted",
        )
        ax.plot(
            perf / 100,
            nevals[:, i, 2],
            color=colors[i],
            linewidth=1,
            linestyle="dotted",
        )
        ax.plot(
            perf / 100,
            nevals[:, i, 0],
            color=colors[i],
            linewidth=4,
            # label=hue[i] + f"-{n_mol:d}",
        )
    ax.set_yscale("log")
    xmin = 0
    xmax = 1
    if ymin is None:
        ymin = logymax(np.min(data["n_eval"])) / 10
    if ymax is None:
        ymax = logymax(np.max(data["n_eval"]))
    ax.set_xlim(
        (
            xmin,
            xmax,
        )
    )
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_ylim(
        (
            ymin,
            ymax,
        )
    )
    ax.set_box_aspect(1.414)
    # ax.legend(loc="upper left")


def plot_bar_lines(ax, stats, y, scale=1.0):
    ind = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10])
    hue = ["vpaa1", "vpaa2", "dummy"]
    mask = (stats["coff"] == 1000.0) & (stats["uoff"] == 1000.0)

    cs = [ublu, uorg, ugrb]
    lss = ["-", ":", "-."]
    lws = [2, 3, 2]
    for vsys, c, ls, lw in zip(hue, cs, lss, lws):
        ax.plot(
            ind,
            scale * np.array(stats[mask & (stats["sys"] == vsys)][y + "_q50"]),
            color=c,
            linestyle=ls,
            lw=lw,
        )
    ymax = 0.0
    ymin = 9999.9
    for dx, vsys, c, ls, lw in zip([-0.15, 0.15, 0.0], hue, cs, lss, lws):
        med = scale * np.array(stats[mask & (stats["sys"] == vsys)][y + "_q50"])
        err = [
            med - scale * np.array(stats[mask & (stats["sys"] == vsys)][y + "_q25"]),
            scale * np.array(stats[mask & (stats["sys"] == vsys)][y + "_q75"]) - med,
        ]
        ymax = max(ymax, np.max(med + err[1]))
        ymin = min(ymin, np.min(med - err[0]))
        ax.errorbar(
            ind + dx,
            med,
            color=c,
            yerr=err,
            marker=".",
            mfc="white",
            mew=0,
            elinewidth=5,
            ms=6,
            linestyle="",
            capsize=0,
        )
    xmin = np.min(ind[0])
    xmax = np.max(ind[-1])
    ymin = logymax(ymin) / 10
    ymax = logymax(ymax)
    ax.set_yscale("log")
    ax.set_xticks(ind)
    ax.set_xlim(
        (
            xmin - 0.5,
            xmax + 0.5,
        )
    )
    ax.set_ylim(
        (
            ymin,
            ymax,
        )
    )
    ax.set_box_aspect(2)
    # ax.legend(loc="upper left")


def plot_cputime(stats):
    print(f"plot cputimes")
    fig, ax1 = plt.subplots(figsize=(5.0, 5.0 * 1.414), dpi=100)
    plot_bar_lines(ax1, stats, "t")
    # plt.show()
    fig.savefig("cputime.png")
    fig.savefig("cputime.eps")
    plt.clf()
    plt.close()


def plot_ub(data):
    print(f"plot ubs")
    vmin = 1.0
    vmax = 10000.0
    fig = plt.figure(figsize=(2.8 * 4.0, 4.0), dpi=200)
    for i, vsys in enumerate(["vpaa1", "vpaa2"]):
        for j, n_mol in enumerate([2, 3, 4, 5, 6, 7, 8, 9, 10]):
            ax = fig.add_subplot(
                2,
                10,
                1 + 10 * i + j,
                xticks=[],
                yticks=[],
            )
            m = (data["sys"] == vsys) & (data["n_mol"] == n_mol) & (data["coff"] == 0.0)
            h, x, y = np.histogram2d(
                data[m]["ulim"],
                data[m]["llim"],
                bins=[70, 40],
                range=[[0.0, 1.4], [0.0, 0.8]],
            )
            im = ax.imshow(
                h[::-1],
                cmap=cmap,
                norm=LogNorm(),
                interpolation="none",
            )
            ax.set_box_aspect(14 / 8)
            im.set_clim(vmin, vmax)
    ax = fig.add_subplot(1, 10, 10, xticks=[], yticks=[])
    ax.set_frame_on(False)
    im = ax.imshow(
        [[vmax]],
        cmap=cmap,
        norm=LogNorm(),
        interpolation="none",
    )
    ax.set_box_aspect(50)
    im.set_clim(vmin, vmax)
    fig.colorbar(im)
    plt.tight_layout()
    fig.savefig("ub.png")
    fig.savefig("ub.eps")
    # plt.show()
    plt.clf()
    plt.close()


def plot_gap(data):
    print(f"plot gaps")
    vmin = 1.0
    vmax = 10000.0
    fig1 = plt.figure(figsize=(4.0, 2 * 4.0), dpi=200)
    fig2 = plt.figure(figsize=(4.0, 2 * 4.0), dpi=200)
    llist = np.unique(data["coff"])[::-1][1:-1]
    ulist = np.unique(data["uoff"])[::-1][2:]
    print(llist)
    print(ulist)
    exit()
    for i, vsys in enumerate(["vpaa1", "vpaa2"]):
        for j, n_mol in enumerate([2, 3, 4, 5, 6, 7, 8, 9, 10]):
            for fig, vcut, rcut, clist in zip(
                [fig1, fig2], ["coff", "uoff"], ["uoff", "coff"], [llist, ulist]
            ):
                ax = fig.add_subplot(
                    9,
                    2,
                    1 + 2 * j + i,
                    # yticks=[0, 1, 2],
                    # yticks=[0, 2, 4, 6],
                    yticks=[-4, -2, 0, 2, 4],
                )
                m = (
                    (data["sys"] == vsys)
                    & (data["n_mol"] == n_mol)
                    & (data[rcut] == 1000.0)
                )
                m0 = m & (data[vcut] == 1000.0)

                for coff, c, lw in zip(
                    clist,
                    get_color_gradient(uylw, uppl, len(clist)),
                    np.linspace(3.0, 1.0, len(clist)),
                ):

                    mm = m & (data[vcut] == coff)
                    h, x, y = np.histogram2d(
                        data[m0]["rmsd"],
                        np.array(data[mm]["ulim"]) - np.array(data[m0]["rmsd"]),
                        bins=[80, 2400],
                        range=[[0.0, 0.8], [0.0, 2.4]],
                    )
                    xx = 10 * 0.5 * (x[1:] + x[:-1])
                    yy = 10 * 0.5 * (y[1:] + y[:-1])
                    h0 = np.sum(h, 1)
                    gap = np.sum(yy * h, 1) / np.where(h0 < 0.5, 1.0, h0)
                    ax.plot(
                        xx,
                        gap,
                        color=c,
                        lw=lw,
                        label=f"{coff}",
                    )
                    h, x, y = np.histogram2d(
                        data[m0]["rmsd"],
                        np.array(data[mm]["llim"]) - np.array(data[m0]["rmsd"]),
                        bins=[80, 2400],
                        range=[[0.0, 0.8], [-2.4, 0.0]],
                    )
                    xx = 10 * 0.5 * (x[1:] + x[:-1])
                    yy = 10 * 0.5 * (y[1:] + y[:-1])
                    h0 = np.sum(h, 1)
                    gap = np.sum(yy * h, 1) / np.where(h0 < 0.5, 1.0, h0)
                    ax.plot(
                        xx,
                        gap,
                        color=c,
                        lw=lw,
                        label=f"{coff}",
                    )
                ax.set_xlim([0, 8])
                # ax.set_ylim([0, 2])
                ax.set_ylim([-4, 4])
                ax.set_xticks(
                    [0, 1, 2, 3, 4, 5, 6, 7, 8],
                    labels=["", "", "", "", "", "", "", "", ""],
                )
                ax.set_box_aspect(1 / 1.618)
                # ax.legend()
    fig1.tight_layout()
    fig2.tight_layout()
    fig1.savefig("gap_coff.png")
    fig1.savefig("gap_coff.eps")
    fig2.savefig("gap_uoff.png")
    fig2.savefig("gap_uoff.eps")
    # plt.show()
    plt.clf()
    plt.close()


def plot_violin(data):
    print(f"plot violins")
    fig, ax = plt.subplots()
    plot_violins(ax, data, "n_eval")
    fig.savefig("violin.png")
    fig.savefig("violin.eps")
    # plt.show()
    plt.clf()
    plt.close()


def plot_corr_map_eps(data):
    coffs = np.unique(data[(data["sys"] == "vpaa1") | (data["sys"] == "vpaa2")]["coff"])
    n = len(coffs)
    fig, axes = plt.subplots(n, 2, figsize=(4.0, n * 4.0 * 1.618), dpi=100)
    for coff, ax in zip(coffs[::-1], axes):
        plot_corr_map(ax[0], data, 10, coff, ymax=10000000)
        plot_corr_lines(ax[1], data, 10, coff)
    plt.tight_layout()
    # plt.show()
    fig.savefig("corr_map.png")
    fig.savefig("corr_map.eps")
    plt.clf()
    plt.close()


def plot_corr_maps_eps(data):
    print(f"plot corr_maps")
    fig = plt.figure(figsize=(8.0, 6.0 * 1.414), dpi=100)
    ax2 = plt.subplot2grid((3, 3), (0, 0))
    ax3 = plt.subplot2grid((3, 3), (0, 1))
    ax4 = plt.subplot2grid((3, 3), (0, 2))
    ax5 = plt.subplot2grid((3, 3), (1, 0))
    ax6 = plt.subplot2grid((3, 3), (1, 1))
    ax7 = plt.subplot2grid((3, 3), (1, 2))
    ax8 = plt.subplot2grid((3, 3), (2, 0))
    ax9 = plt.subplot2grid((3, 3), (2, 1))
    ax0 = plt.subplot2grid((3, 3), (2, 2))
    plot_corr_map(ax2, data, 2, 1000.0, xmax=1.0, ymax=16.0)
    plot_corr_map(ax3, data, 3, 1000.0, xmax=1.0, ymax=80.0)
    plot_corr_map(ax4, data, 4, 1000.0, xmax=1.0, ymax=400.0)
    plot_corr_map(ax5, data, 5, 1000.0, xmax=1.0, ymax=1600.0)
    plot_corr_map(ax6, data, 6, 1000.0, xmax=1.0, ymax=10000.0)
    plot_corr_map(ax7, data, 7, 1000.0, xmax=1.0, ymax=40000.0)
    plot_corr_map(ax8, data, 8, 1000.0, xmax=1.0, ymax=250000.0)
    plot_corr_map(ax9, data, 9, 1000.0, xmax=1.0, ymax=1000000.0)
    plot_corr_map(ax0, data, 10, 1000.0, xmax=1.0, ymax=6000000.0)

    plt.tight_layout()
    plt.savefig("corr_map.png")
    plt.savefig("corr_map.eps")
    plt.clf()
    plt.close()


def plot_corr_lines_eps(data):
    print(f"plot corr_lines")
    fig = plt.figure(figsize=(8.0, 6.0 * 1.414), dpi=100)
    ax2 = plt.subplot2grid((3, 3), (0, 0))
    ax3 = plt.subplot2grid((3, 3), (0, 1))
    ax4 = plt.subplot2grid((3, 3), (0, 2))
    ax5 = plt.subplot2grid((3, 3), (1, 0))
    ax6 = plt.subplot2grid((3, 3), (1, 1))
    ax7 = plt.subplot2grid((3, 3), (1, 2))
    ax8 = plt.subplot2grid((3, 3), (2, 0))
    ax9 = plt.subplot2grid((3, 3), (2, 1))
    ax0 = plt.subplot2grid((3, 3), (2, 2))

    plot_corr_lines(ax2, data, 2, 1000.0, ymin=5.0, ymax=11.0)
    plot_corr_lines(ax3, data, 3, 1000.0, ymin=10.0, ymax=100.0)
    plot_corr_lines(ax4, data, 4, 1000.0, ymin=10.0, ymax=1000.0)
    plot_corr_lines(ax5, data, 5, 1000.0, ymin=10.0, ymax=1000.0)
    plot_corr_lines(ax6, data, 6, 1000.0, ymin=10.0, ymax=10000.0)
    plot_corr_lines(ax7, data, 7, 1000.0, ymin=100.0, ymax=100000.0)
    plot_corr_lines(ax8, data, 8, 1000.0, ymin=100.0, ymax=100000.0)
    plot_corr_lines(ax9, data, 9, 1000.0, ymin=100.0, ymax=1000000.0)
    plot_corr_lines(ax0, data, 10, 1000.0, ymin=100.0, ymax=2000000.0)

    plt.tight_layout()
    plt.savefig("corr_lines.png")
    plt.savefig("corr_lines.eps")
    plt.clf()


def get_stats(data):
    stat = []
    for vsys in pd.unique(data["sys"]):
        m1 = data["sys"] == vsys
        slist = lambda vsys, n_mol, coff, uoff, data, mask: [
            [vsys, n_mol, coff, uoff]
            + [np.mean(data[mask]["n_eval"])]
            + [np.var(data[mask]["n_eval"])]
            + np.percentile(data[mask]["n_eval"], [25, 50, 75]).tolist()
            + [np.mean(0.001 * data[mask]["times"])]
            + [np.var(0.001 * data[mask]["times"])]
            + np.percentile(0.001 * data[mask]["times"], [25, 50, 75]).tolist()
        ]
        for n_mol in pd.unique(data["n_mol"]):
            print(f"get stats for {vsys} n_mol={n_mol:2d}")
            m2 = m1 & (data["n_mol"] == n_mol)
            for coff in pd.unique((data[m2])["coff"]):
                mask = m2 & (data["coff"] == coff) & (data["uoff"] == 1000.0)
                stat += slist(vsys, n_mol, coff, 1000.0, data, mask)
            for uoff in pd.unique((data[m2])["uoff"]):
                if uoff == 1000.0:
                    continue
                mask = m2 & (data["coff"] == 1000.0) & (data["uoff"] == uoff)
                stat += slist(vsys, n_mol, 1000.0, uoff, data, mask)

    return pd.DataFrame(
        stat,
        columns=[
            "sys",
            "n_mol",
            "coff",
            "uoff",
            "n_mean",
            "n_var",
            "n_q25",
            "n_q50",
            "n_q75",
            "t_mean",
            "t_var",
            "t_q25",
            "t_q50",
            "t_q75",
        ],
    )
    return


def plot_cutoff_score_lines(ax, data, vsys, key, y, vcut, clist):

    ind = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10])
    mask = (
        lambda data, vcut, cut, vsys: (data[vcut] == cut)
        & (data[("coff" if vcut == "uoff" else "uoff")] == 1000.0)
        & (data["sys"] == vsys)
    )
    ax.plot(
        ind,
        data[mask(data, vcut, 1000.0, vsys)][key + "_mean"],
        color=ublk,
        lw=1.0,
        label="full",
    )
    for cut, c, lw in zip(
        clist,
        get_color_gradient(uylw, uppl, len(clist)),
        np.linspace(3.0, 1.0, len(clist)),
    ):
        ms = mask(data, vcut, cut, vsys)
        mean = data[ms][key + "_mean"]
        label = 10 * cut
        ax.plot(
            ind,
            mean,
            color=c,
            lw=lw,
            label=f"{label}",
        )
    ax.plot(
        ind,
        data[mask(data, "coff", 0.0, vsys)][key + "_mean"],
        color=ublk,
        lw=1.0,
        label="greedy",
    )
    xmin = ind[0]
    xmax = ind[-1]
    ymin = logymax(np.min(data[mask(data, vcut, 1000.0, vsys)][key + "_mean"])) / 10
    ymax = logymax(np.max(data[mask(data, vcut, 1000.0, vsys)][key + "_mean"]))
    xrange = xmax - xmin + 1.0
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(ind)
    ax.set_xticks(ind, [str(i) for i in ind])
    ax.set_xlim(
        (
            xmin / 1.1,
            xmax * 1.1,
        )
    )
    ax.set_xticks(ind, [str(i) for i in ind])
    ax.set_ylim(
        (
            ymin,
            ymax,
        )
    )
    ax.set_box_aspect(1.618)
    ax.legend()


def plot_cutoff_score(data, vcut, clist):
    print(f"plot cutoff_score {vcut}")
    fig, axes = plt.subplots(2, 2, figsize=[5.0, 1.414 * 5])
    for vsys, ax in zip(["vpaa1", "vpaa2"], axes[0]):
        plot_cutoff_score_lines(ax, data, vsys, "n", "n_eval", vcut, clist)
    for vsys, ax in zip(["vpaa1", "vpaa2"], axes[1]):
        plot_cutoff_score_lines(ax, data, vsys, "t", "times", vcut, clist)
    plt.tight_layout()
    plt.savefig(f"cutoff_{vcut}.eps")
    plt.savefig(f"cutoff_{vcut}.png")
    # plt.show()
    plt.clf()
    plt.close()


def plot_time_lines(ax, data, vsys, key, y, vcut):
    ind = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10])
    mask = (
        lambda data, vcut, cut, vsys: (data[vcut] == cut)
        & (data[("coff" if vcut == "uoff" else "uoff")] == 1000.0)
        & (data["sys"] == vsys)
    )
    for cut, c, lw in zip(
        np.unique(data[vcut])[::-1],
        get_color_gradient(uylw, uppl, 6),
        [3.0, 2.5, 2.0, 1.75, 1.5, 1.0],
    ):
        ms = mask(data, coff, 1000.0, vsys)
        mean = data[ms][key + "_mean"]
        ax.plot(
            ind,
            mean,
            color=c,
            lw=lw,
            label=f"{cut}",
        )
    xmin = ind[0]
    xmax = ind[-1]
    ymin = logymax(np.min(data[key + "_mean"])) / 10
    ymax = logymax(np.max(data[key + "_mean"]))
    xrange = xmax - xmin + 1.0
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(ind)
    ax.set_xticks(ind, [str(i) for i in ind])
    ax.set_xlim(
        (
            xmin / 1.1,
            xmax * 1.1,
        )
    )
    ax.set_xticks(ind, [str(i) for i in ind])
    ax.set_ylim(
        (
            ymin,
            ymax,
        )
    )
    ax.set_box_aspect(1.618)
    ax.legend()


def plot_time_eps(data, vcut):
    fig, axes = plt.subplots(2, 2, figsize=[5.0, 1.414 * 5])
    for vsys, ax in zip(["vpaa1", "vpaa2"], axes[0]):
        plot_time_lines(ax, data, vsys, "n", "n_eval", vcut)
    for vsys, ax in zip(["vpaa1", "vpaa2"], axes[1]):
        plot_time_lines(ax, data, vsys, "t", "times", vcut)
    plt.tight_layout()
    plt.savefig(f"cutoff_{vcut}.eps")
    plt.savefig(f"cutoff_{vcut}.png")
    plt.show()
    plt.clf()


data = pd.concat([pd.read_pickle(f, compression="gzip") for f in sys.argv[1:]])
print(data)

# plot_gap(data)
# plot_ub(data)
plot_violin(data)
plot_corr_maps_eps(data)
plot_corr_map_eps(data)
plot_corr_lines_eps(data)
stats = get_stats(data)
print(stats)

plot_cputime(stats)
plot_cutoff_score(stats, "coff", [0.3, 0.2, 0.15, 0.1])
plot_time_eps(data, "coff")

# plot_cutoff_score(stats, "uoff", [0.8, 0.6, 0.5, 0.4])
# plot_time_eps(data, "uoff")
