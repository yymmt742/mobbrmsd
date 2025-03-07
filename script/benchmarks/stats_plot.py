import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import colorsets


def ymax(v):
    c = 10 ** np.floor(np.log10(v)) // 2
    return ((v // c) + 1) * c


def logymax(v):
    return 10 ** (np.floor(np.log10(v)) + 1)


n_mols = [3, 4, 5, 6, 7, 8, 9, 10]

df = pd.concat([pd.read_pickle(a, compression="gzip") for a in sys.argv[1:]])
stems = [Path(a).stem[:-4] for a in sys.argv[1:]]


def plot_curve():
    cs = (
        colorsets.get_color_gradient(colorsets.uylw, colorsets.sorg, 4)[1:]
        + colorsets.get_color_gradient(colorsets.sorg, colorsets.uppl, 4)[1:]
    )
    lws = [3.0, 2.5, 2.0, 1.5, 1.25, 1.0]
    for n_mol, stem in zip(n_mols, stems):

        ym = ymax(
            1.05 * np.max(df["n_mean"][(1 == df["n_sym"]) & (n_mol == df["n_mol"])])
        )
        ualph = np.unique(df["alpha"])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for na, c, lw in zip(np.unique(df["n_apm"]), cs, lws):
            plt.plot(
                ualph,
                df["n_mean"][
                    (na == df["n_apm"]) & (1 == df["n_sym"]) & (n_mol == df["n_mol"])
                ],
                label=str(na),
                c=c,
                lw=lw,
            )
        plt.legend()
        ax.set_title(stem)
        ax.set_xlabel("alpha")
        ax.set_ylabel("Neval")
        ax.set_xlim([-0.05, 0.95])
        ax.set_ylim([0, ym])
        # plt.tight_layout()
        ax.set_box_aspect(1.618)
        plt.savefig(stem + ".eps")
        plt.savefig(stem + ".png")
        plt.clf()
        plt.cla()
        plt.close()


def plot_m_n():
    ym = logymax(np.max(df["n_q75"][(1 == df["n_sym"])]) * 2)

    for n_apm in [1, 2, 3, 10, 100]:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for alpha, c, mfc, delta, marker, ms, mew, elinewidth in zip(
            [0.0, 0.9],
            [colorsets.ubrw, colorsets.uorg],
            [colorsets.uwht, colorsets.uwht],
            [-0.12, 0.18],
            ["_", "o"],
            [4, 4],
            [1.0, 0],
            [4, 8],
        ):
            mask = (n_apm == df["n_apm"]) & (alpha == df["alpha"]) & (1 == df["n_sym"])
            n_mol = df["n_mol"][mask].to_numpy()
            q25 = df["n_q25"][mask].to_numpy()
            q50 = df["n_q50"][mask].to_numpy()
            q75 = df["n_q75"][mask].to_numpy()

            eb = [q50 - q25, q75 - q50]
            # ax.plot(n_mol, q50, color=c, lw=1)
            ax.errorbar(
                n_mol + delta,
                q50,
                color=c,
                yerr=eb,
                linestyle="",
                capsize=0,
                label=f"{alpha}",
                marker=marker,
                mfc=mfc,
                mew=mew,
                mec="white",
                elinewidth=elinewidth,
                ms=ms,
            )
        ax.set_title(f"m-neval {n_apm:03d}")
        ax.set_xlabel("M")
        ax.set_ylabel("Neval")
        ax.set_xlim([2.5, 10.5])
        ax.set_xticks(
            [i + 3 for i in range(8)],
        )
        ax.set_ylim([1, ym])
        ax.set_yscale("log")
        plt.legend()
        # plt.tight_layout()
        ax.set_box_aspect(1.618)
        plt.savefig(f"m_neval_{n_apm:03d}.eps")
        plt.savefig(f"m_neval_{n_apm:03d}.png")
        plt.clf()
        plt.cla()
        plt.close()


def plot_s_n():
    for n_mol in [3, 4, 5, 6]:
        ym = logymax(np.max(df["n_q75"][n_mol == df["n_mol"]]) * 2)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for alpha, c, mfc, delta, marker, ms, mew, elinewidth in zip(
            [0.0, 0.9],
            [colorsets.ubrw, colorsets.uorg],
            [colorsets.uwht, colorsets.uwht],
            [-0.006, 0.012],
            ["_", "o", "D"],
            [4, 4],
            [1.0, 0],
            [4, 8],
        ):
            mask = (n_mol == df["n_mol"]) & (10 == df["n_apm"]) & (alpha == df["alpha"])
            log_n_sym = np.log10(df["n_sym"][mask])
            q25 = df["n_q25"][mask].to_numpy()
            q50 = df["n_q50"][mask].to_numpy()
            q75 = df["n_q75"][mask].to_numpy()

            eb = [q50 - q25, q75 - q50]
            ax.errorbar(
                log_n_sym + delta,
                q50,
                color=c,
                yerr=eb,
                linestyle="",
                capsize=0,
                label=f"{alpha}",
                marker=marker,
                mfc=mfc,
                mew=mew,
                mec="white",
                elinewidth=elinewidth,
                ms=ms,
            )
        ax.set_title(f"m-neval {n_mol:03d}")
        ax.set_xlabel("s")
        ax.set_ylabel("Neval")
        ax.set_xlim([-0.03, 1.03])
        ax.set_xticks(
            [np.log10(i + 1) for i in range(10)], labels=[str(i + 1) for i in range(10)]
        )
        ax.set_ylim([1, ym])
        ax.set_yscale("log")
        plt.legend()
        # plt.tight_layout()
        ax.set_box_aspect(1 / 1.618)
        plt.savefig(f"s_neval_{n_mol:03d}.eps")
        plt.savefig(f"s_neval_{n_mol:03d}.png")
        plt.clf()
        plt.cla()
        plt.close()


def plot_cc():
    cs = (
        colorsets.get_color_gradient(colorsets.uylw, colorsets.sorg, 5)[1:]
        + colorsets.get_color_gradient(colorsets.sorg, colorsets.uppl, 5)[1:]
    )
    for n_apm in [1, 2, 3, 10, 100, 1000]:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        mask = (df["n_apm"] == n_apm) & (df["n_sym"] == 1)
        alpha = np.unique(df["alpha"][mask])
        for n_mol, c, lw in zip(
            np.unique(df["n_mol"]),
            cs,
            [3.0, 3.0, 2.75, 2.5, 2.25, 2.0, 1.75, 1.5],
        ):
            print(df["nt_cc"][mask & (df["n_mol"] == n_mol)])
            ax.plot(
                alpha,
                df["nt_cc"][mask & (df["n_mol"] == n_mol)],
                label=f"n_mol {n_mol:3d}",
                c=c,
                lw=lw,
            )
        ax.set_title(f"t-neval cc n_apm {n_apm:03d}")
        ax.set_xlabel("alpha")
        ax.set_ylabel("Pearson cc")
        ax.set_ylim([-0.02, 1.02])
        ax.set_yticks(
            [0.2 * i for i in range(6)],
        )
        plt.legend()
        ax.set_box_aspect(1.618)
        plt.tight_layout()
        plt.savefig(f"cc_n{n_apm:04d}.eps")
        plt.savefig(f"cc_n{n_apm:04d}.png")
        plt.clf()
        plt.cla()
        plt.close()

    for alpha in [0.0, 0.6, 0.9]:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        mask = (df["n_apm"] == 10) & (df["alpha"] == alpha)
        n_sym = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for n_mol, c, lw in zip(
            [3, 4, 5, 6],
            cs[::2],
            [3.0, 2.5, 2.0, 1.5],
        ):
            print(df["nt_cc"][mask & (df["n_mol"] == n_mol)])
            ax.plot(
                n_sym,
                df["nt_cc"][mask & (df["n_mol"] == n_mol)],
                label=f"n_mol {n_mol:3d}",
                c=c,
                lw=lw,
            )
        ax.set_title(f"t-neval cc alpha {alpha:3.1f}")
        ax.set_xlabel("S")
        ax.set_ylabel("Pearson cc")
        ax.set_ylim([-0.02, 1.02])
        ax.set_yticks(
            [0.2 * i for i in range(6)],
        )
        plt.legend()
        ax.set_box_aspect(1.618)
        plt.tight_layout()
        plt.savefig(f"cc_a{alpha:3.1f}.eps")
        plt.savefig(f"cc_a{alpha:3.1f}.png")
        plt.clf()
        plt.cla()
        plt.close()


plot_curve()
plot_m_n()
plot_s_n()
plot_cc()
