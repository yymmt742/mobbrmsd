#!/bin/env python
import sys
import mdtraj as md
import numpy as np
import mobbrmsd
import matplotlib.pyplot as plt
import matplotlib.colors as color
import itertools
import random
import time
import pandas as pd


def brute(xyz):
    ret = np.empty([xyz.shape[0], xyz.shape[0]])
    trj1 = md.Trajectory(xyz, topology=top)
    ret[:, :] = np.inf
    for i in range(xyz.shape[0]):
        for perm in itertools.permutations([i for i in range(xyz.shape[1])]):
            trj0 = md.Trajectory(xyz[i, perm, :], topology=top)
            rmsd = md.rmsd(trj1, trj0)
            ret[i, :] = np.where(rmsd < ret[i], rmsd, ret[i])
    return ret


rng = np.random.default_rng()

n_sample = 50000
stats = []

for n_atom in [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]:
    top = md.Topology()
    res = top.add_residue("DUM", top.add_chain())
    for i in range(n_atom):
        top.add_atom("C" + f"{i:02d}", element=None, residue=res)

    for i in range(n_sample):
        x = rng.standard_normal([n_atom, 3])
        x -= np.mean(x, 0)
        y = rng.standard_normal([n_atom, 3])
        y -= np.mean(y, 0)
        start = time.perf_counter()
        rmsd = md.lprmsd(
            md.Trajectory(x, topology=top), md.Trajectory(y, topology=top)
        )[0]
        t = time.perf_counter() - start
        stats += [[n_atom, t]]

    print(f"n_mol {n_atom:2d}")

pd.DataFrame(
    stats,
    columns=[
        "n_mol",
        "time_lprmsd",
    ],
).to_pickle("vs_lprmsd_time.pkl.gz", compression="gzip")

n_sample = 10000
stats = []

for n_atom in [3, 4, 5, 6, 7, 8, 9, 10]:
    top = md.Topology()
    res = top.add_residue("DUM", top.add_chain())
    for i in range(n_atom):
        top.add_atom("C" + f"{i:02d}", element=None, residue=res)

    mobb = mobbrmsd.mobbrmsd(
        mols=mobbrmsd.dataclass.molecules(
            n_apm=1,
            n_mol=n_atom,
        ),
    )

    """
    brut = brute(xyz)
    """
    mdlp = np.empty([n_sample])
    mobr = np.empty([n_sample])

    for i in range(n_sample):
        x = rng.standard_normal([n_atom, 3])
        x -= np.mean(x, 0)
        y = rng.standard_normal([n_atom, 3])
        y -= np.mean(y, 0)
        mdlp[i] = md.lprmsd(
            md.Trajectory(x, topology=top), md.Trajectory(y, topology=top)
        )[0]
        mobr[i] = mobb.rmsd(x, y)
        stats += [[n_atom, mdlp[i], mobr[i]]]

    diff = mdlp - mobr
    mean = np.mean(diff)
    std = np.sqrt(np.mean(np.power(diff - mean, 2)))
    cong = np.count_nonzero(diff < 1.0e-4)
    print(
        f"n_mol {n_atom:2d} mean error {mean:16.9f} std dev {std:16.9f} count(diff<1.e-4) {cong} / {n_sample}",
        flush=True,
    )

    ix = random.sample(range(n_sample), 300)
    plt.plot(np.linspace(0, 2, 100), np.linspace(0, 2, 100), c="gray")
    plt.scatter(
        mdlp[ix],
        mobr[ix],
        s=8.0,
        facecolors="none",
        edgecolors="r",
        lw=0.1,
    )
    plt.xlim([0, 2.0])
    plt.ylim([0, 2.0])
    plt.gca().set_aspect(1.0)
    plt.savefig(f"vs_lprmsd_{n_atom:02d}.png")
    plt.savefig(f"vs_lprmsd_{n_atom:02d}.eps")
    plt.clf()
    plt.close()

pd.DataFrame(
    stats,
    columns=[
        "n_mol",
        "lprmsd",
        "mobbrmsd",
    ],
).to_pickle("vs_lprmsd.pkl.gz", compression="gzip")
