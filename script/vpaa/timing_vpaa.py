#! /usr/bin/env python3
import random
import numpy as np
import pandas as pd
import mobbrmsd as mo
import mdtraj as md
import time
import sys
import json
import matplotlib.pyplot as plt
import seaborn as sns


c = np.eye(3)
c[2, 2] = -1
loge = np.log(np.e)


def pair_generator(length, number):
    for i in range(length):
        pair = tuple(random.sample(range(number), 2))
        yield pair


def min_rmsd(mrmsd, x, y):
    ret1 = mrmsd.run(x, y)
    ret2 = mrmsd.run(x, y @ c)
    if ret1.rmsd < ret2.rmsd:
        return ret1
    else:
        return ret2


def run_with_time(js, mrmsd, xyz, n_sample, cutoff, ub_cutoff):
    tims = []
    rmsd = []
    ulim = []
    llim = []
    coff = []
    uoff = []
    n_eval = []
    ratio = []
    n_traj = xyz.shape[0]
    tt0 = time.perf_counter_ns()
    for k, ij in enumerate(pair_generator(n_sample, xyz.shape[0])):
        for cut in cutoff:
            start = time.perf_counter_ns()
            ret = mrmsd.run(xyz[ij[0]], xyz[ij[1]], cutoff=cut)
            tims += [(time.perf_counter_ns() - start) * 1.0e-6]
            rmsd += [ret.rmsd]
            ulim += [np.sqrt(max(0.0, ret.autocorr + 2 * ret.bounds[0]) / mrmsd.natom)]
            llim += [np.sqrt(max(0.0, ret.autocorr + 2 * ret.bounds[1]) / mrmsd.natom)]
            coff += [cut]
            uoff += [1000.0]
            n_eval += [ret.n_eval]
            ratio += [loge * ret.log_eval_ratio]
        for cut in ub_cutoff:
            start = time.perf_counter_ns()
            ret = mrmsd.run(xyz[ij[0]], xyz[ij[1]], ub_cutoff=cut)
            tims += [(time.perf_counter_ns() - start) * 1.0e-6]
            rmsd += [ret.rmsd]
            ulim += [np.sqrt(max(0.0, ret.autocorr + 2 * ret.bounds[0]) / mrmsd.natom)]
            llim += [np.sqrt(max(0.0, ret.autocorr + 2 * ret.bounds[1]) / mrmsd.natom)]
            coff += [1000.0]
            uoff += [cut]
            n_eval += [ret.n_eval]
            ratio += [loge * ret.log_eval_ratio]

        if k % 20000 == 0:
            tps = np.mean(tims)
            tim = (time.perf_counter_ns() - tt0) * 1.0e-9
            print(
                js["sys"],
                "-",
                js["n_mol"],
                f"{k:8d} /{n_sample:d}",
                f"{tim:9.3f} sec /{tps:9.6f} msec/run",
                flush=True,
            )
    tps = np.mean(tims)
    tim = (time.perf_counter_ns() - tt0) * 1.0e-9
    print(
        js["sys"],
        "-",
        js["n_mol"],
        f"{tim:9.3f} sec / {tps:9.6f} msec/run : sample with {n_sample:d} samples from {n_traj:d} frames",
        flush=True,
    )
    return tims, rmsd, n_eval, ratio, ulim, llim, coff, uoff


def run_json(arg, n_sample):
    data = {
        "rmsd": [],
        "ulim": [],
        "llim": [],
        "times": [],
        "n_eval": [],
        "ratio": [],
        "sys": [],
        "n_mol": [],
        "coff": [],
        "uoff": [],
    }
    with open(arg, "r") as f:
        js = json.load(f)
        prm = md.load_prmtop(js["top"])
        sub = prm.subset(prm.select(js["mask"]))
        dat = md.load_netcdf(js["trj"], top=sub)
        vsys = js["sys"]
        n_mol = js["n_mol"]
        xyz = np.array(dat.xyz, dtype=np.float64)
        mrmsd = mo.mobbrmsd(
            molecules=mo.DataclassMolecule(
                n_apm=js["n_apm"],
                n_mol=js["n_mol"],
                sym=js["sym"],
            )
        )
        tims, rmsd, n_eval, ratio, ulim, llim, coff, uoff = run_with_time(
            js,
            mrmsd,
            xyz,
            n_sample,
            cutoff=[1000.0, 0.3, 0.2, 0.15, 0.1, 0.0],
            ub_cutoff=[0.8, 0.6, 0.5, 0.4],
        )
        data["rmsd"] += rmsd
        data["ulim"] += ulim
        data["llim"] += llim
        data["times"] += tims
        data["n_eval"] += n_eval
        data["ratio"] += ratio
        data["sys"] += [vsys for i in range(len(rmsd))]
        data["n_mol"] += [n_mol for i in range(len(rmsd))]
        data["coff"] += coff
        data["uoff"] += uoff
        del mrmsd
        pd.DataFrame(data).to_pickle(
            f"timing_data_{vsys}_{n_mol:02d}.pkl.gz", compression="gzip"
        )


run_json(sys.argv[1], sys.argv[2])
