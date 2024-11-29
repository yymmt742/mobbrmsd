#! /usr/bin/env python3
import random
import numpy as np
import pandas as pd
import mobbrmsd as mo
import time
import sys
import json
import matplotlib.pyplot as plt
import seaborn as sns


loge = np.log(np.e)


def run_with_time_random(mrmsd, cogen, n_mol, n_apm, a, b, n_sample):
    tims = []
    rmsd = []
    n_eval = []
    ratio = []
    rng = np.random.default_rng()
    for i in range(n_sample):
        x = cogen.generate(n_apm, n_mol, a, b).reshape([-1, 3])
        y = cogen.generate(n_apm, n_mol, a, b).reshape([-1, 3])
        start = time.perf_counter_ns()
        ret = mrmsd.run(x, y)
        tims += [(time.perf_counter_ns() - start) * 1.0e-6]
        rmsd += [ret.rmsd()]
        n_eval += [ret.n_eval()]
        ratio += [loge * ret.log_eval_ratio()]
    return tims, rmsd, n_eval, ratio


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


cogen = mo.coord_generator()
for arg in sys.argv[1:]:
    with open(arg, "r") as f:
        js = json.load(f)
        n_sample = js["n_sample"] // 10
        mrmsd = mo.mobbrmsd(
            molecules=mo.DataclassMolecule(
                n_apm=js["n_apm"],
                n_mol=js["n_mol"],
                sym=js["sym"],
            )
        )
        tims, rmsd, n_eval, ratio = run_with_time_random(
            mrmsd, cogen, js["n_mol"], js["n_apm"], 0.0, 0.0, n_sample
        )
        tim = np.sum(tims) * 1.0e-3
        tps = np.mean(tims)
        print(
            js["sys"],
            "-",
            js["n_mol"],
            f"{tim:9.3f} sec / {tps:9.6f} msec/sample with {n_sample:d} samples",
            flush=True,
        )
        data["rmsd"] += rmsd
        data["ulim"] += rmsd
        data["llim"] += rmsd
        data["times"] += tims
        data["n_eval"] += n_eval
        data["ratio"] += ratio
        data["sys"] += ["dummy" for i in range(n_sample)]
        data["n_mol"] += [js["n_mol"] for i in range(n_sample)]
        data["coff"] += [1000.0 for i in range(n_sample)]
        data["uoff"] += [1000.0 for i in range(n_sample)]
    del mrmsd

pd.DataFrame(data).to_pickle("timing_dummy_data.pkl.gz", compression="gzip")
