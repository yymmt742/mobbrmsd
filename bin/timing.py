#! /usr/bin/env python3
import random
import numpy as np
import mobbrmsd as mo
import mdtraj as md
import time
import sys
import json


def pair_generator(length, number):
    used_pairs = set()
    for i in range(length):
        while True:
            pair = random.sample(range(number), 2)
            pair = tuple(sorted(pair))
            if pair not in used_pairs:
                used_pairs.add(pair)
                break
        yield pair


def run_with_time(xyz, nsample):
    start = time.perf_counter_ns()
    for i, j in pair_generator(nsample, xyz.shape[0]):
        ret = mrmsd.run(xyz[i], xyz[j])
    end = time.perf_counter_ns() - start
    return end * 1.0e-9


with open(sys.argv[1], "r") as f:
    j = json.load(f)
    dat = md.load_netcdf(j["trj"], top=j["top"])
    mrmsd = mo.mobbrmsd()
    mrmsd.add_molecule(j["n_apm"], dat.n_residues, swp=j["sym"])
    xyz = np.array(dat.xyz, dtype=np.float64)
    n_sample = j["n_sample"]
    time = run_with_time(xyz, n_sample)
    tps = 1.0e3 * time / n_sample
    print(f"{time:9.3f} sec with {n_sample:d} samples {tps:16.9f} msec/sample")
