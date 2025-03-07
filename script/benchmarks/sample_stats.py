import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import colorsets


def linfit(x, y):
    sxx = np.sum(np.power(x, 2))
    sxy = np.sum(x * y)
    sx = np.sum(x)
    mx = np.mean(x)
    sy = np.sum(y)
    my = np.mean(y)
    a = (sxy - sx * my) / (sxx - sx * mx)
    return a, my - a * mx


stats = []

for a in sys.argv[1:]:
    x = np.load(a, allow_pickle=True)
    m = np.mean(x, 0)
    s = np.std(x, 0)
    q = [
        np.percentile(x[:, 0], [25, 50, 75]),
        np.percentile(x[:, 1], [25, 50, 75]),
        np.percentile(x[:, 2], [25, 50, 75]),
    ]
    t = [[mi] + [si] + qi.tolist() for mi, si, qi in zip(m, s, q)]
    sn = (x[:, 1] - m[1]) / s[1]
    st = (x[:, 2] - m[2]) / s[2]
    c = np.mean(sn * st)
    ca, cb = linfit(x[:, 1], x[:, 2])
    stem = Path(a).stem.split("_")
    stats += [
        [int(stem[1]), int(stem[2]), int(stem[3]) + 1, float(stem[4])]
        + t[0]
        + t[1]
        + t[2]
        + [c, ca, cb]
    ]

pd.DataFrame(
    stats,
    columns=[
        "n_mol",
        "n_apm",
        "n_sym",
        "alpha",
        "r_mean",
        "r_var",
        "r_q25",
        "r_q50",
        "r_q75",
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
        "nt_cc",
        "nt_ca",
        "nt_cb",
    ],
).to_pickle(Path(sys.argv[1]).stem[:-11] + "stats.pkl.gz", compression="gzip")
