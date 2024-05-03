import numpy as np
import sys
from pathlib import Path
import matplotlib.pyplot as plt


def percentiles(x, q):
    c = 0.0
    for i, xi in enumerate(x):
        ci = c + xi
        if c > q:
            return i + (q - c) / xi
        c = ci


n = int(len(sys.argv[1:]) / 2)
for f1, f2 in zip(sys.argv[1:], sys.argv[n + 1 :]):
    a = np.loadtxt(f1)
    b = np.loadtxt(f2)
    x = np.sum(a, 0)[::-1] + np.sum(b, 0)[::-1]
    s = np.sum(x)
    p25 = percentiles(x, s / 4)
    p50 = percentiles(x, s / 2)
    p75 = percentiles(x, 3 * s / 4)
    mean = np.sum(x * np.linspace(1, x.size, num=x.size)) / s
    print(f1, f2, f"{p25:16.9f} {p50:16.9f} {p75:16.9f} {mean:16.9f}")
