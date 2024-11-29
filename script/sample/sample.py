from mobbrmsd import mobbrmsd, dataclass, coord_generator
import numpy as np
import time
import itertools
import random

cogen = coord_generator()


def run_test(n_apm, n_mol, sym, n_test, a, b):
    mrmsd = mobbrmsd(mols=dataclass.molecules(n_apm=n_apm, n_mol=n_mol, sym=sym))
    dat = np.empty([n_test, 3])
    ttime = 0.0
    for i in range(n_test):
        x = cogen.generate(n_apm, n_mol, a, b).reshape([-1, 3])
        y = cogen.generate(n_apm, n_mol, a, b).reshape([-1, 3])
        start = time.perf_counter()
        ret = mrmsd.run(x, y)
        end = time.perf_counter()
        dat[i, 0] = ret.rmsd()
        dat[i, 1] = ret.n_eval()
        dat[i, 2] = end - start
    return dat


n_apms = [1, 2, 3, 10, 100, 1000]
n_mols = [3, 4, 5, 6, 7, 8, 9, 10]
n_tests = [20000, 20000, 20000, 20000, 20000, 10000, 10000, 10000]
alphas = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
b = 0.0
n_sym = 0

for n_apm in n_apms:
    for n_mol, n_test in zip(n_mols, n_tests):
        for a in alphas:
            print(f"sample_{n_mol:02d}_{n_apm:03d}_{n_sym:03d}_{a:3.1f}")
            d = run_test(n_apm, n_mol, None, n_test, a, b)
            np.save(f"sample_{n_mol:02d}_{n_apm:03d}_{n_sym:03d}_{a:3.1f}", d)

n_apm = 10
n_mols = [3, 4, 5, 6]
S = 10
n_test = 20000
alphas = [0.0, 0.6, 0.9]
b = 0.0
prms = list(itertools.permutations(range(n_apm)))[1:]
for n_mol in n_mols:
    sym = []
    for i, i_sym in enumerate(random.sample(prms, S - 1)):
        sym += [list(i_sym)]
        if i < 4:
            continue
        n_sym = i + 1
        for a in alphas:
            print(f"sample_{n_mol:02d}_{n_apm:03d}_{n_sym:03d}_{a:3.1f}")
            dat = sample.run_test(n_apm, n_mol, sym, n_test, a, b)
            np.save(f"sample_{n_mol:02d}_{n_apm:03d}_{n_sym:03d}_{a:3.1f}", dat)
