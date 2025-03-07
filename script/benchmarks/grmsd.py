import mobbrmsd as mobb
import time
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path

# Must be execute in X-GRMSD/build.

rng = np.random.default_rng()


def run(n_mol):
    x = rng.standard_normal((n_mol, 3))
    y = rng.standard_normal((n_mol, 3))
    mrmsd = mobb.mobbrmsd(mols=mobb.dataclass.molecules(n_mol=n_mol, n_apm=1))

    def gen_csv(path, x):
        with open(path, "w") as f:
            f.write(f"{n_mol},{n_mol},{n_mol},{n_mol}")
            for xi in x:
                f.write(f"\n{xi[0]},{xi[1]},{xi[2]},0")

    gen_csv("model.csv", x)
    gen_csv("data.csv", y)
    start = time.perf_counter()
    subprocess.call("./GRMSD IsometryOpt model.csv data.csv 0 out 0", shell=True)
    time_grmsd = time.perf_counter() - start
    grmsd = float(
        (
            (
                subprocess.Popen(
                    "head -1 out", stdout=subprocess.PIPE, shell=True
                ).communicate()[0]
            )
            .decode("utf-8")
            .split(" ")[-1]
        )
    )
    start = time.perf_counter()
    ret = mrmsd.run(x, y)
    time_mobbrmsd = time.perf_counter() - start
    return [n_mol, time_grmsd, time_mobbrmsd, grmsd, ret.rmsd()]


n_samples = [
    10000,  # 2
    10000,  # 3
    10000,  # 4
    10000,  # 5
    10000,  # 6
    10000,  # 7
    1000,  # 8
    1000,  # 9
    1000,  # 10
    100,  # 10
    100,  # 11
    100,  # 12
    100,  # 13
    100,  # 14
    20,  # 15
    20,  # 16
    20,  # 17
    20,  # 18
    20,  # 19
    10,  # 20
]

for i, n_sample in enumerate(n_samples):
    stats = []
    n_mol = i + 2
    for j in range(n_sample):
        stats += [run(n_mol)]
    print(f"save : vs_grmsd_{n_mol:02d}.pkl.gz")
    pd.DataFrame(
        stats, columns=["n_mol", "time_grmsd", "time_mobbrmsd", "grmsd", "mobbrmsd"]
    ).to_pickle(Path(f"vs_grmsd_{n_mol:02d}.pkl.gz"), compression="gzip")
