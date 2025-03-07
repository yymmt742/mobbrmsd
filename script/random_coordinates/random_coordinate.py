import os
import time
import mobbrmsd
import subprocess
import mdtraj as md
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
import pandas as pd
import numpy as np
from pathlib import Path

# Must be execute in X-GRMSD/build for test IsometryOpt and MatchFastOpt.


def multiple_mol(temp, nmol):
    ret = temp
    for i in range(nmol - 1):
        ret = Chem.rdmolops.CombineMols(ret, temp)
    return ret


def setconf(temp, xyz):
    conf = temp.GetConformer()
    for i, x in enumerate(xyz):
        conf.SetAtomPosition(i, Point3D(x[0], x[1], x[2]))


rng = np.random.default_rng()

grmsd_is_exist = os.path.isfile("./GRMSD")
rdkit_upperbound = 8


def run(n_mol):
    x = rng.standard_normal((n_mol, 3))
    y = rng.standard_normal((n_mol, 3))
    mrmsd = mobbrmsd.mobbrmsd(mols=mobbrmsd.dataclass.molecules(n_mol=n_mol, n_apm=1))

    top = md.Topology()
    dum = top.add_residue("DUM", top.add_chain())
    for i in range(n_mol):
        top.add_atom("X" + f"{i:02d}", element=None, residue=dum)
    x_ = md.Trajectory(x, topology=top)
    y_ = md.Trajectory(y, topology=top)

    def run_grmsd(method):
        start = time.perf_counter()
        subprocess.call(f"./GRMSD {method} model.csv data.csv 0 out 0", shell=True)
        t = time.perf_counter() - start
        return (
            float(
                (
                    (
                        subprocess.Popen(
                            "head -1 out", stdout=subprocess.PIPE, shell=True
                        ).communicate()[0]
                    )
                    .decode("utf-8")
                    .split(" ")[-1]
                )
            ),
            t,
        )

    def gen_csv(path, x):
        with open(path, "w") as f:
            f.write(f"{n_mol},{n_mol},{n_mol},{n_mol}")
            for xi in x:
                f.write(f"\n{xi[0]},{xi[1]},{xi[2]},0")

    if grmsd_is_exist:
        gen_csv("model.csv", x)
        gen_csv("data.csv", y)
        grmsd, time_grmsd = run_grmsd("IsometryOpt")
        frmsd, time_frmsd = run_grmsd("MatchFastOpt")
    else:
        grmsd, time_grmsd, frmsd, time_frmsd = (None, None, None, None)

    start = time.perf_counter()
    lprmsd = md.lprmsd(y_, x_)[0]
    time_lprmsd = time.perf_counter() - start
    start = time.perf_counter()
    mobb = mrmsd.rmsd(x, y)
    time_mobbrmsd = time.perf_counter() - start

    if n_mol <= rdkit_upperbound:
        reftemp = multiple_mol(Chem.MolFromMolFile("dumm.mol", removeHs=False), n_mol)
        trgtemp = multiple_mol(Chem.MolFromMolFile("dumm.mol", removeHs=False), n_mol)
        setconf(reftemp, x)
        setconf(trgtemp, y)
        start = time.perf_counter()
        rdkit = AllChem.GetBestRMS(reftemp, trgtemp, maxMatches=100000000)
        time_rdkit = time.perf_counter() - start
    else:
        rdkit, time_rdkit = (None, None)

    return [
        "random",
        n_mol,
        time_grmsd,
        time_frmsd,
        time_rdkit,
        time_lprmsd,
        time_mobbrmsd,
        grmsd,
        frmsd,
        rdkit,
        lprmsd,
        mobb,
    ]


n_samples = [
    10,  # 2
    10,  # 3
    10,  # 4
    10,  # 5
    10,  # 6
    10,  # 7
    10,  # 8
    10,  # 9
    10,  # 10
]

print("     n_mol        n_sample")
for i, n_sample in enumerate(n_samples):
    stats = []
    n_mol = i + 2
    print(f"{i:2d}{n_mol:8d}{n_sample:16d}")
    for j in range(n_sample):
        stats += [run(n_mol)]
    pd.DataFrame(
        stats,
        columns=[
            "system",
            "n_mol",
            "time_grmsd",
            "time_frmsd",
            "time_rdkit",
            "time_lprmsd",
            "time_mobbrmsd",
            "grmsd",
            "frmsd",
            "rdkit",
            "lprmsd",
            "mobbrmsd",
        ],
    ).to_pickle(Path(f"random_coordinates_{n_mol:02d}.pkl.gz"), compression="gzip")
