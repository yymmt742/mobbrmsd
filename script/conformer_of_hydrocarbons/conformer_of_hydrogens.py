#!/bin/env python
import sys
import time
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
import mdtraj as md
import numpy as np
import mobbrmsd
import pandas as pd


def so3():
    v = np.random.randn(3)
    n = np.sum(np.power(v, 2))
    v /= np.sqrt(n)
    w1 = np.array([[0.0, v[2], -v[1]], [-v[2], 0.0, v[0]], [v[1], -v[0], 0.0]])
    w2 = -np.eye(3) + v.reshape([3, 1]) @ v.reshape([1, 3])
    return np.eye(3) + np.sin(n) * w1 + (1.0 - np.cos(n)) * w2


def setconf(temp, xyz):
    conf = temp.GetConformer()
    for i, x in enumerate(xyz):
        conf.SetAtomPosition(i, Point3D(x[0], x[1], x[2]))


def run(
    system,
    template,
    molsys,
    permutation=None,
    permute_groups=None,
    skip_GetBestRMS=False,
):
    print(system)
    dum = Chem.MolFromMolFile(template, removeHs=False)
    top = md.Topology()
    res = top.add_residue("DUM", top.add_chain())
    for i, a in enumerate(dum.GetAtoms()):
        top.add_atom("X" + f"{i:02d}", element=None, residue=res)

    refs = md.Trajectory(
        [xyz @ so3() * 10 for xyz in md.load_pdb(system + "_ref.pdb").xyz], topology=top
    )  # to angs.
    trgs = md.Trajectory(
        [xyz @ so3() * 10 for xyz in md.load_pdb(system + "_trg.pdb").xyz], topology=top
    )  # to angs.
    mobb = mobbrmsd.mobbrmsd(molsys)

    def runner(method, title, system):
        rmsds = []
        times = []
        i = 0
        print(f"# {title} ({system})")
        print("           wall-clock time (sec)    RMSD")
        for tref, ttrg in zip(refs, trgs):
            r, t = method(tref, ttrg)
            rmsds += [r]
            times += [t]
            i += 1
            print(f"{i:8d}{t:24.12f}{r:8.3f}")
        print()
        return rmsds, times

    def run_GetBestRMS(tref, ttrg):
        reftemp = Chem.MolFromMolFile(template, removeHs=False)
        trgtemp = Chem.MolFromMolFile(template, removeHs=False)
        conf = reftemp.GetConformer()
        setconf(reftemp, np.array(tref.xyz[0], dtype=np.float64))
        setconf(trgtemp, np.array(ttrg.xyz[0], dtype=np.float64))
        start = time.perf_counter()
        r = AllChem.GetBestRMS(reftemp, trgtemp, maxMatches=100000000)
        t = time.perf_counter() - start
        return r, t

    def run_lp_with_ph(tref, ttrg):
        start = time.perf_counter()
        r = md.lprmsd(ttrg, tref, permute_groups=permute_groups)[0]
        t = time.perf_counter() - start
        return r, t

    def run_lp_none_ph(tref, ttrg):
        start = time.perf_counter()
        r = md.lprmsd(ttrg, tref)[0]
        t = time.perf_counter() - start
        return r, t

    def run_mobbrmsd(tref, ttrg):
        x_ = (
            np.array(tref.xyz[0], dtype=np.float64)
            if permutation is None
            else np.array(tref.xyz[0], dtype=np.float64)[permutation]
        )
        y_ = (
            np.array(ttrg.xyz[0], dtype=np.float64)
            if permutation is None
            else np.array(ttrg.xyz[0], dtype=np.float64)[permutation]
        )
        start = time.perf_counter()
        r = mobb.rmsd(x_, y_)
        t = time.perf_counter() - start
        return r, t

    if skip_GetBestRMS:
        rd = [None] * len(refs)
        rd_time = [None] * len(refs)
    else:
        rd, rd_time = runner(run_GetBestRMS, "GetBestRMS", system)

    pp, pp_time = runner(run_lp_with_ph, "LP-RMSD (with permute_groups)", system)
    lp, lp_time = runner(run_lp_none_ph, "LP-RMSD (none permute_groups)", system)
    mb, mb_time = runner(run_mobbrmsd, "mobbRMSD", system)

    stats = []
    for tr, tp, tl, tm, rr, rp, rl, rm in zip(
        rd_time, pp_time, lp_time, mb_time, rd, pp, lp, mb
    ):
        stats += [[system, tr, tp, tl, tm, rr, rp, rl, rm]]
    pd.DataFrame(
        stats,
        columns=[
            "system",
            "time_rdkit",
            "time_lprmsd_pg",
            "time_lprmsd",
            "time_mobbrmsd",
            "rdkit",
            "lprmsd_pg",
            "lprmsd",
            "mobbrmsd",
        ],
    ).to_pickle(Path(f"conformer_of_hydrogens_{system}.pkl.gz"), compression="gzip")


with open(sys.argv[1], "r") as f:
    js = json.load(f)
skip_GetBestRMS = int(sys.argv[2]) == 0

run(
    js["name"],
    js["template"],
    molsys=mobbrmsd.dataclass.molecular_system(
        [
            mobbrmsd.dataclass.molecules(n_apm=n_apm, n_mol=n_mol, sym=sym)
            for n_apm, n_mol, sym in zip(js["n_apm"], js["n_mol"], js["sym"])
        ],
    ),
    permutation=js["permutation"],
    permute_groups=js["permute_groups"],
    skip_GetBestRMS=skip_GetBestRMS,
)
