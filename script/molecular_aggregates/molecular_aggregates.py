#!/bin/env python
import sys
import time
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
import mdtraj as md
import pandas as pd
import json
import mobbrmsd


def so3():
    v = np.random.randn(3)
    n = np.sum(np.power(v, 2))
    v /= np.sqrt(n)
    w1 = np.array([[0.0, v[2], -v[1]], [-v[2], 0.0, v[0]], [v[1], -v[0], 0.0]])
    w2 = -np.eye(3) + v.reshape([3, 1]) @ v.reshape([1, 3])
    return np.eye(3) + np.sin(n) * w1 + (1.0 - np.cos(n)) * w2


def multiple_mol(temp, n_mol):
    ret = temp
    for i in range(n_mol - 1):
        ret = Chem.rdmolops.CombineMols(ret, temp)
    return ret


def mol2xyz(temp, n_mol, nsample):
    conf = temp.GetConformer()
    xyz_ = []
    for i, a in enumerate(temp.GetAtoms()):
        p = conf.GetAtomPosition(i)
        xyz_ += [[p.x, p.y, p.z]]
    xyz = np.array(xyz_)
    xyz -= np.mean(xyz, 0)
    ret = np.empty([nsample, n_mol, len(xyz), 3])
    n13 = np.power(n_mol, 1 / 3)
    for i in range(nsample):
        for j in range(n_mol):
            ret[i, j, :] = xyz @ so3() + n13 * np.random.randn(3).reshape([1, 3])

    return ret.reshape([nsample, -1, 3])


def setconf(temp, xyz):
    conf = temp.GetConformer()
    for i, x in enumerate(xyz):
        conf.SetAtomPosition(i, Point3D(x[0], x[1], x[2]))


def get_permute_groups(js, n_mol):
    ret = []
    for pg in js["permute_groups"]:
        ph = []
        for j in range(n_mol):
            ph += [js["n_apm"] * j + i for i in pg]
        ret += [ph]
    return ret


def run(
    system,
    template,
    molsys,
    n_mol=1,
    nsample=10,
    permute_groups=None,
    skip_GetBestRMS=False,
):
    print(system, n_mol)

    temp = Chem.MolFromMolFile(template, removeHs=False)
    refs = mol2xyz(temp, n_mol, nsample)
    trgs = mol2xyz(temp, n_mol, nsample)
    mobb = mobbrmsd.mobbrmsd(molsys)

    top = md.Topology()
    dum = top.add_residue("DUM", top.add_chain())
    for i in range(n_mol * len(temp.GetAtoms())):
        top.add_atom("X" + f"{i:02d}", element=None, residue=dum)

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

    def run_GetBestRMS(ref, trg):
        reftemp = multiple_mol(temp, n_mol)
        trgtemp = multiple_mol(temp, n_mol)
        conf = reftemp.GetConformer()
        setconf(reftemp, ref)
        setconf(trgtemp, trg)
        start = time.perf_counter()
        r = AllChem.GetBestRMS(reftemp, trgtemp, maxMatches=100000000)
        t = time.perf_counter() - start
        return r, t

    def run_lp_with_ph(ref, trg):
        x = md.Trajectory(ref, topology=top)
        y = md.Trajectory(trg, topology=top)
        start = time.perf_counter()
        r = md.lprmsd(y, x, permute_groups=permute_groups)[0]
        t = time.perf_counter() - start
        return r, t

    def run_lp_none_ph(ref, trg):
        x = md.Trajectory(ref, topology=top)
        y = md.Trajectory(trg, topology=top)
        start = time.perf_counter()
        r = md.lprmsd(y, x)[0]
        t = time.perf_counter() - start
        return r, t

    def run_mobbrmsd(ref, trg):
        start = time.perf_counter()
        r = mobb.rmsd(ref, trg)
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
        stats += [[system, n_mol, tr, tp, tl, tm, rr, rp, rl, rm]]
    pd.DataFrame(
        stats,
        columns=[
            "system",
            "n_mol",
            "time_rdkit",
            "time_lprmsd_pg",
            "time_lprmsd",
            "time_mobbrmsd",
            "rdkit",
            "lprmsd_pg",
            "lprmsd",
            "mobbrmsd",
        ],
    ).to_pickle(Path(f"{system}_{n_mol:02d}.pkl.gz"), compression="gzip")


with open(sys.argv[1], "r") as f:
    js = json.load(f)

for n_mol, n_sample, skip_GetBestRMS in zip(
    js["n_mols"], js["n_samples"], js["skip_GetBestRMS"]
):
    run(
        js["name"],
        js["template"],
        molsys=mobbrmsd.dataclass.molecules(
            n_apm=js["n_apm"], n_mol=n_mol, sym=js["sym"]
        ),
        n_mol=n_mol,
        permute_groups=get_permute_groups(js, n_mol),
        skip_GetBestRMS=skip_GetBestRMS,
    )
