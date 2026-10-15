"""Accuracy and timing comparison: mobbrmsd vs RDKit GetBestRMS.

Run with::

    pytest tests/test_perf_rdkit.py --run-perf -v -s

Requires: ``pip install rdkit``

Test strategy
-------------
For ``n_apm=1, n_mol=k`` (k identical disconnected dummy atoms), both
``mobbrmsd`` and ``AllChem.GetBestRMS`` explore the full k! permutation space
and must agree on the optimal (rotation + permutation) RMSD.  The dummy-atom
molecule is built programmatically—no .mol file required.
"""
import time

import numpy
import pytest

pytestmark = pytest.mark.perf

_N_SAMPLES = 50  # enough for timing estimates, fast for accuracy checks


# ---------------------------------------------------------------------------
# RDKit molecule helpers
# ---------------------------------------------------------------------------

def _make_dummy_mol(n: int):
    """n disconnected dummy atoms — all equivalent, GetBestRMS tries all n! perms."""
    from rdkit import Chem

    mol = Chem.RWMol()
    for _ in range(n):
        at = Chem.Atom(0)       # atomic number 0 = wildcard/dummy
        at.SetNoImplicit(True)
        mol.AddAtom(at)
    conf = Chem.Conformer(n)
    for i in range(n):
        conf.SetAtomPosition(i, (0.0, 0.0, 0.0))
    mol.AddConformer(conf, assignId=True)
    return mol.GetMol()


def _set_coords(mol, xyz) -> None:
    """Update conformer coordinates in-place (xyz shape: (n, 3))."""
    from rdkit.Geometry import Point3D

    conf = mol.GetConformer()
    for i, (x, y, z) in enumerate(xyz):
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))


def _rdkit_best_rmsd(mol_ref, mol_trg, xyz_ref, xyz_trg) -> float:
    from rdkit.Chem import AllChem

    _set_coords(mol_ref, xyz_ref)
    _set_coords(mol_trg, xyz_trg)
    return AllChem.GetBestRMS(mol_ref, mol_trg, maxMatches=100_000_000)


# ---------------------------------------------------------------------------
# Accuracy tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_mol", [2, 3, 4, 5, 6])
def test_accuracy_vs_rdkit(n_mol):
    """mobbrmsd and RDKit GetBestRMS agree within 1e-4 for n_apm=1, varying n_mol."""
    pytest.importorskip("rdkit")
    from mobbrmsd.dataclass import molecules
    from mobbrmsd.mobbrmsd import mobbrmsd

    rng = numpy.random.default_rng(seed=n_mol * 1000 + 7)
    mrmsd = mobbrmsd(mols=molecules(n_apm=1, n_mol=n_mol))
    mol_ref = _make_dummy_mol(n_mol)
    mol_trg = _make_dummy_mol(n_mol)

    max_diff = 0.0
    for _ in range(_N_SAMPLES):
        x = rng.standard_normal((n_mol, 3))
        y = rng.standard_normal((n_mol, 3))
        r_mobb = mrmsd.rmsd(x, y)
        r_rd   = _rdkit_best_rmsd(mol_ref, mol_trg, x, y)
        max_diff = max(max_diff, abs(r_mobb - r_rd))

    assert max_diff < 1e-4, (
        f"n_mol={n_mol}: max |mobbrmsd − rdkit| = {max_diff:.2e}  (threshold 1e-4)"
    )


# ---------------------------------------------------------------------------
# Timing comparison (informational — no RMSD assertion)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_mol", [2, 3, 4, 5, 6, 7, 8])
def test_timing_vs_rdkit(n_mol, capsys):
    """Print wall-time per pair for mobbrmsd vs RDKit (no pass/fail threshold)."""
    pytest.importorskip("rdkit")
    from mobbrmsd.dataclass import molecules
    from mobbrmsd.mobbrmsd import mobbrmsd

    rng = numpy.random.default_rng(seed=42)
    mrmsd = mobbrmsd(mols=molecules(n_apm=1, n_mol=n_mol))
    mol_ref = _make_dummy_mol(n_mol)
    mol_trg = _make_dummy_mol(n_mol)

    xs = [rng.standard_normal((n_mol, 3)) for _ in range(_N_SAMPLES)]
    ys = [rng.standard_normal((n_mol, 3)) for _ in range(_N_SAMPLES)]

    # warm-up
    _rdkit_best_rmsd(mol_ref, mol_trg, xs[0], ys[0])
    mrmsd.rmsd(xs[0], ys[0])

    t0 = time.perf_counter()
    for x, y in zip(xs, ys):
        _rdkit_best_rmsd(mol_ref, mol_trg, x, y)
    t_rd = (time.perf_counter() - t0) / _N_SAMPLES * 1e6   # µs/pair

    t0 = time.perf_counter()
    for x, y in zip(xs, ys):
        mrmsd.rmsd(x, y)
    t_mb = (time.perf_counter() - t0) / _N_SAMPLES * 1e6

    with capsys.disabled():
        print(
            f"\n  n_mol={n_mol:3d}:  mobbrmsd={t_mb:8.2f} µs/pair"
            f"  rdkit={t_rd:8.2f} µs/pair"
            f"  ratio rdkit/mobb={t_rd / max(t_mb, 1e-9):.1f}x"
        )
