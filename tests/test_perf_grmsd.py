"""Accuracy and timing comparison: mobbrmsd vs GRMSD binary.

Run with::

    GRMSD_BIN=/path/to/GRMSD pytest tests/test_perf_grmsd.py --run-perf -v -s

Requires: GRMSD binary compiled from https://github.com/striwata/GRMSD
          Set the ``GRMSD_BIN`` environment variable to the binary path,
          or place the binary on PATH as ``GRMSD`` (``GRMSD.exe`` on Windows).
          Tests are skipped automatically when the binary is not found.

CSV format (from grmsd.py in mobbrmsd_testcodes)
-------------------------------------------------
First line: ``n,n,n,n``  (n = number of atoms)
Following lines: ``x,y,z,0``  per atom (atom type 0 = unlabelled)
Output: first line of the ``out`` file; last whitespace-delimited token = RMSD.
"""
import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy
import pytest

pytestmark = pytest.mark.perf

_N_SAMPLES = 20   # each sample spawns a subprocess; keep low to avoid long runtime


# ---------------------------------------------------------------------------
# GRMSD binary helpers
# ---------------------------------------------------------------------------

def _find_grmsd_bin() -> str:
    """Return path to GRMSD binary, or call pytest.skip if not found."""
    candidate = os.environ.get("GRMSD_BIN", "GRMSD")
    if shutil.which(candidate) or (os.path.isfile(candidate) and os.access(candidate, os.X_OK)):
        return candidate
    # Windows: also try with .exe suffix
    if shutil.which(candidate + ".exe"):
        return candidate + ".exe"
    pytest.skip(
        f"GRMSD binary not found (checked PATH and GRMSD_BIN={candidate!r}). "
        "Build from https://github.com/striwata/GRMSD and set GRMSD_BIN."
    )


def _write_grmsd_csv(path: Path, xyz) -> None:
    """Write coordinate array to GRMSD CSV format."""
    n = len(xyz)
    with open(path, "w") as f:
        f.write(f"{n},{n},{n},{n}")
        for x, y, z in xyz:
            f.write(f"\n{x:.10f},{y:.10f},{z:.10f},0")


def _read_grmsd_result(out_path: Path) -> float:
    """Read RMSD value from GRMSD output file (first line, last token)."""
    with open(out_path, "r") as f:
        line = f.readline()
    return float(line.split()[-1])


def _run_grmsd(xyz_ref, xyz_trg, workdir: Path, grmsd_bin: str) -> float:
    """Write CSVs, run GRMSD binary, return RMSD."""
    model = workdir / "model.csv"
    data  = workdir / "data.csv"
    out   = workdir / "out"
    _write_grmsd_csv(model, xyz_ref)
    _write_grmsd_csv(data,  xyz_trg)
    cmd = [grmsd_bin, "IsometryOpt", "model.csv", "data.csv", "0", "out", "0"]
    subprocess.run(cmd, cwd=workdir, capture_output=True, check=False)
    return _read_grmsd_result(out)


# ---------------------------------------------------------------------------
# Accuracy tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_mol", [2, 3, 4, 5, 6, 7, 8])
def test_accuracy_vs_grmsd(n_mol, tmp_path):
    """mobbrmsd and GRMSD agree within 1e-4 for n_apm=1, varying n_mol."""
    from mobbrmsd.dataclass import molecules
    from mobbrmsd.mobbrmsd import mobbrmsd

    grmsd_bin = _find_grmsd_bin()
    rng = numpy.random.default_rng(seed=n_mol * 1000 + 13)
    mrmsd = mobbrmsd(mols=molecules(n_apm=1, n_mol=n_mol))

    max_diff = 0.0
    for _ in range(_N_SAMPLES):
        x = rng.standard_normal((n_mol, 3))
        y = rng.standard_normal((n_mol, 3))
        r_mobb  = mrmsd.rmsd(x, y)
        r_grmsd = _run_grmsd(x, y, tmp_path, grmsd_bin)
        max_diff = max(max_diff, abs(r_mobb - r_grmsd))

    assert max_diff < 1e-4, (
        f"n_mol={n_mol}: max |mobbrmsd − grmsd| = {max_diff:.2e}  (threshold 1e-4)"
    )


# ---------------------------------------------------------------------------
# Timing comparison (informational — no pass/fail threshold)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_mol", [2, 4, 6, 8, 10])
def test_timing_vs_grmsd(n_mol, tmp_path, capsys):
    """Print wall-time per pair for mobbrmsd vs GRMSD (subprocess overhead included)."""
    from mobbrmsd.dataclass import molecules
    from mobbrmsd.mobbrmsd import mobbrmsd

    grmsd_bin = _find_grmsd_bin()
    rng = numpy.random.default_rng(seed=42)
    mrmsd = mobbrmsd(mols=molecules(n_apm=1, n_mol=n_mol))

    xs = [rng.standard_normal((n_mol, 3)) for _ in range(_N_SAMPLES)]
    ys = [rng.standard_normal((n_mol, 3)) for _ in range(_N_SAMPLES)]

    t0 = time.perf_counter()
    for x, y in zip(xs, ys):
        _run_grmsd(x, y, tmp_path, grmsd_bin)
    t_gr = (time.perf_counter() - t0) / _N_SAMPLES * 1e3   # ms/pair

    t0 = time.perf_counter()
    for x, y in zip(xs, ys):
        mrmsd.rmsd(x, y)
    t_mb = (time.perf_counter() - t0) / _N_SAMPLES * 1e6   # µs/pair

    with capsys.disabled():
        print(
            f"\n  n_mol={n_mol:3d}:  mobbrmsd={t_mb:8.2f} µs/pair"
            f"  grmsd={t_gr:8.2f} ms/pair (incl. subprocess)"
            f"  ratio grmsd/mobb={t_gr * 1e3 / max(t_mb, 1e-9):.0f}x"
        )
