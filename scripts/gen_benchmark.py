#!/usr/bin/env python
"""Generate benchmark plots comparing mobbrmsd with RDKit and GRMSD.

Outputs SVG figures to ``docs/benchmark/img/`` and a JSON summary to
``docs/benchmark/results.json``.  Run this script before ``mkdocs build``
to keep the benchmark page up to date.

Usage::

    # RDKit only (requires: pip install rdkit)
    python scripts/gen_benchmark.py

    # Include GRMSD binary comparison
    GRMSD_BIN=/path/to/GRMSD python scripts/gen_benchmark.py

    # More samples for smoother statistics
    python scripts/gen_benchmark.py --n-samples 500
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── paths ──────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parents[1]
IMG_DIR = REPO_ROOT / "docs" / "benchmark" / "img"
OUT_JSON = REPO_ROOT / "docs" / "benchmark" / "results.json"

# ── colour palette (matching testcodes/sample/colorsets.py) ────────────────
C_MOBB  = "#FF9900"   # orange
C_RD    = "#9A0079"   # purple
C_GR    = "#35A16B"   # green
C_GRAY  = "#7F878F"

plt.rcParams.update({
    "font.size": 9,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
})


# ── RDKit helpers ───────────────────────────────────────────────────────────

def _make_dummy_mol(n: int):
    from rdkit import Chem
    mol = Chem.RWMol()
    for _ in range(n):
        at = Chem.Atom(0)
        at.SetNoImplicit(True)
        mol.AddAtom(at)
    conf = Chem.Conformer(n)
    for i in range(n):
        conf.SetAtomPosition(i, (0.0, 0.0, 0.0))
    mol.AddConformer(conf, assignId=True)
    return mol.GetMol()


def _set_rdkit_coords(mol, xyz) -> None:
    from rdkit.Geometry import Point3D
    conf = mol.GetConformer()
    for i, (x, y, z) in enumerate(xyz):
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))


def rdkit_best_rmsd(mol_ref, mol_trg, xyz_ref, xyz_trg) -> float:
    from rdkit.Chem import AllChem
    _set_rdkit_coords(mol_ref, xyz_ref)
    _set_rdkit_coords(mol_trg, xyz_trg)
    return AllChem.GetBestRMS(mol_ref, mol_trg, maxMatches=100_000_000)


# ── GRMSD helpers ───────────────────────────────────────────────────────────

def find_grmsd_bin() -> str | None:
    candidate = os.environ.get("GRMSD_BIN", "GRMSD")
    if shutil.which(candidate) or (
        os.path.isfile(candidate) and os.access(candidate, os.X_OK)
    ):
        return candidate
    if shutil.which(candidate + ".exe"):
        return candidate + ".exe"
    return None


def _find_cmake() -> str | None:
    """Return cmake executable path, including the venv's Scripts directory."""
    cmake = shutil.which("cmake")
    if cmake:
        return cmake
    # cmake Python package installs into the venv's Scripts dir, which may not
    # be on PATH when the script is invoked as `python scripts/gen_benchmark.py`
    venv_scripts = Path(sys.executable).parent
    for name in ("cmake.exe", "cmake"):
        candidate = venv_scripts / name
        if candidate.exists():
            return str(candidate)
    return None


def build_grmsd(build_dir: Path) -> str | None:
    """Download and build X-GRMSD via cmake; return path to built binary."""
    cmake = _find_cmake()
    if not cmake:
        print("cmake not found — cannot build GRMSD", file=sys.stderr)
        return None

    src_dir = REPO_ROOT / "third_party" / "grmsd"
    if not (src_dir / "CMakeLists.txt").exists():
        print(f"third_party/grmsd/CMakeLists.txt not found", file=sys.stderr)
        return None

    build_dir.mkdir(parents=True, exist_ok=True)
    print(f"Configuring X-GRMSD (build dir: {build_dir}) …")

    # On Windows without MSVC, cmake needs an explicit generator and make tool.
    cfg_args = [cmake, "-S", str(src_dir), "-B", str(build_dir),
                "-DCMAKE_BUILD_TYPE=Release"]
    cxx = shutil.which("g++") or shutil.which("c++")
    make = shutil.which("mingw32-make") or shutil.which("make") or shutil.which("gmake")
    if sys.platform == "win32" and cxx:
        cfg_args += ["-G", "MinGW Makefiles",
                     f"-DCMAKE_CXX_COMPILER={cxx}"]
        if make:
            cfg_args += [f"-DCMAKE_MAKE_PROGRAM={make}"]

    r = subprocess.run(cfg_args, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"cmake configure failed:\n{r.stderr}", file=sys.stderr)
        return None

    print("Building X-GRMSD (this downloads and compiles the source; may take a minute) …")
    r = subprocess.run(
        [cmake, "--build", str(build_dir), "--target", "xgrmsd"],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        print(f"cmake build failed:\n{r.stderr}", file=sys.stderr)
        return None

    # INSTALL_COMMAND copies binary to CMAKE_BINARY_DIR (= build_dir)
    exe_suffix = ".exe" if sys.platform == "win32" else ""
    binary = build_dir / f"GRMSD{exe_suffix}"
    if binary.exists():
        print(f"X-GRMSD built: {binary}")
        return str(binary)

    # Fallback: search the entire build tree
    for candidate in build_dir.rglob(f"GRMSD{exe_suffix}"):
        print(f"X-GRMSD found: {candidate}")
        return str(candidate)

    print("GRMSD binary not found after build", file=sys.stderr)
    return None


def grmsd_rmsd(xyz_ref, xyz_trg, workdir: Path, grmsd_bin: str) -> float:
    n = len(xyz_ref)
    model = workdir / "model.csv"
    data  = workdir / "data.csv"
    out   = workdir / "out"
    for path, xyz in [(model, xyz_ref), (data, xyz_trg)]:
        with open(path, "w") as f:
            f.write(f"{n},{n},{n},{n}")
            for x, y, z in xyz:
                f.write(f"\n{x:.10f},{y:.10f},{z:.10f},0")
    subprocess.run(
        [grmsd_bin, "IsometryOpt", "model.csv", "data.csv", "0", "out", "0"],
        cwd=workdir, capture_output=True, check=False,
    )
    with open(out) as f:
        return float(f.readline().split()[-1])


# ── benchmark runners ────────────────────────────────────────────────────────

def run_rdkit_benchmark(n_mols: list[int], n_samples: int) -> dict:
    """Returns timing and accuracy stats for RDKit vs mobbrmsd."""
    from mobbrmsd.dataclass import molecules
    from mobbrmsd.mobbrmsd import mobbrmsd as MobbRMSD

    results = {}
    for n_mol in n_mols:
        print(f"  RDKit  n_mol={n_mol} ... ", end="", flush=True)
        rng = np.random.default_rng(seed=n_mol * 17 + 3)
        mrmsd = MobbRMSD(mols=molecules(n_apm=1, n_mol=n_mol))
        mol_ref = _make_dummy_mol(n_mol)
        mol_trg = _make_dummy_mol(n_mol)

        xs = [rng.standard_normal((n_mol, 3)) for _ in range(n_samples)]
        ys = [rng.standard_normal((n_mol, 3)) for _ in range(n_samples)]

        rdkit_rmsds, mobb_rmsds = [], []
        t_rd_list, t_mb_list = [], []

        for x, y in zip(xs, ys):
            t0 = time.perf_counter()
            r_rd = rdkit_best_rmsd(mol_ref, mol_trg, x, y)
            t_rd_list.append((time.perf_counter() - t0) * 1e6)

            t0 = time.perf_counter()
            r_mb = mrmsd.rmsd(x, y)
            t_mb_list.append((time.perf_counter() - t0) * 1e6)

            rdkit_rmsds.append(r_rd)
            mobb_rmsds.append(r_mb)

        results[n_mol] = {
            "t_mobb_median": float(np.median(t_mb_list)),
            "t_mobb_q25":    float(np.percentile(t_mb_list, 25)),
            "t_mobb_q75":    float(np.percentile(t_mb_list, 75)),
            "t_rdkit_median": float(np.median(t_rd_list)),
            "t_rdkit_q25":    float(np.percentile(t_rd_list, 25)),
            "t_rdkit_q75":    float(np.percentile(t_rd_list, 75)),
            "max_rmsd_diff":  float(np.max(np.abs(np.array(rdkit_rmsds) - np.array(mobb_rmsds)))),
            "rdkit_rmsds":    [float(v) for v in rdkit_rmsds[:200]],
            "mobb_rmsds":     [float(v) for v in mobb_rmsds[:200]],
        }
        ratio = results[n_mol]["t_rdkit_median"] / max(results[n_mol]["t_mobb_median"], 1e-9)
        print(f"mobb={results[n_mol]['t_mobb_median']:.1f} µs  rdkit={results[n_mol]['t_rdkit_median']:.1f} µs  ratio={ratio:.1f}x")
    return results


def run_grmsd_benchmark(n_mols: list[int], n_samples: int, grmsd_bin: str, workdir: Path) -> dict:
    """Returns timing and accuracy stats for GRMSD vs mobbrmsd.

    All times are stored in µs for consistency with the RDKit results.
    GRMSD subprocess overhead is large (~10-1000 ms), so ``t_grmsd_*``
    values are typically 10,000–1,000,000 µs.
    """
    from mobbrmsd.dataclass import molecules
    from mobbrmsd.mobbrmsd import mobbrmsd as MobbRMSD

    results = {}
    for n_mol in n_mols:
        print(f"  GRMSD  n_mol={n_mol} ... ", end="", flush=True)
        rng = np.random.default_rng(seed=n_mol * 31 + 7)
        mrmsd = MobbRMSD(mols=molecules(n_apm=1, n_mol=n_mol))

        xs = [rng.standard_normal((n_mol, 3)) for _ in range(n_samples)]
        ys = [rng.standard_normal((n_mol, 3)) for _ in range(n_samples)]

        grmsd_rmsds, mobb_rmsds = [], []
        t_gr_list, t_mb_list = [], []

        for x, y in zip(xs, ys):
            t0 = time.perf_counter()
            r_gr = grmsd_rmsd(x, y, workdir, grmsd_bin)
            t_gr_list.append((time.perf_counter() - t0) * 1e6)   # µs

            t0 = time.perf_counter()
            r_mb = mrmsd.rmsd(x, y)
            t_mb_list.append((time.perf_counter() - t0) * 1e6)   # µs

            grmsd_rmsds.append(r_gr)
            mobb_rmsds.append(r_mb)

        results[n_mol] = {
            "t_mobb_median":  float(np.median(t_mb_list)),
            "t_mobb_q25":     float(np.percentile(t_mb_list, 25)),
            "t_mobb_q75":     float(np.percentile(t_mb_list, 75)),
            "t_grmsd_median": float(np.median(t_gr_list)),
            "t_grmsd_q25":    float(np.percentile(t_gr_list, 25)),
            "t_grmsd_q75":    float(np.percentile(t_gr_list, 75)),
            "max_rmsd_diff":  float(np.max(np.abs(np.array(grmsd_rmsds) - np.array(mobb_rmsds)))),
            "grmsd_rmsds":    [float(v) for v in grmsd_rmsds[:50]],
            "mobb_rmsds":     [float(v) for v in mobb_rmsds[:50]],
        }
        t_gr_ms = results[n_mol]["t_grmsd_median"] / 1e3
        print(f"mobb={results[n_mol]['t_mobb_median']:.1f} µs  grmsd={t_gr_ms:.1f} ms (incl. subprocess)")
    return results


# ── plotting ─────────────────────────────────────────────────────────────────

def _errbar(ax, xs, medians, q25s, q75s, *, color, marker, label):
    errs = [
        [m - q25 for m, q25 in zip(medians, q25s)],
        [q75 - m for m, q75 in zip(medians, q75s)],
    ]
    ax.errorbar(
        xs, medians, yerr=errs,
        color=color, marker=marker, ms=5, mew=1.2, mfc="white",
        lw=1.0, elinewidth=1.0, capsize=3, label=label,
    )


def plot_timing(n_mols_rd, rd_data, n_mols_gr, gr_data, *, out_path: Path) -> None:
    """Timing comparison: median ± IQR per pair vs n_mol (log scale)."""
    fig, ax = plt.subplots(figsize=(5.5, 4.0))

    # mobbrmsd: use rd_data where available (200 samples, more precise), gr_data for larger N
    all_n = sorted(set(n_mols_rd) | set(n_mols_gr))
    mobb_n, mobb_med, mobb_q25, mobb_q75 = [], [], [], []
    for n in all_n:
        src = rd_data if (rd_data and n in rd_data) else (gr_data if (gr_data and n in gr_data) else None)
        if src:
            mobb_n.append(n)
            mobb_med.append(src[n]["t_mobb_median"])
            mobb_q25.append(src[n]["t_mobb_q25"])
            mobb_q75.append(src[n]["t_mobb_q75"])
    if mobb_n:
        _errbar(ax, mobb_n, mobb_med, mobb_q25, mobb_q75, color=C_MOBB, marker="o", label="mobbrmsd (B&B)")

    if rd_data:
        rds  = [rd_data[n]["t_rdkit_median"] for n in n_mols_rd]
        q25r = [rd_data[n]["t_rdkit_q25"]    for n in n_mols_rd]
        q75r = [rd_data[n]["t_rdkit_q75"]    for n in n_mols_rd]
        _errbar(ax, n_mols_rd, rds, q25r, q75r, color=C_RD, marker="s", label="RDKit GetBestRMS")

    if gr_data:
        grs  = [gr_data[n]["t_grmsd_median"] for n in n_mols_gr]
        q25g = [gr_data[n]["t_grmsd_q25"]    for n in n_mols_gr]
        q75g = [gr_data[n]["t_grmsd_q75"]    for n in n_mols_gr]
        _errbar(ax, n_mols_gr, grs, q25g, q75g, color=C_GR, marker="^", label="X-GRMSD Hungarian (incl. subprocess)")

    ax.set_yscale("log")
    ax.set_xlabel("Number of molecules $N$")
    ax.set_ylabel("Time per pair (µs)  [log scale]")
    ax.set_xticks(all_n)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="y", lw=0.4, alpha=0.5)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out_path.relative_to(REPO_ROOT)}")


def plot_speedup(n_mols, data, *, ref_key: str, ref_label: str, out_path: Path) -> None:
    """Speedup ratio (reference / mobbrmsd) vs n_mol."""
    fig, ax = plt.subplots(figsize=(4.5, 3.0))

    ratios = [data[n][ref_key] / max(data[n]["t_mobb_median"], 1e-9) for n in n_mols]
    bar_color = C_RD if "rdkit" in ref_key else C_GR
    ax.bar(n_mols, ratios, color=bar_color, alpha=0.8, width=0.6, edgecolor="none")
    ax.axhline(1, color=C_GRAY, lw=0.8, ls="--")
    ax.set_xlabel("Number of molecules $N$")
    ax.set_ylabel(f"Speedup: {ref_label} / mobbrmsd  (×)")
    ax.set_xticks(n_mols)
    ax.grid(axis="y", lw=0.4, alpha=0.5)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out_path.relative_to(REPO_ROOT)}")


def plot_accuracy(n_mols, data, *, ref_key: str, ref_label: str, out_path: Path) -> None:
    """Scatter: reference RMSD vs mobbrmsd RMSD (should lie on y=x)."""
    n_panels = len(n_mols)
    ncols = min(4, n_panels)
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.5 * ncols, 2.5 * nrows), squeeze=False)

    for idx, n_mol in enumerate(n_mols):
        ax = axes[idx // ncols][idx % ncols]
        ref_rmsds  = data[n_mol][ref_key]
        mobb_rmsds = data[n_mol]["mobb_rmsds"]
        lim = max(max(ref_rmsds), max(mobb_rmsds)) * 1.05
        ax.plot([0, lim], [0, lim], color=C_GRAY, lw=0.8, ls="--")
        ax.scatter(ref_rmsds, mobb_rmsds, s=4, alpha=0.5, color=C_MOBB, edgecolors="none")
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_title(f"$N={n_mol}$", fontsize=8)
        ax.set_box_aspect(1)
        if idx % ncols == 0:
            ax.set_ylabel("mobbrmsd", fontsize=7)
        if idx // ncols == nrows - 1:
            ax.set_xlabel(ref_label, fontsize=7)

    # hide unused panels
    for idx in range(n_panels, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out_path.relative_to(REPO_ROOT)}")


# ── main ─────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n-samples",  type=int, default=200, metavar="N",
                   help="Number of random pairs per n_mol (default: 200)")
    p.add_argument("--n-mol-max",  type=int, default=8,   metavar="M",
                   help="Maximum n_mol for RDKit comparison (default: 8)")
    p.add_argument("--grmsd-bin",  default=None, metavar="PATH",
                   help="Path to GRMSD binary (overrides GRMSD_BIN env var)")
    p.add_argument("--build-grmsd", action="store_true",
                   help="Download and build X-GRMSD via cmake before benchmarking")
    p.add_argument("--grmsd-build-dir", type=Path,
                   default=REPO_ROOT / "build" / "grmsd", metavar="DIR",
                   help="cmake build directory for X-GRMSD (default: build/grmsd)")
    p.add_argument("--grmsd-n-mol-max", type=int, default=15, metavar="M",
                   help="Maximum n_mol for GRMSD comparison (default: 15)")
    p.add_argument("--grmsd-n-samples", type=int, default=20, metavar="N",
                   help="Samples for GRMSD (subprocess overhead; default: 20)")
    return p.parse_args()


def main():
    # Force UTF-8 output on Windows (avoids cp932 errors for µ, →, …)
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    if hasattr(sys.stderr, "reconfigure"):
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")

    args = parse_args()
    IMG_DIR.mkdir(parents=True, exist_ok=True)

    grmsd_bin = args.grmsd_bin or find_grmsd_bin()
    if args.build_grmsd and not grmsd_bin:
        grmsd_bin = build_grmsd(args.grmsd_build_dir)
    results = {}

    # ── RDKit ──
    try:
        import rdkit  # noqa: F401
    except ImportError:
        print("rdkit not found — skipping RDKit benchmark", file=sys.stderr)
        rd_data, n_mols_rd = None, []
    else:
        n_mols_rd = list(range(2, args.n_mol_max + 1))
        print(f"Running RDKit benchmark (n_mol={n_mols_rd}, n_samples={args.n_samples}) …")
        rd_data = run_rdkit_benchmark(n_mols_rd, args.n_samples)
        results["rdkit"] = rd_data

    # ── GRMSD ──
    if grmsd_bin:
        import tempfile
        n_mols_gr = list(range(2, args.grmsd_n_mol_max + 1))
        print(f"Running GRMSD benchmark (n_mol={n_mols_gr}, n_samples={args.grmsd_n_samples}) …")
        with tempfile.TemporaryDirectory() as tmp:
            gr_data = run_grmsd_benchmark(n_mols_gr, args.grmsd_n_samples, grmsd_bin, Path(tmp))
        results["grmsd"] = gr_data
    else:
        gr_data, n_mols_gr = None, []
        print("GRMSD binary not found — skipping GRMSD benchmark (set GRMSD_BIN to enable)")

    if not results:
        print("Nothing to plot.", file=sys.stderr)
        sys.exit(1)

    # ── save raw data ──
    OUT_JSON.write_text(json.dumps(results, indent=2))
    print(f"  → {OUT_JSON.relative_to(REPO_ROOT)}")

    # ── plots ──
    print("Generating plots …")
    plot_timing(
        n_mols_rd, rd_data or {},
        n_mols_gr, gr_data or {},
        out_path=IMG_DIR / "timing.svg",
    )
    if rd_data:
        plot_speedup(
            n_mols_rd, rd_data,
            ref_key="t_rdkit_median", ref_label="RDKit",
            out_path=IMG_DIR / "speedup_rdkit.svg",
        )
        plot_accuracy(
            n_mols_rd, rd_data,
            ref_key="rdkit_rmsds", ref_label="RDKit GetBestRMS",
            out_path=IMG_DIR / "accuracy_rdkit.svg",
        )
    if gr_data:
        plot_speedup(
            n_mols_gr, gr_data,
            ref_key="t_grmsd_median", ref_label="GRMSD",
            out_path=IMG_DIR / "speedup_grmsd.svg",
        )
        plot_accuracy(
            n_mols_gr, gr_data,
            ref_key="grmsd_rmsds", ref_label="GRMSD IsometryOpt",
            out_path=IMG_DIR / "accuracy_grmsd.svg",
        )

    print("Done.")


if __name__ == "__main__":
    main()
