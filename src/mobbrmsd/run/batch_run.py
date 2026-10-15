import sys
import json
import numpy
import numpy.typing as npt
import mdtraj
from ..dataclass import molecules, molecular_system
from ..mobbrmsd import mobbrmsd
from pathlib import Path
from typing import Union


def _numbered_path(base, i, n):
    """Insert zero-padded frame index before the suffix; omit when n == 1."""
    if n <= 1:
        return base
    p = Path(base)
    width = max(3, len(str(n)))
    return str(p.parent / f"{p.stem}.{i + 1:0{width}d}{p.suffix}")


def _save_rotated_frames(frames, atom_order, topology, output, overwrite=False):
    """Stack *frames* (mobbrmsd atom order), apply inverse permutation, save."""
    if Path(output).exists() and not overwrite:
        print(
            f"Warning: {output} already exists; skipping. Use -O to overwrite.",
            file=sys.stderr,
        )
        return
    if atom_order is not None:
        inv = numpy.argsort(atom_order)
        xyz = numpy.stack([f[inv] for f in frames])
    else:
        xyz = numpy.stack(frames)
    mdtraj.Trajectory(xyz / 10.0, topology).save(output)


def run(mols, refxyz, trgxyz, output=None, trgmdt=None, refmdt=None, overwrite=False, **kwargs) -> None:
    mrmsd = mobbrmsd(mols=mols)

    if output is not None:
        # Output mode: Python loop calling run(rotate_y=True) per pair.
        if trgxyz is None:
            # Self-comparison: fit every frame to each reference frame.
            n_ref = refxyz.shape[0]
            src_top = refmdt.topology if refmdt is not None else None
            if src_top is None:
                print(
                    "Warning: -o requires coordinates loaded from a file; skipping output.",
                    file=sys.stderr,
                )
            for i in range(n_ref):
                frames = []
                for j in range(n_ref):
                    y = refxyz[j].copy()
                    ret = mrmsd.run(refxyz[i], y, rotate_y=True, **kwargs)
                    print(f"{i+1:8d}{j+1:8d}{ret.rmsd():24.9f}")
                    frames.append(y)
                print()
                if src_top is not None:
                    out_path = _numbered_path(output, i, n_ref)
                    _save_rotated_frames(frames, mrmsd._atom_order, src_top, out_path, overwrite=overwrite)
                    print(f"  Saved to {out_path}")
        else:
            # Cross: fit every target frame to each reference frame.
            n_ref = refxyz.shape[0]
            n_trg = trgxyz.shape[0]
            src_top = trgmdt.topology if trgmdt is not None else None
            if src_top is None:
                print(
                    "Warning: -o requires target coordinates loaded from a file; skipping output.",
                    file=sys.stderr,
                )
            for i in range(n_ref):
                frames = []
                for j in range(n_trg):
                    y = trgxyz[j].copy()
                    ret = mrmsd.run(refxyz[i], y, rotate_y=True, **kwargs)
                    print(f"{i+1:8d}{j+1:8d}{ret.rmsd():24.9f}")
                    frames.append(y)
                print()
                if src_top is not None:
                    out_path = _numbered_path(output, i, n_ref)
                    _save_rotated_frames(frames, mrmsd._atom_order, src_top, out_path, overwrite=overwrite)
                    print(f"  Saved to {out_path}")
        return

    # No output: use fast Fortran batch_run.
    if trgxyz is None:
        ret = mrmsd.batch_run(refxyz, **kwargs)
        for i, ri in enumerate(ret):
            for j, rij in enumerate(ri[i + 1 :]):
                print(f"{i+1:8d}{i+j+2:8d}{rij:24.9f}")
            print()
    else:
        ret = mrmsd.batch_run(
            refxyz,
            trgxyz,
            **kwargs,
        )
        for i, ri in enumerate(ret):
            for j, rij in enumerate(ri):
                print(f"{i+1:8d}{j+1:8d}{rij:24.9f}")
            print()

    return
