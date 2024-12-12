import json
import numpy
import numpy.typing as npt
import mdtraj
from ..dataclass import molecules, molecular_system
from ..mobbrmsd import mobbrmsd
from pathlib import Path
from typing import Union


def run(mols, refxyz, trgxyz, **kwargs) -> None:
    mrmsd = mobbrmsd(mols=mols)
    if trgxyz is None:
        ret = mrmsd.batch_run(refxyz, **kwargs)
        for i, ri in enumerate(ret):
            for j, rij in enumerate(ri[i + 1 :]):
                print(f"{i+1:8d}{i+j+1:8d}{rij:24.9f}")
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
