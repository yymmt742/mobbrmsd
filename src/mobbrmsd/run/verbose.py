import json
import numpy
import numpy.typing as npt
import mdtraj
from ..dataclass import molecules, molecular_system
from ..mobbrmsd import mobbrmsd, mobbrmsd_result
from pathlib import Path
from typing import Union


def print_result(ret: mobbrmsd_result) -> None:
    print(f"  RMSD        :{ret.rmsd():24.9f}")
    print(f"  ----------")
    print(f"  AUTO CORR.  :{ret.autocorr():24.9f}")
    print(
        f"  UPPER BOUND :{ret.upperbound():24.9f}  (as RMSD) :{ret.upperbound_as_rmsd():24.9f}"
    )
    print(
        f"  LOWER BOUND :{ret.lowerbound():24.9f}  (as RMSD) :{ret.lowerbound_as_rmsd():24.9f}"
    )
    print(f"  N EVAL      :{ret.n_eval():24d}  (ratio)   :{ret.eval_ratio():24.18f}")
    if ret.is_finished():
        pass
    else:
        print(f"  !!! BRANCH-AND-BOUND is NOT FINISHED. !!!")

    print(f"  ----------")

    return


def run(mols, refxyz, trgxyz, **kwargs) -> None:
    mrmsd = mobbrmsd(mols=mols)
    if trgxyz is None:
        for i in range(refxyz.shape[0]):
            for j in range(i + 1, refxyz.shape[0]):
                print(f"Reference frame{i+1:8d} | {j+1:8d}")
                ret = mrmsd.run(refxyz[i], refxyz[j], **kwargs)
                print_result(ret)
    else:
        for i in range(refxyz.shape[0]):
            for j in range(trgxyz.shape[0]):
                ret = mrmsd.run(refxyz[i], trgxyz[j], **kwargs)
                print(f"Reference frame{i+1:8d} | Target    frame{j+1:8d}")
                print_result(ret)

    return
