# -*- coding: utf-8 -*-

from typing import Union, Optional, List
import dataclasses


@dataclasses.dataclass(frozen=True)
class molecules:
    """
    Data Class for Identical Molecular Aggregates.

    Parameters:
        n_apm (int): number of atoms per molecule.
        n_mol (int): number of molecules.
        sym (Optional[List[int]]): intramolecular atomic permutation.
        name (Optional[str]): name of chemical species.
    """

    n_apm: int
    n_mol: int
    sym: Optional[List[int]] = None
    name: Optional[str] = None


@dataclasses.dataclass(frozen=True)
class molecular_system:
    """Data Classes for General Molecular Systems.

    Parameters:
        mols (List[molecules]): list of molecules.
        name (Optional[str]): name of system.
    """

    mols: List[molecules]
    name: Optional[str] = None


def load(prms: dict):
    """load

    Args:
        prms : dict. The main keys are follows:\n
               "molecular_system": list[molecules] or str (json string) or dict\n
               "molecules": dict
    """
    import json

    sys = prms.get("molecular_system")
    if sys is None:
        mols = prms.get("molecules")
        name = None
    else:
        if isinstance(sys, str):
            jsys = json.loads(sys)
            mols = jsys.get("mols")
            name = jsys.get("name")
        elif isinstance(sys, dict):
            mols = sys.get("mols")
            name = sys.get("name")
        else:
            mols = sys
            name = None

    if mols is None:
        raise IOError

    if isinstance(mols, str):
        jmols = json.loads(mols)
        if isinstance(jmols, dict):
            return molecules(**jmols)
        else:
            return molecular_system(mols=[molecules(**mol) for mol in jmols], name=name)
    else:
        if isinstance(mols, dict):
            return molecular_system(mols=[molecules(**mols)], name=name)
        else:
            return molecular_system(mols=[molecules(**mol) for mol in mols], name=name)
