# -*- coding: utf-8 -*-

from typing import Union
import dataclasses


@dataclasses.dataclass(frozen=True)
class molecules:
    n_apm: int
    n_mol: int
    sym: Union[None, list] = None
    name: Union[None, str] = None
    """分子集合体のデータクラス

    Args:
        n_apm : number of atoms per molecule.
        n_mol : number of molecules.
        sym : intramolecular atomic permutation. (optional)
        name : name of chemical species. (optional)
    """


@dataclasses.dataclass(frozen=True)
class molecular_system:
    mols: list[molecules]
    name: Union[None, str] = None
    """一般分子系のデータクラス

    Args:
        mols : list of molecules.
        name : name of system. (optional)
    """


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
