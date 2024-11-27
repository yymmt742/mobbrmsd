import dataclasses


##
# @class molecules
# @brief 分子集合体のデータクラス
# @details n_apm : number of atoms per molecule.
#          n_mol : number of molecules.
#          sym : intramolecular atomic permutation. (optional)
#          name : name of chemical species. (optional)


@dataclasses.dataclass(frozen=True)
class molecules:
    n_apm: int
    n_mol: int
    sym: None | list = None
    name: None | str = None


##
# @class molecular_system
# @brief 一般分子系
# @details mols : list of molecules.
#          name : name of system. (optional)


@dataclasses.dataclass(frozen=True)
class molecular_system:
    mols: list[molecules]
    name: None | str = None


##
# @function load
# @brief load
# @details prms : json/molecules list.


def load(prms):
    import json

    sys = prms.get("molecular_system")
    if sys is None:
        mols = prms.get("molecules")
        name = None
    else:
        if isinstance(sys, str):
            try:
                jsys = json.loads(sys)
            except:
                raise IOError
            mols = jsys.get("mols")
            name = jsys.get("name")
        elif isinstance(sys, dict):
            mols = sys.get("mols")
            name = sys.get("name")
        else:
            mols = sys
            name = None

    if isinstance(mols, str):
        try:
            jmols = json.loads(mols)
        except:
            raise IOError
        if isinstance(jmols, dict):
            return molecules(**jmols)
        else:
            return molecular_system(mols=[molecules(**mol) for mol in jmols], name=name)
    else:
        return molecular_system(mols=[molecules(**mol) for mol in mols], name=name)
