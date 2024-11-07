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
    name: str = ""


##
# @class molecular_system
# @brief 一般の分子系
# @details mols : list of molecules.
#          name : name of system. (optional)


@dataclasses.dataclass(frozen=True)
class molecular_system:
    mols: list[molecules]
    name: str = ""
