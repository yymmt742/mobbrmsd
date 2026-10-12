# -*- coding: utf-8 -*-

from typing import Union, Optional, List
import dataclasses
import numpy
import numpy.typing as npt


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


def _sort_species_structurally(
    mols_list: List[molecules],
    atom_order: npt.NDArray,
) -> tuple:
    """Sort species by descending n_apm then ascending n_mol (no coords needed).

    Placing species with more atoms per molecule first tends to establish
    tighter bounds earlier in branch-and-bound.  Among species with equal
    n_apm, fewer molecules (smaller permutation space) come first.
    """
    per_species = []
    pos = 0
    for m in mols_list:
        n = m.n_apm * m.n_mol
        per_species.append(atom_order[pos : pos + n])
        pos += n

    order = sorted(
        range(len(mols_list)),
        key=lambda i: (-mols_list[i].n_apm, mols_list[i].n_mol),
    )
    new_mols = [mols_list[i] for i in order]
    new_atom_order = numpy.concatenate([per_species[i] for i in order])
    return new_mols, new_atom_order


def _sort_species_by_variance(
    mols_list: List[molecules],
    atom_order: npt.NDArray,
    coords: npt.NDArray,
) -> tuple:
    """Sort species by descending COM variance. Returns (new_mols_list, new_atom_order).

    Putting the most spatially dispersed species first tends to improve
    branch-and-bound pruning efficiency.
    """
    if coords.ndim == 2:
        coords = coords[numpy.newaxis]  # [1, natom, d]

    d = coords.shape[-1]
    nframe = coords.shape[0]

    per_species: List[npt.NDArray] = []
    pos = 0
    for m in mols_list:
        n = m.n_apm * m.n_mol
        per_species.append(atom_order[pos : pos + n])
        pos += n

    variances: List[float] = []
    for m, indices in zip(mols_list, per_species):
        sp_coords = coords[:, indices, :]
        coms = sp_coords.reshape(nframe, m.n_mol, m.n_apm, d).mean(axis=2)
        variances.append(float(coms.var()))

    order = numpy.argsort(variances)[::-1]
    new_mols = [mols_list[i] for i in order]
    new_atom_order = numpy.concatenate([per_species[i] for i in order])
    return new_mols, new_atom_order


@dataclasses.dataclass(frozen=True)
class molecular_system:
    """Data Class for General Molecular Systems.

    Parameters:
        mols (List[molecules]): list of molecules.
        name (Optional[str]): name of system.
        atom_order (Optional[npt.NDArray]): permutation index array such that
            ``coords[atom_order]`` gives atoms in the order the core expects
            (species-grouped). Set automatically by :meth:`from_mdtraj` and
            :meth:`sorted_by_variance`.
    """

    mols: List[molecules]
    name: Optional[str] = None
    atom_order: Optional[npt.NDArray] = dataclasses.field(
        default=None, compare=False, hash=False, repr=False
    )

    @classmethod
    def _build_from_top(cls, top) -> tuple:
        """Parse an MDTraj Topology into (mols_list, atom_order) in residue order.

        Raises:
            ValueError: if residues sharing a name have different atom counts.
        """
        from collections import defaultdict

        species_residues: dict = defaultdict(list)
        species_order: List[str] = []
        for res in top.residues:
            atom_indices = [a.index for a in res.atoms]
            if res.name not in species_residues:
                species_order.append(res.name)
            species_residues[res.name].append(atom_indices)

        for sp_name in species_order:
            n_apm = len(species_residues[sp_name][0])
            if any(len(r) != n_apm for r in species_residues[sp_name]):
                raise ValueError(
                    f"Residues named '{sp_name}' have inconsistent atom counts"
                )

        mols_list = [
            molecules(
                n_apm=len(species_residues[sp][0]),
                n_mol=len(species_residues[sp]),
                name=sp,
            )
            for sp in species_order
        ]
        atom_order = numpy.array(
            [
                idx
                for sp in species_order
                for res in species_residues[sp]
                for idx in res
            ],
            dtype=numpy.intp,
        )
        return mols_list, atom_order

    @classmethod
    def from_topology(cls, top, name: Optional[str] = None) -> "molecular_system":
        """Create from an MDTraj Topology using structural species ordering.

        Groups atoms by residue name. Species are sorted by descending
        atoms-per-molecule (n_apm), then ascending molecule count (n_mol).
        This ordering needs no coordinate data and tends to improve
        branch-and-bound pruning by placing heavier, less-numerous species first.

        Args:
            top: ``mdtraj.Topology`` instance.
            name: optional system name.

        Returns:
            ``molecular_system`` with ``atom_order`` set.

        Raises:
            ValueError: if residues of the same name have different atom counts.
        """
        mols_list, atom_order = cls._build_from_top(top)
        mols_list, atom_order = _sort_species_structurally(mols_list, atom_order)
        return cls(mols=mols_list, name=name, atom_order=atom_order)

    @classmethod
    def from_mdtraj(cls, top, coords=None, name: Optional[str] = None) -> "molecular_system":
        """Create from an MDTraj Topology, grouping atoms by residue name.

        Each unique residue name becomes one molecular species. All residues of
        the same name must have the same number of atoms.

        Species ordering:

        * When *coords* are provided, species are sorted by descending COM
          variance (most spatially dispersed first) for better pruning.
        * When *coords* are absent, structural ordering is used instead:
          descending n_apm, then ascending n_mol.

        Args:
            top: ``mdtraj.Topology`` instance.
            coords: optional reference coordinates ``[natom, d]`` or
                ``[nframe, natom, d]``.
            name: optional system name.

        Returns:
            ``molecular_system`` with ``atom_order`` set.

        Raises:
            ValueError: if residues of the same name have different atom counts.
        """
        mols_list, atom_order = cls._build_from_top(top)

        if coords is not None:
            mols_list, atom_order = _sort_species_by_variance(
                mols_list, atom_order, numpy.asarray(coords)
            )
        else:
            mols_list, atom_order = _sort_species_structurally(mols_list, atom_order)

        return cls(mols=mols_list, name=name, atom_order=atom_order)

    def to_dict(self) -> dict:
        """Serialise to a dict compatible with the JSON input format.

        Returns a dict of the form::

            {
                "molecular_system": {
                    "molecules": [{"n_apm": ..., "n_mol": ..., ...}, ...],
                    "name": ...,        # omitted when None
                    "atom_order": [...] # omitted when None
                }
            }

        The returned dict can be passed to :func:`json.dumps` directly or
        used as input to :func:`load`.
        """
        mol_list = []
        for m in self.mols:
            entry: dict = {"n_apm": m.n_apm, "n_mol": m.n_mol}
            if m.sym is not None:
                entry["sym"] = m.sym
            if m.name is not None:
                entry["name"] = m.name
            mol_list.append(entry)

        inner: dict = {"molecules": mol_list}
        if self.name is not None:
            inner["name"] = self.name
        if self.atom_order is not None:
            inner["atom_order"] = self.atom_order.tolist()
        return {"molecular_system": inner}

    def sorted_by_variance(self, coords: npt.NDArray) -> "molecular_system":
        """Return a new molecular_system with species sorted by descending COM variance.

        Args:
            coords: reference coordinates ``[natom, d]`` or ``[nframe, natom, d]``.

        Returns:
            New ``molecular_system`` whose species list and ``atom_order`` are
            reordered so the most spatially dispersed species comes first.
        """
        coords = numpy.asarray(coords)
        if self.atom_order is not None:
            current_order = self.atom_order
        else:
            current_order = numpy.arange(sum(m.n_apm * m.n_mol for m in self.mols))

        new_mols, new_order = _sort_species_by_variance(
            list(self.mols), current_order, coords
        )
        return molecular_system(mols=new_mols, name=self.name, atom_order=new_order)


def _coerce_json(v):
    """Parse v as JSON if it is a string; otherwise return as-is."""
    import json

    return json.loads(v) if isinstance(v, str) else v


def load(prms: dict) -> "molecular_system":
    """Build a molecular_system from a parameter dict.

    Args:
        prms (dict): Keys recognised:

            ``"molecular_system"`` (preferred)
                A dict (or JSON string) with sub-keys:

                * ``"molecules"`` — list of molecule-spec dicts (or a single dict).
                * ``"name"`` — optional system name.
                * ``"atom_order"`` — optional permutation index list.

            ``"molecules"`` (legacy top-level key)
                A list of molecule-spec dicts, or a single dict.

            Values may be pre-parsed Python objects or raw JSON strings.
    """
    sys_val = _coerce_json(prms.get("molecular_system"))
    if isinstance(sys_val, dict):
        raw_mols = _coerce_json(sys_val.get("molecules"))
        name = sys_val.get("name")
        raw_atom_order = sys_val.get("atom_order")
    elif sys_val is not None:
        raw_mols = sys_val
        name = None
        raw_atom_order = None
    else:
        raw_mols = _coerce_json(prms.get("molecules"))
        name = None
        raw_atom_order = None

    if raw_mols is None:
        raise IOError("No molecular system information found in prms.")

    atom_order = (
        numpy.array(raw_atom_order, dtype=numpy.intp)
        if raw_atom_order is not None
        else None
    )

    if isinstance(raw_mols, dict):
        return molecular_system(
            mols=[molecules(**raw_mols)], name=name, atom_order=atom_order
        )
    return molecular_system(
        mols=[molecules(**m) for m in raw_mols], name=name, atom_order=atom_order
    )
