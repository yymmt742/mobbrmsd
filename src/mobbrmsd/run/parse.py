import json
import numpy
import mdtraj
from ..dataclass import molecules, molecular_system
from pathlib import Path
from typing import Optional


def _coerce_json(v):
    """Parse v as JSON if it is a string; otherwise return as-is."""
    return json.loads(v) if isinstance(v, str) else v


# ---------------------------------------------------------------------------
# TopologyBuilder ↔ molecular_system conversions
# ---------------------------------------------------------------------------

def to_molecular_system(builder) -> molecular_system:
    """Convert a :class:`~mobbrmsd.run.topology_builder.TopologyBuilder` to
    a frozen :class:`~mobbrmsd.dataclass.molecular_system`.

    This is the implementation that :meth:`TopologyBuilder.build` delegates to.
    """
    if not builder._species:
        raise ValueError("No species defined — call add() at least once.")
    _BUILDER_ONLY = {"atom_names"}
    mols_list = [
        molecules(**{k: v for k, v in sp.items() if k not in _BUILDER_ONLY})
        for sp in builder._species
    ]
    atom_order = (
        numpy.array(builder._atom_order, dtype=numpy.intp)
        if builder._atom_order is not None
        else None
    )
    return molecular_system(mols=mols_list, name=builder._name, atom_order=atom_order)


def builder_from_topology(top, name: Optional[str] = None):
    """Create a :class:`~mobbrmsd.run.topology_builder.TopologyBuilder` from
    an MDTraj Topology using structural species ordering (n_mol asc, n_apm desc).
    """
    from .topology_builder import TopologyBuilder
    mols = molecular_system.from_topology(top, name=name)
    return TopologyBuilder.from_molecular_system(mols)


def builder_from_mdtraj(top, coords=None, name: Optional[str] = None):
    """Create a :class:`~mobbrmsd.run.topology_builder.TopologyBuilder` from
    an MDTraj Topology.

    When *coords* are provided, species are sorted by descending COM variance;
    otherwise structural ordering (n_mol asc, n_apm desc) is applied.
    """
    from .topology_builder import TopologyBuilder
    mols = molecular_system.from_mdtraj(top, coords=coords, name=name)
    return TopologyBuilder.from_molecular_system(mols)


def load(prms: dict) -> Optional[molecular_system]:
    """Build a molecular_system from a parameter dict.

    Accepts both ``"molecular_system"`` (preferred) and ``"molecules"`` (legacy)
    keys.  Values may be pre-parsed Python objects (dict / list) or raw JSON
    strings, e.g. when supplied via the command line.

    Within ``"molecular_system"``, recognised sub-keys are:

    * ``"molecules"`` — list of molecule-spec dicts (or a single dict).
    * ``"name"`` — optional system name.
    * ``"atom_order"`` — optional permutation index list.
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
        return None

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


class parser:
    def __init__(self, args, mols=None):
        prms = {**args.params}
        if args.inp is not None:
            try:
                with open(args.inp[0], "r") as f:
                    try:
                        prms = {**json.load(f), **prms}
                    except:
                        print(f"Warning : load json [{args.inp[0]}] is failed.")
                        pass
            except:
                print(f"Warning : open json [{args.inp[0]}] is failed.")
                pass

        def parse_coordinate(obj, dtype, top=None, mask=None):
            """Return (coords, top, mdt_traj_or_None).

            mdt_traj_or_None is an mdtraj.Trajectory when obj is a file path,
            sliced to the masked atoms when a mask is applied.  It is None when
            the coordinates come from a JSON literal or list.
            """
            if obj is None:
                return None, top, None

            if isinstance(obj, str):
                try:  # Assume obj is path and load with MDtraj.
                    mdt = mdtraj.load(obj, top=top)
                    top = mdt.topology if isinstance(top, str) else top
                    if (mask is not None) & (top is not None):
                        mx = top.select(mask)
                        if dtype != numpy.float32:
                            ref = (
                                numpy.array(mdt.xyz[:, mx], dtype=dtype) * 10.0
                            )  # to angs.
                        else:
                            ref = mdt.xyz[:, mx] * 10.0  # to angs.
                        mdt_out = mdt.atom_slice(mx)
                    else:
                        if dtype != numpy.float32:
                            ref = numpy.array(mdt.xyz, dtype=dtype) * 10.0  # to angs.
                        else:
                            ref = mdt.xyz * 10.0  # to angs.
                        mdt_out = mdt
                    return ref, top, mdt_out
                except:
                    try:  # Assume obj is json string and load with json.loads.
                        ref = numpy.array(json.loads(obj), dtype=dtype)
                    except:
                        raise Exception(f"Coordinate load error in : {obj}")
            else:  # Assume obj is list.
                ref = numpy.array(obj, dtype=dtype)
            return ref, top, None

        self.dtype = (
            numpy.float32 if prms.get("precision") == "single" else numpy.float64
        )
        self.mask = prms.get("mask")
        self.refxyz, self.top, self.refmdt = parse_coordinate(
            prms.get("reference"),
            self.dtype,
            top=prms.get("topology"),
            mask=self.mask,
        )
        if self.refxyz is None:
            raise Exception(f"No coordinate information was found from the input.")

        self.trgxyz, _, self.trgmdt = parse_coordinate(
            prms.get("target"),
            self.dtype,
            top=self.top,
            mask=self.mask,
        )
        self.mols = mols if mols is not None else load(prms)
        if self.mols is None:
            self.mols = molecular_system([molecules(n_apm=self.refxyz.shape[1], n_mol=1)])
        self.prms = prms
        self.runtype = prms.get("runtype")
        if self.runtype is None:
            nref = self.refxyz.shape[0]
            ntrg = 0 if self.trgxyz is None else self.trgxyz.shape[0]

            if ntrg == 0:
                # Single file: symmetric batch_run → n_frame × n_frame matrix
                from .batch_run import run
            elif nref == 1 and ntrg == 1:
                # Two single-frame files: verbose single-pair run
                from .verbose_run import run
            else:
                # Multi-frame: cross batch_run → n_ref × n_trg matrix
                from .batch_run import run

            self.run = run
