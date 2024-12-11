import json
import numpy
import mdtraj
from ..dataclass import molecules, molecular_system
from pathlib import Path
from typing import Union


def load(prms: dict) -> Union[molecules, molecular_system, None]:

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
        return

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


class parser:
    def __init__(self, args):
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

            if obj is None:
                return None, top

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
                    else:
                        if dtype != numpy.float32:
                            ref = numpy.array(mdt.xyz, dtype=dtype) * 10.0  # to angs.
                        else:
                            ref = mdt.xyz * 10.0  # to angs.
                except:
                    try:  # Assume obj is json string and load with json.loads.
                        ref = numpy.array(json.loads(obj), dtype=dtype)
                    except:
                        raise Exception(f"Coordinate load error in : {obj}")
            else:  # Assume obj is list.
                ref = numpy.array(obj, dtype=dtype)
            return ref, top

        self.dtype = (
            numpy.float32 if prms.get("precision") == "single" else numpy.float64
        )
        self.mask = prms.get("mask")
        self.refxyz, self.top = parse_coordinate(
            prms.get("reference"),
            self.dtype,
            top=prms.get("topology"),
            mask=self.mask,
        )
        if self.refxyz is None:
            raise Exception(f"No coordinate information was found from the input.")

        self.trgxyz, _ = parse_coordinate(
            prms.get("target"),
            self.dtype,
            top=self.top,
            mask=self.mask,
        )
        self.mols = load(prms)
        if self.mols is None:
            self.mols = molecules(n_apm=self.refxyz.shape[1], n_mol=1)
        self.prms = prms
        self.runtype = prms.get("runtype")
        if self.runtype is None:
            nref = self.refxyz.shape[0]
            ntrg = 1 if self.trgxyz is None else self.trgxyz.shape[0]

            if nref * ntrg == 1:
                self.runtype = "verbose"
            else:
                self.runtype = "batch"
