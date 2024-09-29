from . import __version__
from . import coord_generator
from . import mobbrmsd
from .coord_generator import coord_generator
from ._mobbrmsd import *
from . import _demo
import sys
import numpy
import pprint
import networkx
import matplotlib.pyplot as plt


class __demo__(_demo._demo):
    def __init__(self, **kwarg):
        super().__init__(title="batch run", **kwarg)

    def read_input(self):
        while True:
            n_mol = _demo.readinp(
                "input number of molecules",
                6,
                check=lambda n_mol: (n_mol > 0) if isinstance(n_mol, int) else False,
            )

            n_reference = _demo.readinp(
                "input number of referece structures",
                5,
                check=lambda n_reference: (
                    (n_reference > 0) if isinstance(n_reference, int) else False
                ),
            )

            n_target = _demo.readinp(
                "input number of target structures",
                10,
                check=lambda n_target: (
                    (n_target > 0) if isinstance(n_target, int) else False
                ),
            )

            if n_mol > 8 or (n_target * n_reference) > 1000:
                if not self.yes_or_no(
                    "This parameter may take time to compute. May this be run ?"
                ):
                    continue

            break
        print()
        return {"n_mol": n_mol, "n_reference": n_reference, "n_target": n_target}

    def demo(
        self,
        n_apm=3,
        n_mol=6,
        n_reference=5,
        n_target=10,
        n_sym=2,
        a=0.5,
        b=1.0,
        r=0.9,
        **kwarg,
    ):
        n_mol_ = int(n_mol)
        n_apm_ = int(n_apm)
        sym = _demo.generate_sym_indices(n_apm_, int(n_sym))
        n_reference_ = int(n_reference)
        n_target_ = int(n_target)
        a_ = float(a)
        b_ = float(b)
        r_ = float(r)
        cogen = coord_generator()

        x = numpy.array(
            [
                cogen.generate(n_apm_, n_mol_, a_, b_, dtype=self.prec).reshape(
                    [n_apm_ * n_mol_, 3]
                )
                for i in range(n_reference_)
            ]
        )
        for i in range(n_reference_ - 1):
            x[i + 1] = (1.0 - r_) * x[i + 1] + r_ * x[i]
        y = numpy.array(
            [
                cogen.generate(n_apm_, n_mol_, a_, b_, dtype=self.prec).reshape(
                    [n_apm_ * n_mol_, 3]
                )
                for i in range(n_target_)
            ]
        )
        for i in range(n_target_ - 1):
            y[i + 1] = (1.0 - r_) * y[i + 1] + r_ * y[i]

        molecules = DataclassMolecule(n_apm=n_apm_, n_mol=n_mol_, sym=sym)
        _demo.print_system(
            [molecules], title="Demonstration of batch mobbrmsd run (with OpenMP)"
        )

        mrmsd = mobbrmsd(molecules=molecules)
        rmsds = mrmsd.batch_run(x, y)
        del mrmsd

        sep1 = "  ------------------------------------------------------------------------------"
        print(sep1)
        print("       i   j         RMSD")
        print(sep1)
        for i, ri in enumerate(rmsds):
            for j, rij in enumerate(ri):
                print(f"    {i+1:4d}{j+1:4d}{rij:16.9f}")
            print()
        print(sep1)
        return {"mat": rmsds}

    def after(self, mat):

        if self.yes_or_no("Show graph ? (Open matplotlib window)"):
            plt.imshow(mat)
            plt.colorbar()
            plt.xlabel("target")
            plt.ylabel("reference")
            plt.show()
            plt.clf()
