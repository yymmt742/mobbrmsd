from . import _demo
from ._coord_generator import coord_generator
from ..dataclass import molecules, molecular_system
from ..mobbrmsd import mobbrmsd
import sys
import numpy
import pprint
import networkx
import matplotlib.pyplot as plt


class __demo(_demo._demo):
    def __init__(self, **kwarg):
        super().__init__(title="mobbrmsd triu matrix", **kwarg)

    def read_input(self):
        while True:
            n_mol = _demo.readinp(
                "input number of molecules",
                6,
                check=lambda n_mol: (n_mol > 0) if isinstance(n_mol, int) else False,
            )

            n_target = _demo.readinp(
                "input number of target structures",
                10,
                check=lambda n_target: (
                    (n_target > 0) if isinstance(n_target, int) else False
                ),
            )

            if n_mol > 8 or n_target > 30:
                if not self.yes_or_no(
                    "This parameter may take time to compute. May this be run ?"
                ):
                    continue

            break

        print()
        return {"n_mol": n_mol, "n_target": n_target}

    def demo(
        self,
        n_apm=3,
        n_mol=6,
        n_sym=2,
        n_target=10,
        a=0.5,
        b=1.0,
        r=0.9,
        **kwarg,
    ):
        def print_ret(i, j, ret):
            ev, rm, er, ub, lb, df = (
                ret.n_eval,
                ret.eval_ratio,
                ret.rmsd,
                ret.bounds[0],
                ret.bounds[1],
                ret.bounds[0] - ret.bounds[1],
            )

            pre = ""
            post = ""

            if ub > 1.0e10:
                print(
                    pre,
                    f" {j:4d}{i:4d}{ev:12d} {rm:12.6f}          +Infty{lb:16.6f}  +Infty",
                    post,
                )
            else:
                print(
                    pre,
                    f" {j:4d}{i:4d}{ev:12d} {rm:12.6f}{ub:16.6f}{lb:16.6f}{df:8.3f}",
                    post,
                )

        n_mol_ = int(n_mol)
        n_apm_ = int(n_apm)
        sym = _demo.generate_sym_indices(n_apm_, int(n_sym))
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
                for i in range(n_target_)
            ]
        )
        for i in range(n_target_ - 1):
            x[i + 1] = (1.0 - r_) * x[i + 1] + r_ * x[i]

        mols = molecules(n_apm=n_apm_, n_mol=n_mol_, sym=sym)
        _demo.print_system(
            molecular_system([mols]),
            title="Demonstration of batch mobbrmsd triangle run (with OpenMP)",
        )

        mrmsd = mobbrmsd(mols=mols)
        rmsds = mrmsd.batch_run(x)
        del mrmsd

        sep1 = "  ------------------------------------------------------------------------------"
        for i, ri in enumerate(rmsds):
            print(sep1)
            print("       i   j         RMSD")
            print(sep1)
            for j, rij in enumerate(ri):
                print(f"    {i+1:4d}{j+1:4d}{rij:16.9f}")
        print(sep1)
        return {"mat": rmsds}

    def after(self, mat=[[0.0]], path=None, **kwargs):
        if self.yes_or_no("Show graph ? (Open matplotlib window)"):
            plt.imshow(mat)
            plt.colorbar()
            plt.xlabel("target")
            plt.ylabel("reference")
            if path is None:
                plt.show()
            else:
                for p in path.split(","):
                    plt.savefig(p)
            plt.clf()
