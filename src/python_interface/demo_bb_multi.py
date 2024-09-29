from . import __version__
from . import _demo
from .coord_generator import coord_generator
from ._mobbrmsd import *
import sys
import numpy


class __demo__(_demo._demo):
    def __init__(self, **kwarg):
        super().__init__(title="Multi component system", **kwarg)

    def read_input(self):
        import itertools
        import math

        def msg():
            print(
                "    number of molecular symmetry must be less than factrial to atoms per molecule."
            )
            return False

        n_mols = []
        n_apms = []
        n_syms = []
        for i in range(2):
            n_mol = _demo.readinp(
                f"input number of molecules [{i}]",
                2 * i + 3,
                check=lambda n_mol: (n_mol > 0) if isinstance(n_mol, int) else False,
            )
            if n_mol > 8:
                if not self.yes_or_no(
                    "This parameter may take time to compute. May this be run?"
                ):
                    continue

            n_apm = _demo.readinp(
                f"input number of atoms per molecule [{i}]",
                i * 4 + 2,
                check=lambda n_apm: ((n_apm > 0) if isinstance(n_apm, int) else False),
            )

            n_sym = _demo.readinp(
                "input number of molecular symmetry",
                n_apm,
                check=lambda n_sym: (
                    (n_sym > 0) & (True if n_sym <= math.factorial(n_apm) else msg())
                    if isinstance(n_sym, int)
                    else False
                ),
            )

            sym = _demo.generate_sym_indices(n_apm_, int(n_sym))
            n_mols += [n_mol]
            n_apms += [n_mol]
            n_syms += [n_mol]

        return {"n_mols": n_mols, "n_apms": n_apms, "n_syms": n_syms}

    def demo(self, n_mols=[2, 4], n_apms=[8, 3], n_syms=[1, 2], a=0.5, b=1.0, **kwargs):

        def print_ret(ret, post="", end="\n", to_console: bool = False):
            ev, er, ub, lb, sd = (
                ret.n_eval,
                ret.eval_ratio,
                2 * ret.bounds[0] + ret.autocorr,
                2 * ret.bounds[1] + ret.autocorr,
                ret.rmsd,
            )
            if sys.stdout.isatty():
                print(
                    f"\r  {ev:12d} {er:12.8f}{ub:16.6f}{lb:16.6f}{sd:12.6f}  ",
                    post,
                    end=end,
                )
            else:
                if to_console:
                    return
                print(f"  {ev:12d} {er:12.8f}{ub:16.6f}{lb:16.6f}{sd:12.6f}")

        molecules = []
        n_mols_ = n_mols.split(",") if isinstance(n_mols, str) else n_mols
        n_apms_ = n_apms.split(",") if isinstance(n_apms, str) else n_apms
        n_syms_ = n_syms.split(",") if isinstance(n_syms, str) else n_syms
        for n_mol, n_apm, n_sym in zip(n_mols_, n_apms_, n_syms_):
            n_mol_ = int(n_mol)
            n_apm_ = int(n_apm)
            n_sym_ = int(n_sym)
            print(n_mol, n_apm, n_sym)
            sym = _demo.generate_sym_indices(n_apm_, n_sym_)
            molecules += [DataclassMolecule(n_apm=n_apm_, n_mol=n_mol_, sym=sym)]
        a_ = float(a)
        b_ = float(b)

        cogen = coord_generator()
        x = numpy.vstack(
            [
                cogen.generate(mol.n_apm, mol.n_mol, a_, b_, dtype=self.prec).reshape(
                    [-1, 3]
                )
                for mol in molecules
            ]
        )
        y = numpy.vstack(
            [
                cogen.generate(mol.n_apm, mol.n_mol, a_, b_, dtype=self.prec).reshape(
                    [-1, 3]
                )
                for mol in molecules
            ]
        )
        x -= numpy.mean(x, 0)
        y -= numpy.mean(y, 0)
        z = y.copy()

        _demo.print_system(
            molecules,
            title="Demonstration of mobbRMSD with multi component system",
        )
        mrmsd = mobbrmsd(molecules=molecules)

        ub, lb = numpy.inf, -numpy.inf
        ret = mrmsd.run(x, y, maxeval=0)

        sep1 = "  ------------------------------------------------------------------------------"
        sep2 = "  ---------------------------------------|--------|-------------------|---------"

        print("        N_eval   Eval_ratio      Upperbound      Lowerbound      RMSD")
        print(sep1)
        i = 0
        xtra = ["__-¯¯", "-__-¯", "¯-__-", "¯¯-__", "-¯¯-_", "_-¯¯-"]
        erace = "     "
        while not ret.is_finished:
            print_ret(ret, post=xtra[int(i / 5000) % 6], end="", to_console=True)
            if ub > ret.bounds[0] or lb < ret.bounds[1]:
                print_ret(ret, post=erace)
            ub, lb = ret.bounds[0], ret.bounds[1]
            ret.restart(maxeval=0)
            i += 1
        ret.restart(maxeval=0, get_rotation=True)  # get rotation matrix
        print()

        z = ret.rotate_y(y)
        print(sep1)
        print("      -- Final results --")
        print("        N_eval   Eval_ratio      Upperbound      Lowerbound      RMSD")
        print_ret(ret)
        print(sep1, "\n")
        print(sep1)
        print(
            "       Reference     | Target (original) |  disp. |  Target (rotate)  |  disp."
        )
        d1, d2 = 0.0, 0.0
        for xi, yi, zi in zip(x, y, z):
            d1_, d2_ = numpy.sum(numpy.power(xi - yi, 2)), numpy.sum(
                numpy.power(xi - zi, 2)
            )
            d1 += d1_
            d2 += d2_
            print(
                f"  {xi[0]:6.2f}{xi[1]:6.2f}{xi[2]:6.2f}",
                f"|{yi[0]:6.2f}{yi[1]:6.2f}{yi[2]:6.2f} |{d1_:7.2f}",
                f"|{zi[0]:6.2f}{zi[1]:6.2f}{zi[2]:6.2f} |{d2_:7.2f}",
            )
        d1 /= x.shape[0]
        d2 /= x.shape[0]
        print(sep2)
        print(
            f"          mean squared deviation         |{d1:7.2f} |                   |{d2:7.2f}"
        )
        d1, d2 = numpy.sqrt(d1), numpy.sqrt(d2)
        print(
            f"        root mean squared deviation      |{d1:7.2f} |                   |{d2:7.2f}"
        )
        print(sep1)
