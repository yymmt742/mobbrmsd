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

    def demo(
        self, n_mols=[2, 4], n_apms=[8, 3], n_syms=[1, 2], alpha=0.5, beta=1.0, **kwargs
    ):

        molecules = []
        n_mols_ = (
            [int(n_mol) for n_mol in n_mols.split(",")]
            if isinstance(n_mols, str)
            else n_mols
        )
        n_apms_ = (
            [int(n_apm) for n_apm in n_apms.split(",")]
            if isinstance(n_apms, str)
            else n_apms
        )
        n_syms_ = (
            [int(n_sym) for n_sym in n_syms.split(",")]
            if isinstance(n_syms, str)
            else n_syms
        )
        for n_mol, n_apm, n_sym in zip(n_mols_, n_apms_, n_syms_):
            sym = _demo.generate_sym_indices(n_apm, n_sym)
            molecules += [DataclassMolecule(n_apm=n_apm, n_mol=n_mol, sym=sym)]
        a_ = float(alpha)
        b_ = float(beta)

        cogen = coord_generator()
        x, y = cogen.generate_pair(
            n_apms_, n_mols_, alpha=a_, beta=b_, dtype=self.prec, remove_com=False
        )
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

        xtra = ["__-¯¯", "-__-¯", "¯-__-", "¯¯-__", "-¯¯-_", "_-¯¯-"]
        erace = "     "
        i = 0
        _demo.print_ret(ret, end="", header=True)
        while not ret.is_finished():
            _demo.print_ret(ret, post=xtra[int(i / 5000) % 6], end="", to_console=True)
            if ub > ret.upperbound() or lb < ret.lowerbound():
                _demo.print_ret(ret, post=erace)
            ub, lb = ret.upperbound(), ret.lowerbound()
            ret.restart(maxeval=0)
            i += 1
        ret.restart(maxeval=0, get_rotation=True)  # get rotation matrix
        _demo.print_ret(ret, post=erace, header=True, footer=True)

        z = ret.rotate_y(y)

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
