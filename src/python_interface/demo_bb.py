from . import __version__
from . import _demo
from .coord_generator import coord_generator
from ._mobbrmsd import *
import sys
import numpy


class _demo_bb(_demo._demo):
    def __init__(self, **kwarg):
        super().__init__(title="mobbRMSD with random coordinate", **kwarg)

    def read_input(self):
        import itertools
        import math
        import pprint

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

        while True:
            n_mol = _demo.readinp(
                "input number of molecules",
                6,
                check=lambda n_mol: (n_mol > 0) if isinstance(n_mol, int) else False,
            )
            if n_mol > 8:
                if not _demo.yes_or_no(
                    "This parameter may take time to compute. May this be run ?"
                ):
                    continue

            n_apm = _demo.readinp(
                "input number of atoms per molecule",
                3,
                check=lambda n_apm: ((n_apm > 0) if isinstance(n_apm, int) else False),
            )

            def msg():
                print(
                    "    number of molecular symmetry must be less than factrial to atoms per molecule."
                )
                return False

            n_sym = _demo.readinp(
                "input number of molecular symmetry",
                n_apm,
                check=lambda n_sym: (
                    (n_sym > 0) & (True if n_sym <= math.factorial(n_apm) else msg())
                    if isinstance(n_sym, int)
                    else False
                ),
            )

            cost = math.factorial(n_mol) * n_sym**n_mol
            if cost > 10000000:
                if not _demo.yes_or_no(
                    "This parameter may take time to compute (cost is {cost:d}). May this be run ?"
                ):
                    continue

            break

        print()

        per = itertools.permutations(range(n_apm))
        next(per)
        sym = [next(per) for i in range(n_sym - 1)]

        return {"n_apm": n_apm, "n_mol": n_mol, "sym": sym}

    def demo(self, n_apm=3, n_mol=8, sym=((1, 2, 0), (2, 0, 1)), a=0.5, b=1.0):

        cogen = coord_generator()
        x = cogen.generate(n_apm, n_mol, a, b, self.prec).reshape([-1, 3])
        y = cogen.generate(n_apm, n_mol, a, b, self.prec).reshape([-1, 3])
        x -= numpy.mean(x, 0)
        y -= numpy.mean(y, 0)
        z = y.copy()

        sep1 = "  ------------------------------------------------------------------------------"
        sep2 = "  ---------------------------------------|--------|-------------------|---------"

        print(sep1)
        print("                Demonstration of mobbRMSD with random coordinates")
        print(sep1)
        print("      --System settings--")
        print(
            f"    Atoms per molecule :{n_apm:6d}",
        )
        print(f"    Number of molecule :{n_mol:6d}")
        print("    Molecular symmetry :     0 ")
        pp = pprint.pformat(tuple([i for i in range(n_apm)]), width=74, compact=True)
        for l in pp.split("\n"):
            print("      ", l)
        for i, s in enumerate(sym):
            print(f"                        {i+1:6d}")
            pp = pprint.pformat(s, width=74, compact=True)
            for l in pp.split("\n"):
                print("      ", l)
        print()

        molecules = DataclassMolecule(n_apm=n_apm, n_mol=n_mol, sym=sym)
        mrmsd = mobbrmsd(molecules=molecules)

        ub, lb = numpy.inf, -numpy.inf
        ret = mrmsd.run(x, y, maxeval=0, get_rotation=True)

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
            ret.restart(maxeval=0, get_rotation=True)
            i += 1

        print_ret(ret, post=erace)
        y = ret.rotate_y(y)
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
        for xi, yi, zi in zip(
            x.reshape([n_mol, n_apm, 3]),
            y.reshape([n_mol, n_apm, 3]),
            z.reshape([n_mol, n_apm, 3]),
        ):
            print(sep2)
            for xij, yij, zij in zip(xi, yi, zi):
                d1_ = numpy.sum(numpy.power(xij - zij, 2))
                d2_ = numpy.sum(numpy.power(xij - yij, 2))
                d1 += d1_
                d2 += d2_
                print(
                    f"  {xij[0]:6.2f}{xij[1]:6.2f}{xij[2]:6.2f}",
                    f"|{yij[0]:6.2f}{yij[1]:6.2f}{yij[2]:6.2f} |{d1_:7.2f}",
                    f"|{zij[0]:6.2f}{zij[1]:6.2f}{zij[2]:6.2f} |{d2_:7.2f}",
                )
        print(sep2)
        d1 *= 1.0 / (n_mol * n_apm)
        d2 *= 1.0 / (n_mol * n_apm)
        print(
            f"          mean squared deviation         |{d1:7.2f} |                   |{d2:7.2f}"
        )
        d1, d2 = numpy.sqrt(d1), numpy.sqrt(d2)
        print(
            f"        root mean squared deviation      |{d1:7.2f} |                   |{d2:7.2f}"
        )
        print(sep1)
