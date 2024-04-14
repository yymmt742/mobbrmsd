from . import __version__
from . import coord_generator
from . import mobbrmsd
import sys
import numpy


def print_ret(ret, post="", end="\n", to_console: bool = False):
    ev, er, ub, lb, sd = (
        ret.n_eval,
        ret.eval_ratio,
        ret.bounds[0],
        ret.bounds[1],
        ret.rmsd,
    )
    if sys.stdout.isatty():
        print(
            f"\r  {ev:12d} {er:12.8f}{ub:16.6f}{lb:16.6f}{sd:12.6f}       ",
            post,
            end=end,
        )
    else:
        if to_console:
            return
        print(f"  {ev:12d} {er:12.8f}{ub:16.6f}{lb:16.6f}{sd:12.6f}")


def main(n_apm=3, n_mol=10, sym=[[1, 2, 0], [2, 0, 1]], a=0.5, b=1.0):
    cogen = coord_generator()
    x = cogen.generate(n_apm, n_mol, a, b).reshape([n_apm * n_mol, 3])
    y = cogen.generate(n_apm, n_mol, a, b).reshape([n_apm * n_mol, 3])
    x -= numpy.mean(x, 0)
    y -= numpy.mean(y, 0)
    z = y.copy()

    sep1 = "  ------------------------------------------------------------------------------"
    sep2 = "  ---------------------------------------|--------|-------------------|---------"

    print(sep1)
    print("                 Molecular-oriented RMSD for Branch-and-bound")
    print(sep1)
    print("      --System settings--")
    print(
        f"    Atoms per molecule :{n_apm:2d}    molecular symmetry : 0 ",
        [i for i in range(3)],
    )
    print(f"    Number of molecule :{n_mol:2d}                         1 ", sym[0])
    for i in range(len(sym) - 1):
        print(
            f"                                                  {i+2:2d} ",
            [si for si in sym[i + 1]],
        )
    print()

    mrmsd = mobbrmsd()
    mrmsd.add_molecule(n_apm, n_mol, sym)

    ub, lb = numpy.inf, 0.0
    ret = mrmsd.run(x, y, maxeval=0)

    print("        N_eval   Eval_ratio      Upperbound      Lowerbound      RMSD")
    print(sep1)
    i = 0
    xtra = ["|    ", " /   ", "  -  ", "   \\ ", "    |", "   \\ ", "  -  ", " /   "]
    erace = "     "
    while not ret.is_finished:
        print_ret(ret, post=xtra[int(i / 10000) % 8], end="", to_console=True)
        if ub > ret.bounds[0] or lb < ret.bounds[1]:
            print_ret(ret, post=erace)
        ub, lb = ret.bounds[0], ret.bounds[1]
        ret = mrmsd.restart(ret, maxeval=0)
        i += 1

    print_ret(ret, post=erace)
    ret = mrmsd.restart(ret, maxeval=0, Y=y)
    print(sep1)
    print("      -- Final results --")
    print("        N_eval   Eval_ratio      Upperbound      Lowerbound      RMSD")
    print_ret(ret)
    print(sep1, "\n")
    print(sep1)
    print(
        "       Reference     | Target (original) |  disp. |  Target (rotate)  |  disp."
    )
    for xi, yi, zi in zip(
        x.reshape([n_mol, n_apm, 3]),
        y.reshape([n_mol, n_apm, 3]),
        z.reshape([n_mol, n_apm, 3]),
    ):
        print(sep2)
        for xij, yij, zij in zip(xi, yi, zi):
            d1 = numpy.sum(numpy.power(xij - zij, 2))
            d2 = numpy.sum(numpy.power(xij - yij, 2))
            print(
                f"  {xij[0]:6.2f}{xij[1]:6.2f}{xij[2]:6.2f}",
                f"|{yij[0]:6.2f}{yij[1]:6.2f}{yij[2]:6.2f} |{d1:7.2f}",
                f"|{zij[0]:6.2f}{zij[1]:6.2f}{zij[2]:6.2f} |{d2:7.2f}",
            )
    d1, d2 = numpy.sum(numpy.power(x - z, 2)), numpy.sum(numpy.power(x - y, 2))
    print(sep2)
    print(
        f"             squared deviation           |{d1:7.2f} |                   |{d2:7.2f}"
    )
    d1 /= n_apm * n_mol
    d2 /= n_apm * n_mol
    print(
        f"          mean squared deviation         |{d1:7.2f} |                   |{d2:7.2f}"
    )
    d1, d2 = numpy.sqrt(d1), numpy.sqrt(d2)
    print(
        f"        root mean squared deviation      |{d1:7.2f} |                   |{d2:7.2f}"
    )
    print(sep1)


if __name__ == "__main__":
    main()
