from . import __version__
from . import coord_generator
from . import mobbrmsd
import sys
import pprint
import numpy

title = "Demonstration of mobbRMSD with random coordinates"


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
            f"\r  {ev:12d} {er:12.8f}{ub:16.6f}{lb:16.6f}{sd:12.6f}  ",
            post,
            end=end,
        )
    else:
        if to_console:
            return
        print(f"  {ev:12d} {er:12.8f}{ub:16.6f}{lb:16.6f}{sd:12.6f}")


def read_input() -> tuple:
    import itertools
    import math

    while True:
        while True:
            inp = input("    input number of molecules (default : 6)          >> ")
            if inp == "":
                n_mol = 6
                break
            if inp[0] == "q" or inp[0] == "Q":
                exit()
            elif inp.isdigit():
                n_mol = int(inp)
                if n_mol > 0:
                    break

        if n_mol > 8:
            while True:
                inp = input(
                    "    This parameter may take time to compute. May this be run ? [Y/n] > "
                )
                if inp == "":
                    continue
                if inp[0] == "q" or inp[0] == "Q":
                    exit()
                elif inp == "y" or inp == "Y" or inp == "n" or inp == "N":
                    break
            if inp == "n" or inp == "N":
                continue

        while True:
            inp = input("    input number of atoms per molecule (default : 3) >> ")
            if inp == "":
                n_apm = 3
                break
            if inp[0] == "q" or inp[0] == "Q":
                exit()
            elif inp.isdigit():
                n_apm = int(inp)
                if n_apm > 1:
                    break

        while True:
            inp = input("    input number of molecular symmetry (default : 3) >> ")
            if inp == "":
                n_sym = 3
                break
            if inp[0] == "q" or inp[0] == "Q":
                exit()
            elif inp.isdigit():
                n_sym = int(inp)
                if n_sym > 0:
                    break

            if n_sym > math.factorial(n_apm):
                print(
                    "    number of molecular symmetry must be less than factrial to atoms per molecule."
                )
                continue

        cost = math.factorial(n_mol) * n_sym**n_mol
        if cost > 10000000:
            while True:
                inp = input(
                    f"    This parameter may take time to compute. ({cost:10d}) May this be run ? [Y/n] > "
                )
                if inp == "":
                    continue
                if inp[0] == "q" or inp[0] == "Q":
                    exit()
                elif inp == "y" or inp == "Y" or inp == "n" or inp == "N":
                    break
            if inp == "n" or inp == "N":
                continue

        per = itertools.permutations(range(n_apm))
        next(per)
        sym = [next(per) for i in range(n_sym - 1)]

        break

    return {"n_apm": n_apm, "n_mol": n_mol, "sym": sym}


def main(n_apm=3, n_mol=8, sym=((1, 2, 0), (2, 0, 1)), a=0.5, b=1.0):
    cogen = coord_generator()
    x = cogen.generate(n_apm, n_mol, a, b).reshape([n_apm * n_mol, 3])
    y = cogen.generate(n_apm, n_mol, a, b).reshape([n_apm * n_mol, 3])
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

    mrmsd = mobbrmsd()
    mrmsd.add_molecule(n_apm, n_mol, sym)

    ub, lb = numpy.inf, 0.0
    ret = mrmsd.run(x, y, maxeval=0)
    mrmsd.clear()

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
    print(sep2)
    d1, d2 = numpy.mean(numpy.power(x - z, 2)), numpy.mean(numpy.power(x - y, 2))
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
