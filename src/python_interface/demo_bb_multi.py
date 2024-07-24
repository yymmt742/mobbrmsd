from . import __version__
from .coord_generator import coord_generator
from ._mobbrmsd import *
import sys
import pprint
import numpy

title = "mobbRMSD with multi component system"


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


def read_input() -> tuple:
    import itertools
    import math

    ret = []
    for i in range(2):
        while True:
            while True:
                inp = input(
                    f"    input number of molecules [{i}] (default : {2*i+3})          >> "
                )
                if inp == "":
                    n_mol = 2 * i + 3
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
                inp = input(
                    f"    input number of atoms per molecule [{i}] (default : {i*4+2}) >> "
                )
                if inp == "":
                    n_apm = i * 4 + 2
                    break
                if inp[0] == "q" or inp[0] == "Q":
                    exit()
                elif inp.isdigit():
                    n_apm = int(inp)
                    if n_apm > 0:
                        break

            while True:
                inp = input(
                    f"    input number of molecular symmetry [{i}] (default : {i+1}) >> "
                )
                if inp == "":
                    n_sym = i + 1
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
            if cost > 100000:
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
        ret += [DataclassMolecule(n_apm=n_apm, n_mol=n_mol, sym=sym)]

    return {"molecules": ret}


def main(molecules=[], a=0.5, b=1.0):
    cogen = coord_generator()
    x = numpy.vstack(
        [
            cogen.generate(mol.n_apm, mol.n_mol, a, b).reshape(
                [mol.n_apm * mol.n_mol, 3]
            )
            for mol in molecules
        ]
    )
    y = numpy.vstack(
        [
            cogen.generate(mol.n_apm, mol.n_mol, a, b).reshape(
                [mol.n_apm * mol.n_mol, 3]
            )
            for mol in molecules
        ]
    )

    sep1 = "  ------------------------------------------------------------------------------"
    sep2 = "  ---------------------------------------|--------|-------------------|---------"

    print(sep1)
    print("            Demonstration of mobbRMSD with multi component system")
    print(sep1)
    print("      --System settings--")
    for k, mol in enumerate(molecules):
        print(
            f"    Atoms per molecule [{k}] :{mol.n_apm:6d}",
        )
        print(f"    Number of molecule     :{mol.n_mol:6d}")
        print("    Molecular symmetry     :     0 ")
        pp = pprint.pformat(
            tuple([i for i in range(mol.n_apm)]), width=74, compact=True
        )
        for l in pp.split("\n"):
            print("      ", l)
        for i, s in enumerate(mol.sym):
            print(f"                            {i+1:6d}")
            pp = pprint.pformat(s, width=74, compact=True)
            for l in pp.split("\n"):
                print("      ", l)
    print()

    mrmsd = mobbrmsd(molecules=molecules)

    ub, lb = numpy.inf, -numpy.inf
    ret = mrmsd.run(x, y, maxeval=0)

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

    print_ret(ret, post=erace)
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


if __name__ == "__main__":
    main()
