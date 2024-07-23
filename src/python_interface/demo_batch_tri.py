from . import __version__
from . import coord_generator
from . import mobbrmsd
from ._mobbrmsd import *
import sys
import numpy
import pprint
import networkx
import matplotlib.pyplot as plt

title = "mobbrmsd triu matrix batch run"


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


def read_input() -> tuple:
    while True:
        while True:
            inp = input("    input number of molecules (default : 6)   >> ")
            if inp == "":
                inp = "6"
            if inp[0] == "q" or inp[0] == "Q":
                exit()
            elif inp.isdigit:
                n_mol = int(inp)
            if n_mol > 0:
                break

        while True:
            inp = input("    input number of structures (default : 10) >> ")
            if inp == "":
                inp = "10"
            if inp[0] == "q" or inp[0] == "Q":
                exit()
            elif inp.isdigit:
                n_target = int(inp)
            if n_target > 1:
                break

        if n_mol > 8 or n_target > 30:
            while True:
                inp = input(
                    "    This parameter may take time to compute. May this be run ? [Y/n] > "
                )
                if inp == "y" or inp == "Y" or inp == "n" or inp == "N":
                    break
            if inp == "n" or inp == "N":
                continue

        break

    print()
    return {"n_mol": n_mol, "n_target": n_target}


def main(n_apm=3, n_mol=6, n_target=10, sym=((1, 2, 0), (2, 0, 1)), a=0.5, b=1.0):
    cogen = coord_generator()
    x = numpy.array(
        [
            cogen.generate(n_apm, n_mol, a, b).reshape([n_apm * n_mol, 3])
            for i in range(n_target)
        ]
    )
    for i in range(n_target - 1):
        # x[i + 1] = i / n_target * x[-1] + (n_target - i) / n_target * x[0]
        x[i + 1] = 0.001 * x[i + 1] + 0.999 * x[i]

    sep1 = "  ------------------------------------------------------------------------------"
    sep2 = "  ---------------------------------------|--------|-------------------|---------"

    print(sep1)
    print("        Demonstration of batch mobbrmsd run")
    print(sep1)
    print("      --System settings--")
    print(
        f"    Atoms per molecule  :{n_apm:6d}",
    )
    print(f"    Number of molecule  :{n_mol:6d}")
    print(f"    Number of structure :{n_target:6d}")

    pp = pprint.pformat(tuple([i for i in range(n_apm)]), width=50, compact=True).split(
        "\n"
    )
    print("    Molecular symmetry  :     1", pp[0])
    for i, l in enumerate(pp[1:]):
        print("                               ", l)
    for i, s in enumerate(sym):
        pp = pprint.pformat(s, width=50, compact=True).split("\n")
        print(f"                        :{i+2:6d}", pp[0])
        for l in pp[1:]:
            print("                               ", l)
    print()

    molecules = DataclassMolecule(n_apm=n_apm, n_mol=n_mol, sym=sym)
    mrmsd = mobbrmsd(molecules=molecules)
    rmsds = mrmsd.batch_run(x)
    del mrmsd

    print(sep1)
    print("     i   j         RMSD")
    print(sep1)
    for i, ri in enumerate(rmsds):
        for j, rij in enumerate(ri):
            print(f"    {i:8d}{j:8d}{rij:16.9f}")
        print()
    print(sep1)
    return rmsds


def show_graph(mat):

    while True:
        inp = input("  Show graph ? (Open matplotlib window) ['y'es, 'n'o, 's'ave] >> ")
        if inp == "":
            continue
        if inp[0] == "n" or inp[0] == "N":
            print()
            return
        elif inp[0] == "q" or inp[0] == "Q":
            exit()
        break

    plt.imshow(mat)
    plt.show()
    plt.clf()
    print()


if __name__ == "__main__":
    main()
