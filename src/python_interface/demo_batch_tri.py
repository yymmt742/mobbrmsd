from . import __version__
from . import coord_generator
from . import mobbrmsd
from . import _demo
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
            if not _demo.yes_or_no(
                "This parameter may take time to compute. May this be run ?"
            ):
                continue

        break

    print()
    return {"n_mol": n_mol, "n_target": n_target}


class _demo_batch_tri(_demo._demo):
    def __init__(self, **kwarg):
        super().__init__(title="mobbrmsd triu matrix batch run", **kwarg)


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

    for i, ri in enumerate(rmsds):
        print(sep1)
        print("       i   j         RMSD")
        print(sep1)
        for j, rij in enumerate(ri):
            print(f"    {i+1:4d}{j+1:4d}{rij:16.9f}")
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
