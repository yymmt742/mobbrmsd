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

title = "mobbrmsd batch run"


def read_input() -> tuple:
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
            if not _demo.yes_or_no(
                "This parameter may take time to compute. May this be run ?"
            ):
                continue

        break
    print()
    return {"n_mol": n_mol, "n_reference": n_reference, "n_target": n_target}


class _demo_batch(_demo._demo):
    def __init__(self, **kwarg):
        super().__init__(title="mobbrmsd batch run", **kwarg)

    def read_input(self):
        return read_input()

    def main(self, **kwarg):
        return main(**kwarg)

    def after(self, **kwarg):
        return show_graph(**kwarg)


def main(
    n_apm=3,
    n_mol=6,
    n_reference=5,
    n_target=10,
    sym=((1, 2, 0), (2, 0, 1)),
    a=0.5,
    b=1.0,
):
    cogen = coord_generator()
    x = numpy.array(
        [
            cogen.generate(n_apm, n_mol, a, b).reshape([n_apm * n_mol, 3])
            for i in range(n_reference)
        ]
    )
    for i in range(n_reference - 1):
        x[i + 1] = 0.01 * x[i + 1] + 0.99 * x[i]
    y = numpy.array(
        [
            cogen.generate(n_apm, n_mol, a, b).reshape([n_apm * n_mol, 3])
            for i in range(n_target)
        ]
    )
    for i in range(n_target - 1):
        y[i + 1] = 0.01 * y[i + 1] + 0.99 * y[i]

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
    rmsds = mrmsd.batch_run(x, y)
    del mrmsd

    print(sep1)
    print("     i   j         RMSD")
    print(sep1)
    for i, ri in enumerate(rmsds):
        for j, rij in enumerate(ri):
            print(f"    {i:8d}{i:8d}{rij:16.9f}")
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
