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
            f"\r  {ev:12d} {er:12.8f}{ub:16.6f}{lb:16.6f}{sd:12.6f}  ",
            post,
            end=end,
        )
    else:
        if to_console:
            return
        print(f"  {ev:12d} {er:12.8f}{ub:16.6f}{lb:16.6f}{sd:12.6f}")


def main(n_apm=3, n_mol=6, n_target=10, sym=[[1, 2, 0], [2, 0, 1]], a=0.5, b=1.0):
    cogen = coord_generator()
    x = numpy.array(
        [
            cogen.generate(n_apm, n_mol, a, b).reshape([n_apm * n_mol, 3])
            for i in range(n_target)
        ]
    )

    sep1 = "  ------------------------------------------------------------------------------"
    sep2 = "  ---------------------------------------|--------|-------------------|---------"

    print(sep1)
    print("                  Molecular-oriented RMSD for Branch-and-bound")
    print(sep1)
    print("      --System settings--")
    print(
        f"    Atoms per molecule :{n_apm:14d}  | Molecular symmetry :   0 ",
        [i for i in range(len(sym[0]))],
    )
    print(f"    Number of molecule :{n_mol:14d}  |                        1 ", sym[0])
    for i in range(len(sym) - 1):
        print(
            f"                                        |                 {i+2:8d} ",
            [si for si in sym[i + 1]],
        )
    print()

    mrmsd = mobbrmsd()
    mrmsd.add_molecule(n_apm, n_mol, sym)

    edges, weights, states =mrmsd.min_span_tree(x)
    print(edges)
    print(weights)


if __name__ == "__main__":
    main()
