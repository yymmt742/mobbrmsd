from . import __version__
from . import coord_generator
from . import mobbrmsd
from . import _demo
from ._mobbrmsd import *
import sys
import numpy
import networkx
import matplotlib.pyplot as plt


class _demo_mst(_demo._demo):
    def __init__(self, **kwarg):
        super().__init__(title="Minimum spanning tree with mobbrmsd", **kwarg)

    def read_input(self):
        while True:
            n_mol = _demo.readinp(
                "input number of molecules",
                6,
                check=lambda n_mol: (n_mol > 0) if isinstance(n_mol, int) else False,
            )
            n_target = _demo.readinp(
                "input number of structures",
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

    def demo(
        self, n_apm=3, n_mol=6, n_target=10, sym=((1, 2, 0), (2, 0, 1)), a=0.5, b=1.0
    ):
        import pprint

        def print_ret(i, j, ret, g):
            ev, er, ub, lb, df = (
                ret.n_eval,
                ret.eval_ratio,
                ret.bounds[0],
                ret.bounds[1],
                ret.bounds[0] - ret.bounds[1],
            )

            pre = ""
            post = ""

            for e in networkx.edges(g):
                if (e[0] == i and e[1] == j) or (e[0] == j and e[1] == i):
                    if sys.stdout.isatty():
                        pre = "\033[36m"
                        post = "**\033[0m"
                        break
                    else:
                        post = "**"

            if ub > 1.0e10:
                print(
                    pre,
                    f" {j:4d}{i:4d}{ev:12d} {er:12.6f}          +Infty{lb:16.6f}  +Infty",
                    post,
                )
            else:
                print(
                    pre,
                    f" {j:4d}{i:4d}{ev:12d} {er:12.6f}{ub:16.6f}{lb:16.6f}{df:8.3f}",
                    post,
                )
                cogen = coord_generator()
                x = numpy.array(
                    [
                        cogen.generate(n_apm, n_mol, a, b, dtype=self.prec).reshape(
                            [-1, 3]
                        )
                        for i in range(n_target)
                    ]
                )

        sep1 = "  ------------------------------------------------------------------------------"
        sep2 = "  ---------------------------------------|--------|-------------------|---------"

        print(sep1)
        print(
            "        Demonstration of minimum spanning tree construction with mobbrmsd"
        )
        print(sep1)
        print("      --System settings--")
        print(
            f"    Atoms per molecule  :{n_apm:6d}",
        )
        print(f"    Number of molecule  :{n_mol:6d}")
        print(f"    Number of structure :{n_target:6d}")

        pp = pprint.pformat(
            tuple([i for i in range(n_apm)]), width=50, compact=True
        ).split("\n")
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
        g = mrmsd.min_span_tree(x, verbose=True)
        del mrmsd
        return {"g": g}

    def after(self, g):

        if _demo.yes_or_no("Show graph ? (Open matplotlib window)"):
            n_target = len(g.nodes())

            vmax = numpy.array(
                [v for k, v in networkx.get_edge_attributes(g, "weight").items()]
            ).max()
            weights = [
                5 * (vmax - 0.95 * v) / vmax
                for k, v in networkx.get_edge_attributes(g, "weight").items()
            ]
            reverse_weights = {
                k: {"reverse_weights": 10 - 9 * v / vmax}
                for k, v in networkx.get_edge_attributes(g, "weight").items()
            }
            edge_labels = {
                k: "{:.1f}".format(v)
                for k, v in networkx.get_edge_attributes(g, "weight").items()
            }

            networkx.set_edge_attributes(g, reverse_weights)

            pos = networkx.spring_layout(g, weight="reverse_weights")
            networkx.draw_networkx_nodes(
                g,
                pos,
                node_size=int(5000 / n_target),
                node_color="white",
                edgecolors="red",
            )
            networkx.draw_networkx_labels(g, pos, font_size=int(50 / n_target) + 5)
            networkx.draw_networkx_edges(g, pos, width=weights, edge_color="tab:red")

            plt.show()
            if _demo.yes_or_no("Save graph ?"):
                path = _demo.readinp(
                    "Enter a file name",
                    "",
                    check=lambda path: (
                        (path != "") if isinstance(path, str) else False
                    ),
                )
                plt.savefig(path)

            plt.clf()
            plt.close()
            print()
