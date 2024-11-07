from . import _demo
from ._coord_generator import coord_generator
from ..dataclass import molecules, molecular_system
from .._mobbrmsd import mobbrmsd
import sys
import numpy
import networkx
import matplotlib.pyplot as plt


class __demo(_demo._demo):
    def __init__(self, **kwarg):
        super().__init__(title="Minimum spanning tree", **kwarg)

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
                if not self.yes_or_no(
                    "This parameter may take time to compute. May this be run ?"
                ):
                    continue
            break

        print()
        return {"n_mol": n_mol, "n_target": n_target}

    def demo(
        self,
        n_apm=3,
        n_mol=6,
        n_sym=2,
        n_target=10,
        a=0.5,
        b=1.0,
        **kwargs,
    ):
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

        n_mol_ = int(n_mol)
        n_apm_ = int(n_apm)
        sym = _demo.generate_sym_indices(n_apm_, int(n_sym))
        n_target_ = int(n_target)
        a_ = float(a)
        b_ = float(b)

        cogen = coord_generator()
        x = numpy.array(
            [
                cogen.generate(n_apm_, n_mol_, a_, b_, dtype=self.prec).reshape([-1, 3])
                for i in range(n_target_)
            ]
        )

        mols = molecules(n_apm=n_apm_, n_mol=n_mol_, sym=sym)
        _demo.print_system(
            molecular_system([mols]),
            "Demonstration of minimum spanning tree construction with mobbrmsd",
        )
        mrmsd = mobbrmsd(mols=mols)
        g = mrmsd.min_span_tree(x, verbose=True)
        del mrmsd
        return {"g": g}

    def after(self, g=None, path=None, **kwargs):

        if g is None:
            return
        if self.yes_or_no("Show graph? (Open matplotlib window)"):
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

            if path is None:
                plt.show()
            else:
                for p in path.split(","):
                    plt.savefig(p)

            plt.clf()
            plt.close()
