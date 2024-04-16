from . import __version__
from . import coord_generator
from . import mobbrmsd
import sys
import numpy
import networkx
import matplotlib.pyplot as plt

title = "Demonstration of minimum spanning tree construction with mobbrmsd"

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

    print(
        pre,
        f" {j:4d}{i:4d}{ev:12d} {er:12.6f}{ub:16.6f}{lb:16.6f}{df:8.3f}",
        post,
    )


def main(n_apm=3, n_mol=6, n_target=20, sym=[[1, 2, 0], [2, 0, 1]], a=0.5, b=1.0):
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
    print("        Demonstration of minimum spanning tree construction with mobbrmsd")
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

    g, states = mrmsd.min_span_tree(x)
    mrmsd.clear()
    print(sep1)

    print("")
    print("     i   j      N_eval   Eval_ratio      Upperbound      Lowerbound    Gap")
    print(sep1)
    for i in range(len(states) - 1):
        for j in range(i + 1, len(states)):
            print_ret(i, j, states[i][j], g)
        print()
    print(sep1)
    return g


def show_graph(g):

    while True:
        print("\r  show graph ? ['y', 'n'] >> ", end="")
        inp = input()
        if inp == "y":
            break
        if inp == "n":
            print()
            return

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
        g, pos, node_size=int(5000 / n_target), node_color="white", edgecolors="red"
    )
    networkx.draw_networkx_labels(g, pos, font_size=int(50 / n_target) + 5)

    networkx.draw_networkx_edges(g, pos, width=weights, edge_color="tab:red")
    # networkx.draw_networkx_edge_labels(
    #    g, pos, edge_labels, font_size=int(50 / n_target) + 5, rotate=False
    # )
    plt.show()
    print()


if __name__ == "__main__":
    main()
