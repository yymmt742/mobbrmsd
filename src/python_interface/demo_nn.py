from . import __version__
from . import coord_generator
from . import mobbrmsd
import sys
import numpy
import networkx
import matplotlib.pyplot as plt


def print_ret(i, j, ret, edges):
    ev, er, ub, lb, df = (
        ret.n_eval,
        ret.eval_ratio,
        ret.bounds[0],
        ret.bounds[1],
        ret.bounds[0] - ret.bounds[1],
    )

    pre  = ''
    post = ''
    for e in edges:
      if (e[0]==i and e[1]==j) or (e[0]==j and e[1]==i):
        if sys.stdout.isatty():
          pre  = '\033[36m'
          post = '**\033[0m'
          break
        else:
          post = '**'

    print(
        pre,
        f" {j:4d}{i:4d}{ev:12d} {er:12.6f}{ub:12.3f}{lb:12.3f}{df:8.3f}",
        post,
    )


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

    g, states =mrmsd.min_span_tree(x)
    print(sep1)

    for e, w in  zip(edges, weights):
      print(f"{e[0]:8d} {e[1]:8d} {w:12.3f}")

    print("")
    print("     i   j      N_eval   Eval_ratio  Upperbound  Lowerbound    Gap")
    print(sep1)
    for i in range(len(states)-1):
      for j in range(i+1, len(states)):
        print_ret(i, j, states[i][j], edges)
      print()
    print(sep1)

    g = networkx.Graph()
    for e, w in  zip(edges, weights):
      g.add_edge(e[0], e[1], weight=w)

    pos = networkx.spring_layout(g)
    edge_labels = {k: '{:.1f}'.format(v) for k, v in networkx.get_edge_attributes(g, "weight").items()}
    networkx.draw_networkx_nodes(g, pos, node_size=500,node_color="white", edgecolors="red")
    networkx.draw_networkx_labels(g, pos)
    networkx.draw_networkx_edges(g, pos, width=5*((weights.max()-weights)/weights.max()), edge_color="tab:red")
    networkx.draw_networkx_edge_labels(g, pos, edge_labels)
    plt.show()


if __name__ == "__main__":
    main()
