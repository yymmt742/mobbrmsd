from . import __version__
from . import _demo
from . import coord_generator
import sys
import numpy
import pprint

title = "coord_generator demo"


def read_input() -> tuple:
    import itertools
    import math

    def isfloat(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    while True:
        n_mol = _demo.readinp(
            f"input number of molecules",
            6,
            check=lambda n_mol: (n_mol > 0) if isinstance(n_mol, int) else False,
        )
        n_apm = _demo.readinp(
            f"input number of atoms per molecule",
            3,
            check=lambda n_apm: ((n_apm > 0) if isinstance(n_apm, int) else False),
        )
        alpha = _demo.readinp(
            f"input alpha",
            0.8,
            check=lambda alpha: isinstance(alpha, float),
        )
        beta = _demo.readinp(
            f"input beta",
            1.0,
            check=lambda beta: isinstance(beta, float),
        )
        break

    return {"n_apm": n_apm, "n_mol": n_mol, "alpha": alpha, "beta": beta}


def main(n_apm=3, n_mol=8, alpha=0.5, beta=1.0):
    cogen = coord_generator()
    x = cogen.generate(n_apm, n_mol, alpha, beta, n_sample=10000)
    m = numpy.mean(x.reshape((-1, 3)), 0)
    s = numpy.std(x.reshape((-1, 3)), 0)
    cx = numpy.cov(x.reshape((-1, n_apm, 3))[:, :, 0].T)
    cy = numpy.cov(x.reshape((-1, n_apm, 3))[:, :, 1].T)
    cz = numpy.cov(x.reshape((-1, n_apm, 3))[:, :, 2].T)
    sep = "  ------------------------------------------------------------------------------"

    print(sep)
    print(f"    Statistics with 10000 samples")
    print(
        f"        n_apm = {n_apm:4d}  n_mol = {n_mol:4d}  alpha = {alpha:6.2f}  beta = {beta:6.2f}"
    )
    print()
    print(f"                              X               Y               Z")
    print(f"        mean vector {m[0]:16.9f}{m[1]:16.9f}{m[2]:16.9f}")
    print(f"        std vector  {s[0]:16.9f}{s[1]:16.9f}{s[2]:16.9f}")
    print(f"        covariance matrices")
    for i in range(n_apm):
        for j in range(n_apm):
            print(
                f"                 {i+1:d}-{j+1:d}{cx[i,j]:16.9f}{cy[i,j]:16.9f}{cz[i,j]:16.9f}"
            )
        print()
    print(sep)

    return {
        "cogen": cogen,
        "n_apm": n_apm,
        "n_mol": n_mol,
        "alpha": alpha,
        "beta": beta,
    }


def show(ret):
    import matplotlib.pyplot as plt

    def onclick(event):
        ax.cla()
        x = ret["cogen"].generate(ret["n_apm"], ret["n_mol"], ret["alpha"], ret["beta"])
        for xi in x:
            ax.plot(xi[:, 0], xi[:, 1], xi[:, 2])
            ax.scatter(xi[:, 0], xi[:, 1], xi[:, 2])
        ax.set_xlim([-2.5, 2.5])
        ax.set_ylim([-2.5, 2.5])
        ax.set_zlim([-2.5, 2.5])
        ax.set_box_aspect([1, 1, 1])
        plt.draw()
        pass

    inp = ""

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    fig.canvas.mpl_connect("key_press_event", onclick)
    while True:
        inp = input("  Show samples ? (Open matplotlib window) ['y'es, 'n'o] >> ")
        if inp == "":
            continue
        break
    if inp[0] == "q" or inp[0] == "Q":
        exit()
    elif inp[0] == "y" or inp[0] == "Y":
        x = ret["cogen"].generate(ret["n_apm"], ret["n_mol"], ret["alpha"], ret["beta"])
        for xi in x:
            ax.plot(xi[:, 0], xi[:, 1], xi[:, 2])
            ax.scatter(xi[:, 0], xi[:, 1], xi[:, 2])
        ax.set_xlim([-2.5, 2.5])
        ax.set_ylim([-2.5, 2.5])
        ax.set_zlim([-2.5, 2.5])
        ax.set_box_aspect([1, 1, 1])
        plt.show()
        plt.clf()
        plt.close()
    print()


if __name__ == "__main__":
    main()
