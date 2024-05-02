from . import __version__
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
        while True:
            inp = input("    input number of molecules (default : 6)          >> ")
            if inp == "":
                n_mol = 6
                break
            if inp[0] == "q" or inp[0] == "Q":
                exit()
            elif inp.isdigit():
                n_mol = int(inp)
                if n_mol > 0:
                    break
        while True:
            inp = input("    input number of atoms per molecule (default : 3) >> ")
            if inp == "":
                n_apm = 3
                break
            if inp[0] == "q" or inp[0] == "Q":
                exit()
            elif inp.isdigit():
                n_apm = int(inp)
                if n_apm > 1:
                    break
        while True:
            inp = input("    input alpha (default : 0.8) >> ")
            if inp == "":
                alpha = 0.8
                break
            if inp[0] == "q" or inp[0] == "Q":
                exit()
            elif isfloat(inp):
                alpha = float(inp)
                break
        while True:
            inp = input("    input beta (default : 1.0) >> ")
            if inp == "":
                beta = 1.0
                break
            if inp[0] == "q" or inp[0] == "Q":
                exit()
            elif isfloat(inp):
                beta = float(inp)
                break
        break

    return {"n_apm": n_apm, "n_mol": n_mol, "alpha": alpha, "beta": beta}


def main(n_apm=3, n_mol=8, alpha=0.5, beta=1.0):
    cogen = coord_generator()
    x = cogen.generate(n_apm, n_mol, alpha, beta, n_sample=1000)
    m = numpy.mean(x.reshape((-1, 3)), 0)
    s = numpy.std(x.reshape((-1, 3)), 0)
    cx = numpy.cov(x.reshape((-1, n_apm, 3))[:, :, 0].T)
    cy = numpy.cov(x.reshape((-1, n_apm, 3))[:, :, 1].T)
    cz = numpy.cov(x.reshape((-1, n_apm, 3))[:, :, 2].T)

    print(f"    Statistics with 1000 samples")
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
        inp = input("  Show coordinate ? (Open matplotlib window) ['y'es, 'n'o] >> ")
        if inp == "":
            continue
        if inp[0] == "n" or inp[0] == "N":
            break
        elif inp[0] == "q" or inp[0] == "Q":
            exit()
        elif inp[0] == "y" or inp[0] == "Y":
            x = ret["cogen"].generate(
                ret["n_apm"], ret["n_mol"], ret["alpha"], ret["beta"]
            )
            for xi in x:
                ax.plot(xi[:, 0], xi[:, 1], xi[:, 2])
                ax.scatter(xi[:, 0], xi[:, 1], xi[:, 2])
            ax.set_xlim([-2.5, 2.5])
            ax.set_ylim([-2.5, 2.5])
            ax.set_zlim([-2.5, 2.5])
            ax.set_box_aspect([1, 1, 1])
            plt.show()
    print()


if __name__ == "__main__":
    main()
