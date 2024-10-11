from . import __version__
from . import _demo
from . import coord_generator
import sys
import numpy
import pprint

title = "coord_generator demo"


class __demo__(_demo._demo):
    def __init__(self, **kwarg):
        super().__init__(title="coord_generator demo", **kwarg)

    def read_input(self):
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

    def demo(self, n_apm=3, n_mol=8, alpha=0.5, beta=1.0, n_sample=1000, **kwargs):

        n_mol_ = int(n_mol)
        n_apm_ = int(n_apm)
        alpha_ = float(alpha)
        beta_ = float(beta)
        n_sample_ = int(n_sample)

        cogen = coord_generator()
        x = cogen.generate(n_apm_, n_mol_, alpha_, beta_, n_sample=n_sample_)
        m = numpy.mean(x.reshape((-1, 3)), 0)
        s = numpy.std(x.reshape((-1, 3)), 0)
        cx = numpy.cov(x.reshape((-1, n_apm_, 3))[:, :, 0].T)
        cy = numpy.cov(x.reshape((-1, n_apm_, 3))[:, :, 1].T)
        cz = numpy.cov(x.reshape((-1, n_apm_, 3))[:, :, 2].T)
        sep = "  ------------------------------------------------------------------------------"

        print(sep)
        print(f"    Statistics with {n_sample_} samples")
        print(
            f"        n_apm = {n_apm_:4d}  n_mol = {n_mol_:4d}  alpha = {alpha_:6.2f}  beta = {beta_:6.2f}"
        )
        print()
        print(f"                              X               Y               Z")
        print(f"        mean vector {m[0]:16.9f}{m[1]:16.9f}{m[2]:16.9f}")
        print(f"        std vector  {s[0]:16.9f}{s[1]:16.9f}{s[2]:16.9f}")
        print(f"        covariance matrices")
        for i in range(n_apm_):
            for j in range(n_apm_):
                print(
                    f"                 {i+1:3d}-{j+1:3d}{cx[i,j]:16.9f}{cy[i,j]:16.9f}{cz[i,j]:16.9f}"
                )
            print()
        print(sep)

        return {
            "cogen": cogen,
            "n_apm": n_apm_,
            "n_mol": n_mol_,
            "alpha": alpha_,
            "beta": beta_,
        }

    def after(
        self,
        cogen=coord_generator(),
        n_apm=3,
        n_mol=8,
        alpha=0.5,
        beta=1.0,
        path=None,
        **kwargs,
    ):

        import matplotlib.pyplot as plt

        def onclick(event):
            ax.cla()
            x = cogen.generate(n_apm, n_mol, alpha, beta)
            for xi in x:
                ax.plot(xi[:, 0], xi[:, 1], xi[:, 2])
                ax.scatter(xi[:, 0], xi[:, 1], xi[:, 2])
            ax.set_xlim([-2.5, 2.5])
            ax.set_ylim([-2.5, 2.5])
            ax.set_zlim([-2.5, 2.5])
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_zlabel("Z")
            plt.draw()
            pass

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        fig.canvas.mpl_connect("key_press_event", onclick)
        if self.yes_or_no("Show samples ? (Open matplotlib window)"):
            x = cogen.generate(n_apm, n_mol, alpha, beta)
            for xi in x:
                ax.plot(xi[:, 0], xi[:, 1], xi[:, 2])
                ax.scatter(xi[:, 0], xi[:, 1], xi[:, 2])
            ax.set_xlim([-2.5, 2.5])
            ax.set_ylim([-2.5, 2.5])
            ax.set_zlim([-2.5, 2.5])
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_zlabel("Z")
            ax.set_box_aspect([1, 1, 1])
            if path is None:
                plt.show()
            else:
                for p in path.split(","):
                    plt.savefig(p)
            plt.clf()
            plt.close()


if __name__ == "__main__":
    main()
