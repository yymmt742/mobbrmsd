from . import _demo
from ._coord_generator import coord_generator
import sys
import numpy
import pprint

title = "coord_generator demo"


class __demo(_demo._demo):
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

        cog = coord_generator()
        x = cog.generate(
            n_apm_,
            n_mol_,
            alpha_,
            beta_,
            n_sample=n_sample_,
            remove_com=False,
        )
        m = numpy.mean(x.reshape([-1, 3]), 0)
        s = numpy.std(x.reshape([-1, 3]), 0)
        c = numpy.cov(x.reshape([-1, 3]).T)
        sep = "  ------------------------------------------------------------------------------"

        print(sep)
        print(f"    Statistics with {n_sample_} samples")
        print(
            f"        n_apm = {n_apm_:4d}  n_mol = {n_mol_:4d}  alpha = {alpha_:6.2f}  beta = {beta_:6.2f}"
        )
        print(f"                           X           Y           Z")
        print(f"        mean vector {m[0]:12.6f}{m[1]:12.6f}{m[2]:12.6f}")
        print(f"        std vector  {s[0]:12.6f}{s[1]:12.6f}{s[2]:12.6f}")
        for i, a in zip(range(3), ["X", "Y", "Z"]):
            print(f"                  {a} {c[i,0]:12.6f}{c[i,1]:12.6f}{c[i,2]:12.6f}")
        print(sep)

        return {
            "cog": cog,
            "n_apm": n_apm_,
            "n_mol": n_mol_,
            "alpha": alpha_,
            "beta": beta_,
        }

    def after(
        self,
        cog=coord_generator(),
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
            x = cog.generate(n_apm, n_mol, alpha, beta).reshape([n_mol, n_apm, 3])
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

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        fig.canvas.mpl_connect("key_press_event", onclick)
        if self.yes_or_no("Show samples ? (Open matplotlib window)"):
            x = cog.generate(n_apm, n_mol, alpha, beta).reshape([n_mol, n_apm, 3])
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
