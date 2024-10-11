from . import __version__
from . import _demo
from .coord_generator import coord_generator
from ._mobbrmsd import *
import sys
import numpy


class __demo__(_demo._demo):
    def __init__(self, **kwarg):
        super().__init__(title="mobbRMSD basic", **kwarg)

    def read_input(self):
        import itertools
        import math

        while True:
            n_mol = _demo.readinp(
                "input number of molecules",
                6,
                check=lambda n_mol: (n_mol > 0) if isinstance(n_mol, int) else False,
            )
            if n_mol > 8:
                if not self.yes_or_no(
                    "This parameter may take time to compute. May this be run?"
                ):
                    continue

            n_apm = _demo.readinp(
                "input number of atoms per molecule",
                3,
                check=lambda n_apm: ((n_apm > 0) if isinstance(n_apm, int) else False),
            )

            def msg():
                print(
                    "    number of molecular symmetry must be less than factrial to atoms per molecule."
                )
                return False

            n_sym = _demo.readinp(
                "input number of molecular symmetry",
                n_apm,
                check=lambda n_sym: (
                    (n_sym > 0) & (True if n_sym <= math.factorial(n_apm) else msg())
                    if isinstance(n_sym, int)
                    else False
                ),
            )

            cost = math.factorial(n_mol) * n_sym**n_mol
            if cost > 10000000:
                if not self.yes_or_no(
                    f"This parameter may take time to compute (cost is {cost}). May this be run ?"
                ):
                    continue

            break

        print()

        return {"n_apm": n_apm, "n_mol": n_mol, "n_sym": n_sym}

    def demo(self, n_apm=3, n_mol=8, n_sym=2, a=0.8, b=1.0, **kwarg):
        import pprint

        n_mol_ = int(n_mol)
        n_apm_ = int(n_apm)
        sym = _demo.generate_sym_indices(n_apm_, int(n_sym))
        a_ = float(a)
        b_ = float(b)

        cogen = coord_generator()
        x = cogen.generate(n_apm_, n_mol_, a_, b_, dtype=self.prec).reshape([-1, 3])
        y = cogen.generate(n_apm_, n_mol_, a_, b_, dtype=self.prec).reshape([-1, 3])
        x -= numpy.mean(x, 0)
        y -= numpy.mean(y, 0)
        z = y.copy()

        molecules = DataclassMolecule(n_apm=n_apm_, n_mol=n_mol_, sym=sym)
        _demo.print_system(
            [molecules], title="Demonstration of mobbRMSD with random coordinates"
        )
        mrmsd = mobbrmsd(molecules=molecules)

        ub, lb = numpy.inf, -numpy.inf
        ret = mrmsd.run(x, y, maxeval=0, get_rotation=True)

        xtra = ["__-¯¯", "-__-¯", "¯-__-", "¯¯-__", "-¯¯-_", "_-¯¯-"]
        erace = "     "

        i = 0
        _demo.print_ret(ret, end="", header=True)
        while not ret.is_finished():
            _demo.print_ret(ret, post=xtra[int(i / 5000) % 6], end="", to_console=True)
            if ub > ret.upperbound() or lb < ret.lowerbound():
                _demo.print_ret(ret, post=erace)
            ub, lb = ret.upperbound(), ret.lowerbound()
            ret.restart(maxeval=0)
            i += 1
        ret.restart(maxeval=0, get_rotation=True)
        _demo.print_ret(ret, post=erace, header=True, footer=True)

        y = ret.rotate_y(y)

        sep1 = "  ------------------------------------------------------------------------------"
        sep2 = "  ---------------------------------------|--------|-------------------|---------"

        print(sep1)
        print(
            "       Reference     | Target (original) |  disp. |  Target (rotate)  |  disp."
        )
        d1, d2 = 0.0, 0.0
        for xi, yi, zi in zip(
            x.reshape([n_mol_, n_apm_, 3]),
            y.reshape([n_mol_, n_apm_, 3]),
            z.reshape([n_mol_, n_apm_, 3]),
        ):
            print(sep2)
            for xij, yij, zij in zip(xi, yi, zi):
                d1_ = numpy.sum(numpy.power(xij - zij, 2))
                d2_ = numpy.sum(numpy.power(xij - yij, 2))
                d1 += d1_
                d2 += d2_
                print(
                    f"  {xij[0]:6.2f}{xij[1]:6.2f}{xij[2]:6.2f}",
                    f"|{yij[0]:6.2f}{yij[1]:6.2f}{yij[2]:6.2f} |{d1_:7.2f}",
                    f"|{zij[0]:6.2f}{zij[1]:6.2f}{zij[2]:6.2f} |{d2_:7.2f}",
                )
        print(sep2)
        d1 *= 1.0 / (n_mol_ * n_apm_)
        d2 *= 1.0 / (n_mol_ * n_apm_)
        print(
            f"          mean squared deviation         |{d1:7.2f} |                   |{d2:7.2f}"
        )
        d1, d2 = numpy.sqrt(d1), numpy.sqrt(d2)
        print(
            f"        root mean squared deviation      |{d1:7.2f} |                   |{d2:7.2f}"
        )
        print(sep1)

        return {
            "x": x.reshape([n_mol_, n_apm_, 3]),
            "y": y.reshape([n_mol_, n_apm_, 3]),
            "z": z.reshape([n_mol_, n_apm_, 3]),
        }

    def after(
        self,
        x=[[[0.0, 0.0, 0.0]]],
        y=[[[0.0, 0.0, 0.0]]],
        z=[[[0.0, 0.0, 0.0]]],
        path=None,
        **kwargs,
    ):
        ang = 0

        if self.yes_or_no("Show samples ? (Open matplotlib window)"):
            import matplotlib.pyplot as plt
            import matplotlib.animation as animation

            cmap = plt.get_cmap("tab20")
            fig = plt.figure(figsize=(8, 8 / 1.618))
            axes = [
                fig.add_subplot(1, 2, 1, projection="3d"),
                fig.add_subplot(1, 2, 2, projection="3d"),
            ]
            for ax, ref, tgt in zip(axes, [x, x], [z, y]):
                for i, xy in enumerate(zip(ref, tgt)):
                    for xi, yi in zip(xy[0], xy[1]):
                        ax.plot(
                            [xi[0], yi[0]],
                            [xi[1], yi[1]],
                            [xi[2], yi[2]],
                            color=cmap(2 * i),
                            ls=":",
                            lw=1.0,
                        )
                for crd, fillstyle, ms, j in zip(
                    [ref, tgt], ["none", "full"], [6, 6], [1, 0]
                ):
                    for i, xi in enumerate(crd):
                        ax.plot(
                            xi[:, 0],
                            xi[:, 1],
                            xi[:, 2],
                            color=cmap(2 * i + j),
                            marker="o",
                            fillstyle=fillstyle,
                            ms=ms,
                        )
                ax.set_xlabel("X")
                ax.set_ylabel("Y")
                ax.set_zlabel("Z")
                ax.set_xlim([-2.5, 2.5])
                ax.set_ylim([-2.5, 2.5])
                ax.set_zlim([-2.5, 2.5])
                ax.view_init(azim=15)
                ax.set_box_aspect([1, 1, 1])
            plt.tight_layout()

            def rot(frame):
                axes[0].view_init(azim=frame)
                axes[1].view_init(azim=frame)

            if path is None:
                ani = animation.FuncAnimation(
                    fig, rot, frames=numpy.arange(15, 374, 2), interval=1
                )
                plt.show()
            else:
                plt.savefig(path)
            plt.clf()
            plt.close()
