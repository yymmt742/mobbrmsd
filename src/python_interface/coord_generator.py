from typing import Any, Iterable
import numpy
import matplotlib.pyplot as plt


class so_generator:
    def __init__(self, d: int = 3):

        self.rng = numpy.random.default_rng()

        def so_1():
            return numpy.array([[1.0]])

        def so_2():
            g = 2 * numpy.pi * self.rng.random(1)
            a = numpy.cos(g[0])
            b = numpy.sin(g[0])
            return numpy.array([[a, -b], [b, a]])

        def so_3():
            q = self.rng.random(4)
            q /= numpy.linalg.norm(q)
            qq = numpy.power(q, 2)
            q02 = 2 * q[0] * q[2]
            q03 = 2 * q[0] * q[3]
            q12 = 2 * q[1] * q[2]
            q13 = 2 * q[1] * q[3]
            q01 = 2 * q[0] * q[1]
            q23 = 2 * q[2] * q[3]

            return numpy.array(
                [
                    [
                        qq[0] + qq[1] - qq[2] - qq[3],
                        q12 - q03,
                        q13 + q02,
                    ],
                    [
                        q12 + q03,
                        qq[0] - qq[1] + qq[2] - qq[3],
                        q23 - q01,
                    ],
                    [
                        q13 - q02,
                        q23 + q01,
                        qq[0] - qq[1] - qq[2] + qq[3],
                    ],
                ]
            )

        def so_g():
            g = self.rng.random([d, d]) @ self.rng.random([d, d])
            u, s, v = numpy.linalg.svd(g)
            det = numpy.linalg.det(u @ v)
            if det < 0.0:
                u[-1] = -u[-1]
            return u @ v

        if d < 2:
            self.generate = so_1
        elif d == 2:
            self.generate = so_2
        elif d == 3:
            self.generate = so_3
        else:
            self.generate = so_g


class coord_generator:
    def __init__(self, d: int = 3):
        self.d = d
        self.sog = so_generator(d)
        self.rng = numpy.random.default_rng()

    def generate(
        self, n: int, m: int, a: float | Iterable[float], b: float | Iterable[float]
    ) -> numpy.ndarray:

        Xstr = self.rng.standard_normal((m, n, self.d))
        temp = self.rng.standard_normal((n, self.d))
        Xtem = numpy.array([temp @ self.sog.generate() for i in range(m)])
        Xvar = self.rng.standard_normal((m, 1, self.d))

        def x_sample(a, b, Xvar, Xtem, Xstr):
            a_ = 0.5 * numpy.pi * a
            b_ = 0.5 * numpy.pi * b

            sa = numpy.sin(a_)
            ca = numpy.cos(a_)
            sb = numpy.sin(b_)
            cb = numpy.cos(b_)

            return ca * (cb * Xstr + sb * Xtem) + sa * Xvar

        if isinstance(a, float):
            if isinstance(b, float):
                return x_sample(a, b, Xvar, Xtem, Xstr)
            elif isinstance(b, Iterable):
                return numpy.array([x_sample(a, bi, Xvar, Xtem, Xstr) for bi in b])
        elif isinstance(a, Iterable):
            if isinstance(b, float):
                return numpy.array([x_sample(ai, b, Xvar, Xtem, Xstr) for ai in a])
            elif isinstance(b, Iterable):
                return numpy.array(
                    [[x_sample(ai, bi, Xvar, Xtem, Xstr) for bi in b] for ai in a]
                )


if __name__ == "__main__":
    import itertools

    cogen = coord_generator(2)
    n = 3
    m = 5
    a = (0.0, 0.8, 0.9, 0.99)
    b = (0.0, 1.0)
    x = cogen.generate(n, m, a, b)

    print("--- coord_generator.py ---")
    print("This is sampled coordinates.\n")
    for i, j in itertools.product((range(len(a))), (range(len(b)))):
        print("    a = ", a[i], "b = ", b[j], "\n")
        print(x[i, j], "\n")

"""
        sa = "{:d}".format(round(a[i] * 100)).zfill(3)
        sb = "{:d}".format(round(b[j] * 100)).zfill(3)
        path = "sample_" + sa + "_" + sb
        plt.xlim(-3.0, 3.0)
        plt.ylim(-3.0, 3.0)
        for xk in x[i, j]:
            plt.plot(xk[:, 0], xk[:, 1])
            plt.scatter(xk[:, 0], xk[:, 1])
        plt.savefig(path + ".eps")
        plt.savefig(path + ".png")
        plt.clf()
"""
