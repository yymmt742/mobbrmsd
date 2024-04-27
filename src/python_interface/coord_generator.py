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
            q = 2 * self.rng.random(4) - 1.0
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
        self,
        n_apm: int,
        n_mol: int,
        a: float | Iterable[float],
        b: float | Iterable[float],
    ) -> numpy.ndarray:

        Xstr = self.rng.standard_normal((n_mol, n_apm, self.d))
        temp = self.rng.standard_normal((n_apm, self.d))
        Xtem = numpy.array([temp @ self.sog.generate() for i in range(n_mol)])
        Xvar = self.rng.standard_normal((n_mol, 1, self.d))

        def x_sample(a, b, Xvar, Xtem, Xstr):
            a_ = 0.5 * numpy.pi * a
            b_ = 0.5 * numpy.pi * b

            sa = numpy.sin(a_)
            ca = numpy.cos(a_)
            sb = numpy.sin(b_)
            cb = numpy.cos(b_)
            ret = ca * (cb * Xstr + sb * Xtem) + sa * Xvar
            return ret - numpy.mean(ret.reshape((n_apm * n_mol, self.d)), 0)

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
