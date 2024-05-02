from typing import Any, Iterable
import numpy


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
            while True:
                q = self.rng.standard_normal(4)
                r = numpy.linalg.norm(q)
                if r != 0:
                    break
            q /= r
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
            while True:
                A = self.rng.standard_normal(size=[d, d])
                if numpy.linalg.det(A) != 0:
                    break
            Q, R = numpy.linalg.qr(A)
            if numpy.linalg.det(Q) < 0.0:
                Q[-1] = -Q[-1]
            if A[0, 0] > 0.0:
                return Q
            else:
                return -Q

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
        n_sample: int = 1,
        temp: None | numpy.ndarray = None,
    ) -> numpy.ndarray:

        if temp == None:
            temp = self.rng.standard_normal((n_apm, self.d))

        def x_sample(
            sa: float,
            ca: float,
            sb: float,
            cb: float,
            temp: numpy.ndarray,
        ) -> numpy.ndarray:
            Xstr = self.rng.standard_normal((n_mol, n_apm, self.d))
            Xtem = numpy.array([temp @ self.sog.generate() for i in range(n_mol)])
            Xvar = self.rng.standard_normal((n_mol, 1, self.d))
            return ca * (cb * Xstr + sb * Xtem) + sa * Xvar

        def x_samples(
            a: float,
            b: float,
            n_sample: int,
            temp: numpy.ndarray,
        ) -> numpy.ndarray:
            a_ = 0.5 * numpy.pi * a
            b_ = 0.5 * numpy.pi * b
            sa = numpy.sin(a_)
            ca = numpy.cos(a_)
            sb = numpy.sin(b_)
            cb = numpy.cos(b_)
            if n_sample == 1:
                return x_sample(sa, ca, sb, cb, temp)
            else:
                return numpy.array(
                    [x_sample(sa, ca, sb, cb, temp) for i in range(n_sample)]
                )

        if isinstance(a, float):
            if isinstance(b, float):
                return x_samples(a, b, n_sample, temp)
            elif isinstance(b, Iterable):
                return numpy.array([x_samples(a, bi, n_sample, temp) for bi in b])
            elif isinstance(a, Iterable):
                if isinstance(b, float):
                    return numpy.array([x_samples(ai, b, n_sample, temp) for ai in a])
                elif isinstance(b, Iterable):
                    return numpy.array(
                        [[x_samples(ai, bi, n_sample, temp) for bi in b] for ai in a]
                    )
