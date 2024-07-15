import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


def teylar_0(x):
    return (
        1 / 6 * (-3 + 3 * np.sqrt(3) - np.sqrt(6))
        + 1 / 12 * (2 - np.sqrt(6)) * x
        + ((3 * np.sqrt(2) - 4) * np.power(x, 2)) / (48 * np.sqrt(3))
        + (2 / 81 - 1 / (16 * np.sqrt(6))) * np.power(x, 3)
        + (5 * (81 * np.sqrt(2) - 112) * np.power(x, 4)) / (20736 * np.sqrt(3))
    )


def spline_quat(x, f0, f1, f2):
    y = f0(x)
    y_ = f1(x)
    dx = x[1:] - x[:-1]
    dy = y[1:] - y[:-1] - dx * y_[:-1]
    m = np.zeros((3 * dx.size, 3 * dx.size))
    p = np.zeros((3 * dy.size))
    j = 0
    for dxi, dyi in zip(dx, dy):
        p[j] = dyi
        m[j, j] = dxi * dxi / 2
        m[j, j + 1] = dxi * dxi * dxi / 6
        m[j, j + 2] = dxi * dxi * dxi * dxi / 24
        if j + 3 < m.shape[1]:
            m[j + 1, j] = dxi
            m[j + 1, j + 1] = dxi * dxi / 2
            m[j + 1, j + 2] = dxi * dxi * dxi / 6
            m[j + 2, j] = 1.0
            m[j + 2, j + 1] = dxi
            m[j + 2, j + 2] = dxi * dxi / 2
            m[j + 2, j + 3] = -1.0
        else:
            m[j + 1, 2] = 1.0
            m[j + 2, j + 2] = 1.0

        j += 3

    w = np.linalg.solve(m, p).reshape([-1, 3])

    def s(t):
        if t < x[0]:
            return 0.0
        for i in range(w.size):
            if t < x[i + 1]:
                return (
                    y[i]
                    + y_[i] * (t - x[i])
                    + (w[i, 0] / 2) * np.power(t - x[i], 2)
                    + (w[i, 1] / 6) * np.power(t - x[i], 3)
                    + (w[i, 2] / 24) * np.power(t - x[i], 4)
                )

        return y[-1]

    return np.vectorize(s)


def spline_cubic(x, f0, f1):
    y = f0(x)
    dx = x[1:] - x[:-1]
    dy = y[1:] - y[:-1]
    m = np.zeros((3 * dx.size, 3 * dx.size))
    l = np.zeros((3 * dx.size, 3 * dx.size))
    p = np.zeros((3 * dy.size))
    j = 0
    for dxi, dyi in zip(dx, dy):
        p[j] = dyi
        m[j, j] = dxi
        m[j, j + 1] = dxi * dxi / 2
        m[j, j + 2] = dxi * dxi * dxi / 6
        if j + 3 < m.shape[1]:
            m[j + 1, j] = 1.0
            m[j + 1, j + 1] = dxi
            m[j + 1, j + 2] = dxi * dxi / 2
            m[j + 2, j + 1] = 1.0
            m[j + 2, j + 2] = dxi
            m[j + 1, j + 3] = -1.0
            m[j + 2, j + 4] = -1.0
        else:
            m[j + 1, 2] = 1.0
            m[j + 2, j + 2] = 1.0
        l[j, j] = dxi
        l[j, j + 1] = -dxi * dxi / 2
        l[j, j + 2] = dxi * dxi * dxi / 6
        if j > 0:
            l[j + 1, j] = 1.0
            l[j + 1, j + 1] = -dxi
            l[j + 1, j + 2] = dxi * dxi / 2
            l[j + 2, j + 1] = 1.0
            l[j + 2, j + 2] = -dxi
            l[j + 1, j - 3] = -1.0
            l[j + 2, j - 2] = -1.0
        else:
            l[1, -1] = 1.0
            l[2, 2] = 1.0
        j += 3

    w = np.linalg.solve(m, p).reshape([-1, 3])
    v = np.linalg.solve(l, p).reshape([-1, 3])

    def s0(t):
        if t < x[0]:
            return 0.0
        for i in range(w.shape[0]):
            if t < x[i + 1]:
                return (
                    y[i + 1]
                    + v[i, 0] * (t - x[i + 1])
                    + (v[i, 1] / 2) * np.power(t - x[i + 1], 2)
                    + (v[i, 2] / 6) * np.power(t - x[i + 1], 3)
                )

        return (
            y[-1]
            + v[-1, 0] * (t - x[-1])
            + (v[-1, 1] / 2) * np.power(t - x[-1], 2)
            + (v[-1, 2] / 6) * np.power(t - x[-1], 3)
        )

    def s1(t):
        for i in range(w.shape[0]):
            if t < x[i + 1]:
                return (
                    y[i]
                    + w[i, 0] * (t - x[i])
                    + (w[i, 1] / 2) * np.power(t - x[i], 2)
                    + (w[i, 2] / 6) * np.power(t - x[i], 3)
                )

        return (
            y[-1]
            + w[-1, 0] * (t - x[-1])
            + (w[-1, 1] / 2) * np.power(t - x[-1], 2)
            + (w[-1, 2] / 6) * np.power(t - x[-1], 3)
        )

    return np.vectorize(lambda t: (s0(t) + s1(t)) / 2)


f0 = np.vectorize(
    lambda x: (
        np.cos(np.arccos(x) / 3.0)
        - (np.sqrt(1 + x)) / np.sqrt(6)
        - 0.5
        + (x + 1) / 18
        - (5 * np.power(x + 1, 3 / 2)) / (108 * np.sqrt(6))
        + 2 / 243 * np.power(x + 1, 2)
        - (77 * np.power(x + 1, 5 / 2)) / (7776 * np.sqrt(6))
        + 14 * np.power(x + 1, 3) / 6561
        if np.abs(x) < 1.0
        else (
            (
                np.cosh(np.log(x + np.sqrt(np.power(x, 2) - 1)) / 3)
                - (np.sqrt(1 + x)) / np.sqrt(6)
                - 0.5
                + (x + 1) / 18
                - (5 * np.power(x + 1, 3 / 2)) / (108 * np.sqrt(6))
                + 2 / 243 * np.power(x + 1, 2)
                - (77 * np.power(x + 1, 5 / 2)) / (7776 * np.sqrt(6))
                + 14 * np.power(x + 1, 3) / 6561
            )
            if x > 0.0
            else 0.0
        )
    )
)
f1 = np.vectorize(
    lambda x: (
        (np.sin(np.arccos(x) / 3) / (3 * np.sqrt(1 - x)) - 1 / (2 * np.sqrt(6)))
        / np.sqrt(x + 1)
        + 14 * np.power(x + 1, 2) / 2187
        - 385 * np.power(x + 1, 5 / 2) / (15552 * np.sqrt(6))
        + 4 * (x + 1) / 243
        - 5 * np.sqrt(x + 1) / (72 + np.sqrt(6))
        + 1 / 18
        if np.abs(x) < 1.0
        else (
            (
                np.sinh(np.log(x + np.sqrt(np.power(x, 2) - 1)) / 3)
                / (3 * np.sqrt(np.power(x, 2) - 1))
                + 14 * 4 / 2187
                - 385 * 2 / (15552 * np.sqrt(3))
                + 8 / 243
                - 5 / (72 * np.sqrt(3))
                - 1 / (4 * np.sqrt(3))
                + 1 / 18
            )
            if x > 0.0
            else 0.0
        )
    )
)
x = np.arange(-2.0, 2.0, 0.1)
# x = np.array([-1.0, -0.75, -0.5, -0.25, 0.0])
# s = spline_quat(x, f0, f1, f2)
s = spline_cubic(x, f0, f1)
t = np.arange(-1.5, 1.5, 0.001)
"""
for xi, yi, wi in zip(x, y, w):
    s = (
        lambda t: yi
        + wi[0] * (t - xi)
        + (wi[1] / 2) * np.power(t - xi, 2)
        + (wi[2] / 6) * np.power(t - xi, 3)
    )
    plt.plot(t, s(t), label="s" + str(xi))
"""


# plt.plot(t, teylar_0(t) - f0(t), label="teylar")
plt.plot(t, s(t) - f0(t), label="S")
# plt.plot(t, s(t), label="S")
# plt.plot(t, f0(t), label="f", ls=":")
# plt.ylim([-1.1, 1.1])
# plt.ylim([-1.0e-6, 1.0e-6])
plt.legend()
plt.show()
exit()

# S = interpolate.interp1d(x, y, kind="cubic")
# c = interpolate.PPoly.from_spline(interpolate.splrep(x, y)).c.T
# print(c)
for ci, xi in zip(c, x):
    s = (
        lambda x: ci[3]
        + ci[2] * (x - xi)
        + ci[1] * np.power(x - xi, 2)
        + ci[0] * np.power(x - xi, 3)
    )
    plt.plot(t, s(t), label="p" + str(ci[0]))
plt.plot(t, f0(t), label="f", ls=":")
# plt.plot(x, p(x) - f0(x), label="S0")
# plt.plot(x, f0(x), label="f0", ls=":")
# plt.plot(x, f1(x), label="f1", ls=":")
# plt.plot(x, f2(x), label="f2", ls=":")
# plt.plot(x, np.sin(x), label="sin", ls=":")
plt.ylim((-0.5, 0.5))
# plt.ylim((-1.1, 1.1))
# plt.ylim((-10.1, 10.1))
plt.legend()
plt.show()
exit()

x = np.arange(-0.99, 0.99, 0.01)
# plt.plot(x, np.power(np.cos(np.arccos(x) / 3), 2))
p = lambda x, m, k: np.power(x, m) * np.power(1 - np.power(x, 2), -k / 2)
g = lambda x, f: (f(x + 1e-8) - f(x - 1e-8)) * 5e7
plt.plot(x, (1 / 3) * np.arccos(x), label="f")
plt.plot(x, (-1 / 3) * p(x, 0, 1), label="p")
# plt.plot(x, g(x, lambda x: (-1 / 3) * np.arccos(x)), ls="--", label="p(num)")
plt.plot(x, (-1 / 3) * p(x, 1, 3), label="p'")
# plt.plot(x, g(x, lambda x: (-1 / 3) * p(x, 0, 1)), ls="--", label="p'(num)")
plt.plot(x, (-1 / 3) * (p(x, 0, 3) + 3 * p(x, 2, 5)), label="p''")
# plt.plot(x, g(x, lambda x: (-1 / 3) * p(x, 1, 3)), ls="--", label="p''(num)")
plt.plot(x, (-1 / 3) * (9 * p(x, 1, 5) + 15 * p(x, 3, 7)), label="p''")
"""
plt.plot(
    x,
    g(x, lambda x: (-1 / 3) * (p(x, 0, 3) + 3 * p(x, 2, 5))),
    ls="--",
    label="p''(num)",
)
"""
plt.plot(
    x, (-1 / 3) * (9 * p(x, 0, 5) + 90 * p(x, 2, 7) + 105 * p(x, 4, 9)), label="p'''"
)
plt.plot(
    x,
    g(x, lambda x: (-1 / 3) * ((9 * p(x, 1, 5) + 15 * p(x, 3, 7)))),
    ls="--",
    label="p''(num)",
)
# plt.plot(x, x)
plt.ylim([-np.pi / 3, np.pi / 3])
# plt.ylim([-5 * np.pi, 5 * np.pi])
plt.legend()
plt.show()
