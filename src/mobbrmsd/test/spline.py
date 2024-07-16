import numpy as np
import matplotlib.pyplot as plt

# from scipy import interpolate


# expand np.cos(np.arccos(x) / 3.0)
def expand(x):
    return (
        0.5
        + (np.sqrt(1 + x)) / np.sqrt(6)
        - (x + 1) / 18
        + (5 * np.power(x + 1, 3 / 2)) / (108 * np.sqrt(6))
        - 2 / 243 * np.power(x + 1, 2)
        + (77 * np.power(x + 1, 5 / 2)) / (7776 * np.sqrt(6))
        - 14 * np.power(x + 1, 3) / 6561
        + (2431 * np.power(x + 1, 7 / 2)) / (839808 * np.sqrt(6))
        - 40 * np.power(x + 1, 4) / 59049
    )


# newton np.cos(np.arccos(x) / 3.0)
def newton_(x):
    y = np.cos(np.arccos(x) / 3)
    s = 1
    while s > 1.0e-16:
        f = 4 * y**3 - 3 * y - x
        df = 12 * y**2 - 3
        s = f / df
        y -= s
    return y


def spline_cubic(x, f0):
    y = f0(x)
    dx = x[1:] - x[:-1]
    dm = (x[1:] + x[:-1]) / 2
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
            m[j + 1, 5] = -1.0
            m[j + 2, j - 1] = 1.0
            m[j + 2, j + 2] = -1.0

            # m[j + 1, 2] = 1.0
            # m[j + 2, j + 2] = 1.0

            # m[j + 1, 1] = 1.0
            # m[j + 2, j + 1] = 1.0
            # p[j + 1] = f1(x[0])
            # p[j + 2] = f1(x[-1])
        j += 3

    w = np.linalg.solve(m, p).reshape([-1, 3])

    def s0(t):
        if t < x[0]:
            return 0.0
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

    def corr2(t):
        for i in range(w.shape[0]):
            if t < x[i + 1]:
                diff = -4 * (f0(dm[i]) - s0(dm[i])) / dx[i] ** 2
                return diff * (t - x[i]) * (t - x[i + 1])
        return 0.0

    def corr3(t):
        for i in range(w.shape[0]):
            if t < x[i + 1]:
                dl = dm[i] - np.sqrt(3) * dx[i] / 6
                du = dm[i] + np.sqrt(3) * dx[i] / 6
                diff = (
                    (f0(dl) - s0(dl) - corr2(dl) - f0(du) + s0(du) + corr2(du))
                    * 18
                    / (np.sqrt(3) * dx[i] ** 3)
                )
                return diff * (t - x[i]) * (t - dm[i]) * (t - x[i + 1])
        return 0.0

    def corr4(t):
        for i in range(w.shape[0]):
            if t < x[i + 1]:
                dl = dm[i] - np.sqrt(2) * dx[i] / 4
                du = dm[i] + np.sqrt(2) * dx[i] / 4
                diff = (
                    -32
                    * (
                        f0(du)
                        - s0(du)
                        - corr2(du)
                        - corr3(du)
                        + f0(dl)
                        - s0(dl)
                        - corr2(dl)
                        - corr3(dl)
                    )
                    / dx[i] ** 4
                )
                return diff * (t - x[i]) * (t - x[i + 1]) * (t - dm[i]) ** 2
        return 0.0

    return (
        np.vectorize(s0),
        np.vectorize(corr2),
        np.vectorize(corr3),
        np.vectorize(corr4),
    )


f0 = np.vectorize(lambda x: (newton_(x) - expand(x) if x > -1.0 else 0.0))

x = 1.0 * np.cos(np.pi * np.linspace(1, 0, 128))
s0, corr2, corr3, corr4 = spline_cubic(x, f0)
f = lambda t: s0(t) + corr4(t) + corr3(t) + corr2(t) + expand(t)
t = np.arange(-1.0, 1.0, 0.001)


"""
S = interpolate.make_interp_spline(
    x,
    f0(x),
    k=5,
    # bc_type="clamped",
)
plt.plot(t, S(t) - f0(t), label="interp1d")
"""
# plt.plot(t, s0(t) - f0(t), label="S")
# plt.plot(t, corr2(t), label="corr2")
# plt.plot(t, s0(t) + corr2(t) - f0(t), label="S+corr2")
# plt.plot(t, corr3(t), label="corr3")
# plt.plot(t, s0(t) + corr2(t) + corr3(t) - f0(t), label="S+corr2+corr3")
# plt.plot(t, corr4(t), label="corr4")
plt.plot(t, s0(t) + corr4(t) + corr3(t) + corr2(t) - f0(t), label="S+corr2+corr3+corr4")
# plt.stem(x, np.ones(x.size) * np.max(S(t) - f0(t)))
# xm = (x[1:] + x[:-1]) / 2
# plt.stem(xm, np.ones(xm.size) * np.max(S(t) - f0(t)), linefmt=":")
# plt.plot(t, s0(t), label="S")
# plt.plot(t, f0(t), label="f", ls=":")
# plt.ylim([-1.1, 1.1])
# plt.ylim([-1.0e-14, 1.0e-14])
plt.plot(t, 4 * np.power(f(t), 3) - 3 * f(t) - t, label="4y(x)^3-3y(x)^3=x")
plt.legend()
plt.show()
