import numpy as np
import matplotlib.pyplot as plt

# from scipy import interpolate


# expand np.cos(np.arccos(x) / 3.0)
def Taylor0(x):
    return (
        np.sqrt(3) / 2
        + x / 6
        - np.power(x, 2) / (12 * np.sqrt(3))
        + 2 * np.power(x, 3) / 81
        - 35 * np.power(x, 4) / (1296 * np.sqrt(3))
    )


def Taylor1(x):
    return (
        1
        + (x - 1) / 9
        - 4 * np.power(x - 1, 2) / 243
        + 28 * np.power(x - 1, 3) / 6561
        - 80 * np.power(x - 1, 4) / 59049
    )


def Puiseux_1(x):
    return (
        1 / 2
        + (np.sqrt(1 + x)) / np.sqrt(6)
        - (x + 1) / 18
        + (5 * np.power(x + 1, 3 / 2)) / (108 * np.sqrt(6))
        - 2 / 243 * np.power(x + 1, 2)
        + (77 * np.power(x + 1, 5 / 2)) / (7776 * np.sqrt(6))
        - 14 * np.power(x + 1, 3) / 6561
        + (2431 * np.power(x + 1, 7 / 2)) / (839808 * np.sqrt(6))
        - 40 * np.power(x + 1, 4) / 59049
        + (1062347 * np.power(x + 1, 9 / 2)) / (1088391168 * np.sqrt(6))
    )


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
def newton_(x, y0, k):
    if x <= -1.0:
        return 1 / 2
    i = 0
    y = y0
    s = 1
    while i < k and s > 1e-17:
        yy = y * y
        f = (4 * yy - 3) * y - x
        df = np.max([12 * y**2 - 3, 1e-18])
        s = f / df
        y -= s
        i += 1
    return y


def spline_linear(x, f):
    dx = x[1:] - x[:-1]
    dm = (x[1:] + x[:-1]) / 2
    y = f(x)
    dy = y[1:] - y[:-1]
    g = dy / dx
    b = (x[1:] * y[:-1] - y[1:] * x[:-1]) / dx

    def s0(t):
        if t < x[0]:
            return 0.0
        for i in range(dx.shape[0]):
            if t < x[i + 1]:
                return b[i] + g[i] * t
        return b[-1] + g[-1] * t

    return np.vectorize(s0)


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
            # not a not
            m[j + 1, 2] = 1.0
            m[j + 1, 5] = -1.0
            m[j + 2, j - 1] = 1.0
            m[j + 2, j + 2] = -1.0

            # natural
            # m[j + 1, 2] = 1.0
            # m[j + 2, j + 2] = 1.0

            # cramped
            # m[j + 1, 1] = 1.0
            # m[j + 2, j + 1] = 1.0
            # p[j + 1] = f1(x[0])
            # p[j + 2] = f1(x[-1])

        j += 3

    w = np.linalg.solve(m, p).reshape([-1, 3])

    z = np.zeros((w.shape[0], 5))
    for xi, yi, zi, wi in zip(x, y, z, w):
        zi[0] = yi - wi[0] * xi + wi[1] * xi * xi / 2 - wi[2] * xi * xi * xi / 6
        zi[1] = wi[0] - wi[1] * xi + wi[2] * xi * xi / 2
        zi[2] = wi[1] / 2 - wi[2] * xi / 2
        zi[3] = wi[2] / 6

    def s0(t):
        if t < x[0]:
            return 0.0
        for i in range(w.shape[0]):
            if t < x[i + 1]:
                return (
                    z[i, 0]
                    + z[i, 1] * t
                    + z[i, 2] * np.power(t, 2)
                    + z[i, 3] * np.power(t, 3)
                )
                """
                return (
                    y[i]
                    + w[i, 0] * (t - x[i])
                    + (w[i, 1] / 2) * np.power(t - x[i], 2)
                    + (w[i, 2] / 6) * np.power(t - x[i], 3)
                )
                """

        return (
            y[-1]
            + w[-1, 0] * (t - x[-1])
            + (w[-1, 1] / 2) * np.power(t - x[-1], 2)
            + (w[-1, 2] / 6) * np.power(t - x[-1], 3)
        )

    return np.vectorize(s0), z


def corr(x, f, s):
    dx = x[1:] - x[:-1]
    dm = (x[1:] + x[:-1]) / 2
    y = f(x)
    dy = y[1:] - y[:-1]
    v2 = np.zeros((dx.shape[0], 5))
    v3 = np.zeros((dx.shape[0], 5))
    v4 = np.zeros((dx.shape[0], 5))

    d2 = 16 * (f(dm) - s(dm)) / np.power(dx, 4)

    for i in range(dx.shape[0]):
        v2[i, 0] = d2[i] * x[i] * x[i] * x[i + 1] * x[i + 1]
        v2[i, 1] = -2 * d2[i] * x[i] * x[i + 1] * (x[i] + x[i + 1])
        v2[i, 2] = d2[i] * (4 * x[i] * x[i + 1] + x[i] * x[i] + x[i + 1] * x[i + 1])
        v2[i, 3] = -2 * d2[i] * (x[i] + x[i + 1])
        v2[i, 4] = d2[i]

    def corr2(t):
        for i in range(dx.shape[0]):
            if t < x[i + 1]:
                return (
                    v2[i, 0]
                    + v2[i, 1] * t
                    + v2[i, 2] * np.power(t, 2)
                    + v2[i, 3] * np.power(t, 3)
                    + v2[i, 4] * np.power(t, 4)
                )
                # return d2[i] * np.power(t - x[i], 2) * np.power(t - x[i + 1], 2)
        return 0.0

    l3 = dm - np.sqrt(3) * dx / 6
    u3 = dm + np.sqrt(3) * dx / 6
    d3 = (
        (
            f(l3)
            - s(l3)
            - np.vectorize(corr2)(l3)
            - f(u3)
            + s(u3)
            + np.vectorize(corr2)(u3)
        )
        * 18
        / (np.sqrt(3) * np.power(dx, 3))
    )

    for i in range(dx.shape[0]):
        v3[i, 0] = -d3[i] * x[i] * x[i + 1] * (x[i] + x[i + 1]) / 2
        v3[i, 1] = d3[i] * (x[i] * x[i] + 4 * x[i] * x[i + 1] + x[i + 1] * x[i + 1]) / 2
        v3[i, 2] = -3 * d3[i] * (x[i] + x[i + 1]) / 2
        v3[i, 3] = d3[i]

    def corr3(t):
        for i in range(dx.shape[0]):
            if t < x[i + 1]:
                # return d3[i] * (t - x[i]) * (t - dm[i]) * (t - x[i + 1])
                return (
                    v3[i, 0]
                    + v3[i, 1] * t
                    + v3[i, 2] * np.power(t, 2)
                    + v3[i, 3] * np.power(t, 3)
                    + v3[i, 4] * np.power(t, 4)
                )
        return 0.0

    l4 = dm - np.sqrt(2) * dx / 4
    u4 = dm + np.sqrt(2) * dx / 4
    d4 = (
        -32
        * (
            f(u4)
            - s(u4)
            - np.vectorize(corr2)(u4)
            - np.vectorize(corr3)(u4)
            + f(l4)
            - s(l4)
            - np.vectorize(corr2)(l4)
            - np.vectorize(corr3)(l4)
        )
        / np.power(dx, 4)
    )

    for i in range(dx.shape[0]):
        v4[i, 0] = d4[i] * (x[i] * x[i + 1] * dm[i] * dm[i])
        v4[i, 1] = -d4[i] * dm[i] * (dm[i] * (x[i] + x[i + 1]) + 2 * x[i] * x[i + 1])
        v4[i, 2] = d4[i] * (x[i] * x[i + 1] + (2 * x[i] + 2 * x[i] + dm[i]) * dm[i])
        v4[i, 3] = -d4[i] * (x[i] + x[i + 1] + 2 * dm[i])
        v4[i, 4] = d4[i]

    def corr4(t):
        for i in range(dx.shape[0]):
            if t < x[i + 1]:
                # return d4[i] * (t - x[i]) * (t - x[i + 1]) * (t - dm[i]) ** 2
                return (
                    v4[i, 0]
                    + v4[i, 1] * t
                    + v4[i, 2] * np.power(t, 2)
                    + v4[i, 3] * np.power(t, 3)
                    + v4[i, 4] * np.power(t, 4)
                )
        return 0.0

    return (
        np.vectorize(corr2),
        np.vectorize(corr3),
        np.vectorize(corr4),
        v2 + v3,  # + v4
    )


n = 16
x0 = np.cos(np.pi * np.linspace(1, 0, n)) / 3 - 2 / 3
x1 = np.cos(np.pi * np.linspace(1, 0, n)) / 3
x2 = np.cos(np.pi * np.linspace(1, 0, n)) / 3 + 2 / 3
f0 = np.vectorize(
    lambda x: (np.cos(np.arccos(x) / 3) - Puiseux_1(x) if x > -1.0 else 0.0)
)
f1 = np.vectorize(lambda x: np.cos(np.arccos(x) / 3) - Taylor0(x))
f2 = np.vectorize(lambda x: np.cos(np.arccos(x) / 3) - Taylor1(x))

# c0 = spline_linear(x0, f0)
# c1 = spline_linear(x1, f1)
# c2 = spline_linear(x2, f2)
c0, z0 = spline_cubic(x0, f0)
c1, z1 = spline_cubic(x1, f1)
c2, z2 = spline_cubic(x2, f2)
r02, r03, r04, v0 = corr(x0, f0, c0)
r12, r13, r14, v1 = corr(x1, f1, c1)
r22, r23, r24, v2 = corr(x2, f2, c2)

w0 = np.array(
    [
        1 / 2 - 1 / 18 - 2 / 243 - 14 / 6561 - 40 / 59049,
        -1 / 18 - 2 * 2 / 243 - 3 * 14 / 6561 - 4 * 40 / 59049,
        -2 / 243 - 3 * 14 / 6561 - 6 * 40 / 59049,
        -14 / 6561 - 4 * 40 / 59049,
        -40 / 59049,
    ]
)
ww0 = np.array(
    [
        1 + 5 / 108 + 77 / 7776 + 2431 / 839808 + 1062347 / 1088391168,
        5 / 108 + 2 * 77 / 7776 + 3 * 2431 / 839808 + 4 * 1062347 / 1088391168,
        77 / 7776 + 3 * 2431 / 839808 + 6 * 1062347 / 1088391168,
        2431 / 839808 + 4 * 1062347 / 1088391168,
        1062347 / 1088391168,
    ]
) / np.sqrt(6)
w1 = np.array(
    [
        np.sqrt(3) / 2,
        1 / 6,
        -1 / (12 * np.sqrt(3)),
        2 / 81,
        -35 / (1296 * np.sqrt(3)),
    ]
)
w2 = np.array(
    [
        1 - 1 / 9 - 4 / 243 - 28 / 6561 - 80 / 59049,
        1 / 9 + 2 * 4 / 243 + 3 * 28 / 6561 + 4 * 80 / 59049,
        -4 / 243 - 3 * 28 / 6561 - 6 * 80 / 59049,
        28 / 6561 + 4 * 80 / 59049,
        -80 / 59049,
    ]
)
z0 += v0
z0 += w0
z1 += v1
z1 += w1
z2 += v2
z2 += w2
print(x0)
print(z0)
print(x1)
print(z1)
print(x2)
print(z2)


def g(t):
    return np.where(
        t < 1 / 3,
        np.where(t < -1 / 3, Puiseux_1(t) + c0(t), Taylor0(t) + c1(t)),
        Taylor1(t) + c2(t),
    )


def r2(t):
    return np.where(t < 1 / 3, np.where(t < -1 / 3, r02(t), r12(t)), r22(t))


def r3(t):
    return np.where(t < 1 / 3, np.where(t < -1 / 3, r03(t), r13(t)), r23(t))


def r4(t):
    return np.where(t < 1 / 3, np.where(t < -1 / 3, r04(t), r14(t)), r24(t))


def h(t, k):
    # return np.vectorize(lambda t: newton_(t, g(t), k))(t)
    # return np.vectorize(lambda t: newton_(t, g(t) + r2(t), k))(t)
    return np.vectorize(lambda t: newton_(t, g(t) + r2(t) + r3(t), k))(t)
    # return np.vectorize(lambda t: newton_(t, g(t) + r2(t) + r3(t) + r4(t), k))(t)


def h_(t, k):

    def mul_(t):
        if t < -1 / 3:
            x = x0
            z = z0 + np.sqrt(t + 1) * ww0
        elif t < 1 / 3:
            x = x1
            z = z1
        else:
            x = x2
            z = z2
        for i in range(z.shape[0]):
            if t < x[i + 1]:
                return newton_(
                    t,
                    z[i, 0]
                    + z[i, 1] * t
                    + z[i, 2] * np.power(t, 2)
                    + z[i, 3] * np.power(t, 3)
                    + z[i, 4] * np.power(t, 4),
                    k,
                )
        return 0.0

    return np.vectorize(mul_)(t)


"""
S = interpolate.make_interp_spline(
    x,
    f0(x),
    k=5,
    # bc_type="clamped",
)
plt.plot(t, S(t) - f0(t), label="interp1d")
"""
t = np.arange(-1.0, 1.0, 0.0005)
# plt.plot(t, g(t) - np.cos(np.arccos(t) / 3), label="g")
# plt.plot(t, r2(t), label="corr 2")# plt.plot(t, r2(t) + r3(t), label="corr 3")
# plt.plot(t, r2(t) + r3(t) + r4(t), label="corr 4")
# plt.plot(t, g(t) + r2(t) + r3(t) + r4(t) - np.cos(np.arccos(t) / 3), label="g+r2+r3")
# plt.plot(t, r2(t), label="r2")
# plt.plot(t, r3(t), label="r3")
# plt.plot(t, r4(t), label="r4")
plt.plot(t, h_(t, 0) - np.cos(np.arccos(t) / 3), label="h 0")
plt.plot(t, h_(t, 1) - np.cos(np.arccos(t) / 3), label="h 1")
plt.plot(t, h_(t, 2) - np.cos(np.arccos(t) / 3), label="h 2")
# plt.plot(t, f1(t), label="Taylor0")
# plt.plot(t, f2(t), label="Taylor1")
# plt.plot(t, s0(t) - f0(t), label="S")
# plt.plot(t, corr2(t), label="corr2")
# plt.plot(t, s0(t) + corr2(t) - f0(t), label="S+corr2")
# plt.plot(t, corr3(t), label="corr3")
# plt.plot(t, s0(t) + corr2(t) + corr3(t) - f0(t), label="S+corr2+corr3")
# plt.plot(t, corr4(t), label="corr4")
# plt.plot(t, l0(t) - f0(t), label="linear")
# plt.plot(t, l0(t) + l2(t) - f0(t), label="linear+l2")
# plt.plot(t, l0(t) + l2(t) + l3(t) - f0(t), label="linear+l2+l3")
# plt.plot(t, l0(t) + l2(t) + l3(t) + l4(t) - f0(t), label="linear+l2+l3+l4")
# plt.plot(t, s1(t) + corr2(t) + corr3(t) + corr4(t) - f0(t), label="S+corr2+corr3+corr4")
# plt.stem(x, np.ones(x.size) * np.max(S(t) - f0(t)))
# xm = (x[1:] + x[:-1]) / 2
# plt.stem(xm, np.ones(xm.size) * np.max(S(t) - f0(t)), linefmt=":")
# plt.plot(t, s0(t), label="S")
# plt.plot(t, f0(t), label="f", ls=":")
# plt.plot(t, 4 * np.power(f(t), 3) - 3 * f(t) - t, label="4y(x)^3-3y(x)^3=x")
# plt.ylim([-1.1, 1.1])
# plt.ylim([-1.0e-8, 1.0e-8])
plt.legend()
plt.show()
