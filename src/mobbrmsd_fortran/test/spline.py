import numpy as np
import matplotlib.pyplot as plt

"""
# from scipy import interpolate
S = interpolate.make_interp_spline(
    x,
    f0(x),
    k=5,
    # bc_type="clamped",
)
plt.plot(t, S(t) - f0(t), label="interp1d")
"""


# expand cosh(log(|x| + SQRT(x^2 -1))/3)
def Taylor1_cosh(x):
    return (
        1
        + (x - 1) / 9
        - 4 / 243 * np.power(x - 1, 2)
        + (28 * np.power(x - 1, 3)) / 6561
        - (80 * np.power(x - 1, 4)) / 59049
    )


def Taylor1_rcosh(x):
    return (
        np.sqrt(3) / 2
        + x / 6
        - np.power(x, 2) / (12 * np.sqrt(3))
        + 2 * np.power(x, 3) / 81
        - 35 * np.power(x, 4) / (1296 * np.sqrt(3))
    )


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


# newton np.cos(np.arccos(x) / 3.0)
def newton_(x, y0, k):
    if x <= -1.0:
        return 1 / 2
    i = 0
    y = y0
    s = 1.0
    while i < k and np.abs(s) > 1e-17:
        yy = y * y
        f = (4 * yy - 3) * y - x
        df = 12 * y**2 - 3
        df = np.sign(df) * np.max([np.abs(df), 1e-18])
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


n = 4 + 1
x0 = np.cos(np.pi * np.linspace(1, 0, n)) / 3 - 2 / 3
x1 = np.cos(np.pi * np.linspace(1, 0, n)) / 3
x2 = np.cos(np.pi * np.linspace(1, 0, n)) / 3 + 2 / 3
n = 4 + 1
x3 = np.cos(np.pi * np.linspace(1, 0, n)) / 2 + 1 / 2

f0 = np.vectorize(
    lambda x: (np.cos(np.arccos(x) / 3) - Puiseux_1(x) if x > -1.0 else 0.0)
)
f1 = np.vectorize(lambda x: np.cos(np.arccos(x) / 3) - Taylor0(x))
f2 = np.vectorize(lambda x: np.cos(np.arccos(x) / 3) - Taylor1(x))

f3 = np.vectorize(
    lambda x: (
        0.0
        if x <= 0.0
        else np.cosh(np.arccosh(1 / x) / 3)
        - (np.power(2 / x, 1 / 3) + 1 / np.power(2 / x, 1 / 3)) / 2
    )
)

# c0 = spline_linear(x0, f0)
# c1 = spline_linear(x1, f1)
# c2 = spline_linear(x2, f2)
c0, z0 = spline_cubic(x0, f0)
c1, z1 = spline_cubic(x1, f1)
c2, z2 = spline_cubic(x2, f2)
c3, z3 = spline_cubic(x3, f3)
r02, r03, r04, v0 = corr(x0, f0, c0)
r12, r13, r14, v1 = corr(x1, f1, c1)
r22, r23, r24, v2 = corr(x2, f2, c2)
r32, r33, r34, v3 = corr(x3, f3, c3)

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
z3 += v3


def h(t, k):

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


def hc(t, k):

    def mul_(t):
        if t < 0.0:
            return 0.0
        for i in range(z3.shape[0]):
            if t < x3[i + 1]:
                return newton_(
                    1 / t,
                    (np.power(2 / t, 1 / 3) + np.power(t / 2, 1 / 3)) / 2
                    + z3[i, 0]
                    + z3[i, 1] * t
                    + z3[i, 2] * np.power(t, 2)
                    + z3[i, 3] * np.power(t, 3)
                    + z3[i, 4] * np.power(t, 4),
                    k,
                )
        return 1.0

    return np.vectorize(mul_)(t)


"""
t = np.arange(-1.0, 1.0, 0.0005)
plt.plot(t, h(t, 0) - np.cos(np.arccos(t) / 3), label="h 0")
plt.plot(t, h(t, 1) - np.cos(np.arccos(t) / 3), label="h 1")
plt.plot(t, h(t, 2) - np.cos(np.arccos(t) / 3), label="h 2")
plt.plot(t, f1(t), label="Taylor0")
plt.plot(t, f2(t), label="Taylor1")
plt.stem(x, np.ones(x.size) * np.max(S(t) - f0(t)))
xm = (x[1:] + x[:-1]) / 2
plt.stem(xm, np.ones(xm.size) * np.max(S(t) - f0(t)), linefmt=":")
plt.plot(t, s0(t), label="S")
plt.plot(t, f0(t), label="f", ls=":")
plt.plot(t, 4 * np.power(f(t), 3) - 3 * f(t) - t, label="4y(x)^3-3y(x)^3=x")
plt.ylim([-1.1, 1.1])
plt.ylim([-1.0e-8, 1.0e-8])
s = np.arange(0.1, 1.0, 0.001)
t = np.arange(1.0, 5.0, 0.001)
nf = hc(s, 1)
plt.plot(
    s,
    4 * np.power(nf, 3) - 3 * nf - 1 / s,
    label="spline n1",
)
nf = hc(s, 2)
plt.plot(
    s,
    4 * np.power(nf, 3) - 3 * nf - 1 / s,
    label="spline n2",
    ls=":",
)
plt.legend()
plt.show()
"""

xl = 0.0
xu = 1.0
print("#:set cos1_sqrt_part_x = [ &")
print(f'&   "{xl:.16e}", "{xu:.16e}", &')
print("& ]")
print("#:set cos1_sqrt_part_w = [ &")
print(
    f'&   "{ww0[4]:.16e}", "{ww0[3]:.16e}", "{ww0[2]:.16e}", "{ww0[1]:.16e}", "{ww0[0]:.16e}", &'
)
print("& ]")

for x, z, t in zip(
    [x0, x1, x2, x3],
    [z0, z1, z2, z3],
    [
        "cos1",
        "cos2",
        "cos3",
        "cosh",
    ],
):
    print("#:set " + t + "_x = [ &")
    for xl, xu in zip(x[:-1], x[1:]):
        print(f'&   ["{xl:.16e}", "{xu:.16e}",], &')
    print("& ]")
    print("#:set " + t + "_w = [ &")
    for zi in z:
        print(
            f'&   ["{zi[4]:.16e}", "{zi[3]:.16e}", "{zi[2]:.16e}", "{zi[1]:.16e}", "{zi[0]:.16e}",], &'
        )
    print("& ]")
