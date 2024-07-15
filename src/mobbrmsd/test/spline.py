import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


def spline_forward(f0, f1, f2, x0, x1, a0):
    dx = x1 - x0
    A = (f0(x1) - f0(x0)) / (dx * dx) - f1(x0) / dx
    B = (f1(x1) - f1(x0)) / dx
    y0 = f0(x0)
    y1 = f1(x0)
    a1 = 12 * A - 6 * B - a0
    # a1 = -12 * A + 6 * B + a0
    b0 = (24 * A - 6 * B - 6 * a0) / dx
    c0 = (-72 * A + 24 * B + 12 * a0) / (dx * dx)
    print(f"y0 {a0:9.3f} {b0:9.3f} {c0:9.3f} {a1:9.3f} ")
    return (
        a1,
        lambda x: (
            y0
            + y1 * (x - x0)
            + (a0 / 2) * np.power(x - x0, 2)
            + (b0 / 6) * np.power(x - x0, 3)
            + (c0 / 24) * np.power(x - x0, 4)
        ),
    )


f0 = lambda x: np.cos(np.arccos(x) / 3.0) - (np.sqrt(1 + x)) / 2 - 0.5
f1 = lambda x: (
    (np.sin(np.arccos(x) / 3) / (3 * np.sqrt(1.0 - x)) - 0.25) / np.sqrt(x + 1.0)
)
f2 = lambda x: (
    -np.cos(np.arccos(x) / 3) / (9 * (1 - np.power(x, 2)))
    + (x * np.sin(np.arccos(x) / 3) / (3 * np.power(np.sqrt(1 - x), 3)) + 0.125)
    / np.power(np.sqrt(1 + x), 3)
)
x = np.array([-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0])
y = f0(x)
# S = interpolate.interp1d(x, y, kind="cubic")
c = interpolate.PPoly.from_spline(interpolate.splrep(x, y)).c.T
print(c)

t = np.arange(-0.99, 0.99, 0.01)
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
