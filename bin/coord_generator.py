import numpy as np

d = 3
rng = np.random.default_rng()

def so3():
    g = rng.random(3)
    a = g / np.linalg.norm(g)
    return np.array([
        [ a[0] * a[0],        a[0] * a[1] - a[2], a[0] * a[2] + a[1] ],
        [ a[0] * a[1] + a[2], a[1] * a[1],        a[1] * a[2] - a[0] ],
        [ a[0] * a[2] - a[1], a[1] * a[2] + a[0], a[2] * a[2]        ] ])

def mol(m, n, s):
    X0 = rng.standard_normal((m, d))
    X1 = rng.standard_normal((n, m, d))
    return X1

def gen(m, n, a, b):

    a0 = np.sin(0.5*np.pi*a)
    a1 = np.cos(0.5*np.pi*a)
    b0 = np.sin(0.5*np.pi*b)
    b1 = np.cos(0.5*np.pi*b)

    x0 = rng.standard_normal((m, d))
    x1 = rng.standard_normal((n, m, d))
    cx = rng.standard_normal((n, d))

    for i in range(n):
        x1[i] = b0 * cx[i].reshape([1, d]) + b1 * (a0 * x1[i] + a1 * x0)@so3()
    x1 -= np.mean(x1.reshape((n * m, d)), 0).reshape([1, d])

    return x1.reshape([n * m, d])

if (__name__=='__main__'):
    for r in np.linspace(0.0, 1.0, 10):
        for s in np.linspace(0.0, 1.0, 10):
            x = gen_pair(4, 5, r, s)
            print(20,'\n')
            for xi in x:
                print('C ',xi[0]*3, xi[1]*3, xi[2]*3)
