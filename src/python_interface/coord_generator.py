import numpy as np

d = 3
rng = np.random.default_rng()

def so2():
    g = rng.random(1)
    a = np.cos(2*np.pi*g[0])
    b = np.sin(2*np.pi*g[0])
    return np.array([ ( a,-b ),
                      ( b, a ) ])

def so3():
    g = rng.random(4)
    a = g[1:] / np.linalg.norm(g[1:])
    b = np.sin(2*np.pi*g[0])
    c = np.cos(2*np.pi*g[0])
    s = 1.0 - c
    return np.array([
        [ a[0] * a[0] * s + c,        a[0] * a[1] * s - a[2] * b, a[0] * a[2] * s + a[1] * b ],
        [ a[0] * a[1] * s + a[2] * b, a[1] * a[1] * s + c,        a[1] * a[2] * s - a[0] * b ],
        [ a[0] * a[2] * s - a[1] * b, a[1] * a[2] * s + a[0] * b, a[2] * a[2] * s + c        ] ])

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
  m = 5
  n = 5
  z = 1.5 / np.sqrt(0.1)
  x0 = rng.standard_normal((m, 2))
  x1 = rng.standard_normal((n, m, 2))
  cx = rng.standard_normal((n, 2))
  x2 = np.zeros((n, m, 2))
  so = np.array([so2() for i in range(n)])

  for s in (0.2, 1.0):
    b0 = np.cos(0.5*np.pi*s)
    b1 = np.sin(0.5*np.pi*s)
    for r in (0.0, 1.0):
      a0 = np.sin(0.5*np.pi*r)
      a1 = np.cos(0.5*np.pi*r)
      for i in range(n):
          x2[i] = b0 * cx[i].reshape([1, 2]) + b1 * (a0 * x1[i] + a1 * x0)@so[i]
      sa = '{:d}'.format(round(r*100)).zfill(3)
      sb = '{:d}'.format(round(s*100)).zfill(3)
      path = 'sample_'+sa+'_'+sb+'.npy'
      np.save(path, x2)
