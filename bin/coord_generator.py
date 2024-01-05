import numpy as np

def gen(d, m, n, r, s):
    X0 = np.random.randn(m, d)
    C0 = np.random.randn(n, d)
    X1 = np.random.randn(n, m, d)
    for i in range(n):
        X1[i] = r * C0[i] + (1.0-r) *( s*X1[i] + (1-s)*X0 )
    return X1.reshape([m*n,d])

d=3
m=4
n=5
for r in np.linspace(0.01,0.99):
    for s in np.linspace(0.01,0.99):
        x = gen(d,m,n,r,s)
        print(20,'\n')
        for xi in x:
            print('C ',xi[0]*3, xi[1]*3, xi[2]*3)
