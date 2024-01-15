import numpy as np
import coord_generator as cg
from symrmsd import driver as sy
import sys

nblock = 20
nrepeat = 100

m = 50
n = 12
s = 1
ss = np.empty(0)

sy.add_molecule(m, n, s, ss)
sy.setup()

def sample(r, s, ns):
    return

def test_block(r, s):
    x = cg.gen(m, n, r, s).flatten()
    y = np.array([cg.gen(m, n, r, s) for i in range(nblock)]).flatten()
    return sy.run(x, y, nblock)

def sample(a, b):
    sa = '{:d}'.format(int(a*100)).zfill(3)
    sb = '{:d}'.format(int(b*100)).zfill(3)
    path = 'sample_'+sa+'_'+sb+'.dat'
    with open(path, 'w') as f:
        for i in range(nrepeat):
            t = test_block(a, b)
            for z in zip(t[0], t[1], t[2]):
                f.write('{:24.9f} {:24.9f} {:24d}\n'.format(z[0], z[1], z[2]))

sample(0.70, 1.00)
for a in np.linspace(0.8, 1.0, 3):
   for b in np.linspace(0.5, 1.0, 6):
       print(a, b)
       sample(a, b)

sy.clear()
