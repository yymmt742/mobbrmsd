import numpy as np
import coord_generator as cg
from symrmsd import driver as sy
import sys

nblock = 20
nrepeat = 500

m = 50
n = 12
s = 1
ss = np.empty(0)

sy.add_molecule(m, n, s, ss)
sy.setup()

def test_block(r, s):
    x = cg.gen(m, n, r, s).flatten()
    y = np.array([cg.gen(m, n, r, s) for i in range(nblock)]).flatten()
    return sy.run(x, y, nblock)

def sample(a, b):
    sa = '{:d}'.format(round(a*100)).zfill(3)
    sb = '{:d}'.format(round(b*100)).zfill(3)
    path = 'sample_'+sa+'_'+sb+'.dat'
    with open(path, 'w') as f:
        for i in range(nrepeat):
          if i%10==0: print('  repeat',i)
          t = test_block(a, b)
          for z in zip(t[0], t[1], t[2]):
            f.write('{:24.9f} {:24.9f} {:24d}\n'.format(z[0], z[1], z[2]))

for a in np.linspace(0.0, 1.0, 11):
   for b in np.linspace(0.2, 0.4, 3):
       print(a, b)
       sample(a, b)

sy.clear()
