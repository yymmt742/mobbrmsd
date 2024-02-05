import numpy as np
import coord_generator as cg
from symrmsd import driver as sy
import time
import sys

nblock = 1

n = int(sys.argv[1])
a = float(sys.argv[2])
b = float(sys.argv[3])
nrepeat = int(sys.argv[4])
m = 50
s = 1
ss = np.empty(0)

sy.add_molecule(m, n, s, ss)
sy.setup()

def time_block(r, s):
    x = cg.gen(m, n, r, s).flatten()
    y = np.array([cg.gen(m, n, r, s) for i in range(nblock)]).flatten()
    start = time.perf_counter_ns()
    for i in range(nrepeat):
      sy.run(x, y, nblock)
    end = time.perf_counter_ns() - start
    return end / (nrepeat*nblock*1000000000)

print(n, nrepeat, time_block(a, b))

sy.clear()

