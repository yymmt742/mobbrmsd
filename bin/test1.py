import numpy as np
from symrmsd import driver as sy

a = np.array([3,2,5,4,1])
sy.add_molecule(5, 3, 2, a)
b = np.array([3,6,2,5,4,7,1])
sy.add_molecule(7, 2, 1, b)

n = 1000
x = np.random.rand((5*3+7*2)*3)
y = np.array([0.4*x + 0.6*np.random.rand((5*3+7*2)*3) for i in range(n)]).flatten()

sy.setup()
rmsd, log_ratio, nsearch = sy.run(x, y, n)
print(rmsd)
print(np.exp(log_ratio))
print(nsearch)

sy.clear()

sy.add_molecule(7, 6, 2, b)

n = 1000
x = np.random.rand((7*6)*3)
y = np.array([0.4*x + 0.6*np.random.rand((7*6)*3) for i in range(n)]).flatten()

sy.setup()
rmsd, log_ratio, nsearch = sy.run(x, y, n)
print(rmsd)
print(np.exp(log_ratio))
print(nsearch)

