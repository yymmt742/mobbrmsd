import numpy as np
import mobbrmsd as sy

print(dir(sy))
sy.add_molecule(5, 3, 2, [3,2,5,4,1])
sy.add_molecule(7, 2, 1)

n = 1000
x = np.random.rand((5*3+7*2)*3)
y = np.array([0.4*x + 0.6*np.random.rand((5*3+7*2)*3) for i in range(n)]).flatten()

rmsd, upper, lower, log_ratio, neval = sy.run(x, y, n)
print(rmsd)
print(np.exp(log_ratio))
print(neval)
exit()

sy.add_molecule(7, 6, 2, [2,3,5,6,7,1,4])

n = 1000
x = np.random.rand((7*6)*3)
y = np.array([0.4*x + 0.6*np.random.rand((7*6)*3) for i in range(n)]).flatten()

sy.maxeval = 100
rmsd, log_ratio, nsearch, rmsd_we = sy.run(x, y, n)
print(rmsd)
print(np.exp(log_ratio))
print(nsearch)

rmsd, log_ratio, nsearch, rmsd_we = sy.run(x, y, n, )
