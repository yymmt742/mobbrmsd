import numpy as np
import mobbrmsd as mob

print(dir(mob))
mob3 = mob.mobbrmsd()
mob3.add_molecule(5, 3, 2, [3,2,5,4,1])
mob3.add_molecule(7, 2, 1)

n = 1000
x = np.random.rand((5*3+7*2)*3)
y = np.array([0.4*x + 0.6*np.random.rand((5*3+7*2)*3) for i in range(n)]).flatten()

rmsd, upper, lower, log_ratio, neval = mob3.run(x, y, n)
print(rmsd)
print(np.exp(log_ratio))
print(neval)
print()

mob2 = mob.mobbrmsd(d=2)
mob2.add_molecule(7, 6, 2, [2,3,5,6,7,1,4])

n = 1000
x = np.random.rand((7*6)*2)
y = np.array([0.4*x + 0.6*np.random.rand((7*6)*2) for i in range(n)]).flatten()

mob2.maxeval = 100
rmsd, upper, lower, log_ratio, neval = mob3.run(x, y, n)
print(rmsd)
print(np.exp(log_ratio))
print(neval)
print()
