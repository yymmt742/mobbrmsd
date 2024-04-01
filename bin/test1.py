import numpy as np
import matplotlib.pyplot as plt
import mobbrmsd as mob

print(dir(mob))
mob3 = mob.mobbrmsd()
mob3.add_molecule(5, 3, [3,2,5,4,1])
mob3.add_molecule(7, 2)

n = 1000
x = np.random.rand((5*3+7*2),3)
y = np.array([0.4*x + 0.6*np.random.rand((5*3+7*2),3) for i in range(n)])
print(x.shape)
print(y.shape)

ret = mob3.run(x, y)
#plt.hist(mob3.rmsd(ret['state']))
#plt.show()
print(mob3.bounds(ret['state'][0:2]))
print(mob3.n_eval(ret['state'][-1:0:-1]))
print(mob3.eval_ratio(ret['state'][9]))
print()
exit()

mob2 = mob.mobbrmsd(d=2)
mob2.add_molecule(7, 6, [2,3,5,6,7,1,4])

n = 1000
x = np.random.rand((7*6)*2)
y = np.array([0.4*x + 0.6*np.random.rand((7*6)*2) for i in range(n)])

ints, floats = mob3.run(x, y, n)
print(ints[0:3])
print(floats[0:3])
print(ints.flags['C_CONTIGUOUS'])
print()
