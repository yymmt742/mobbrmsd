from . import coord_generator
from . import mobbrmsd

cogen = coord_generator()
n_apm = 3
n_mol = 5
a = 0.8
b = 1.0

x = cogen.generate(n_apm, n_mol, a, b).reshape([n_apm * n_mol, 3])
y = cogen.generate(n_apm, n_mol, a, b).reshape([n_apm * n_mol, 3])
z = y.copy()
print(z - y)

print("  Reference coordinate\n")
print(x.reshape([n_apm, n_mol, 3]))
print()
print("  Target coordinate\n")
print(y.reshape([n_apm, n_mol, 3]))
print()

mrmsd = mobbrmsd()
mrmsd.add_molecule(n_apm, n_mol, [[2, 3, 1], [3, 1, 2]])

ret = mrmsd.run(x, y, maxeval=0)
print(ret.rmsd, ret.bounds, ret.n_eval, ret.eval_ratio)
while not ret.is_finished:
    ret = mrmsd.restart(ret, maxeval=0)
    print(ret.rmsd, ret.bounds, ret.n_eval, ret.eval_ratio)
ret, _ = mrmsd.restart(ret, maxeval=0, Y=y)
print(z - y)
import numpy

print(numpy.sqrt(numpy.sum(numpy.power(x - y, 2)) / (n_apm * n_mol)))
