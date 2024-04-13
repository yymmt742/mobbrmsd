from . import coord_generator
from . import mobbrmsd

cogen = coord_generator()
n = 8
m = 5
a = 0.8
b = 1.0

x = cogen.generate(n, m, a, b)
y = cogen.generate(n, m, a, b)

print("  Reference coordinate")
print(x)
print("  Target coordinate")
print(y)

mrmsd = mobbrmsd()
mrmsd.add_molecule(n, m)

ret = mrmsd.run(x.reshape([m * n, 3]), y.reshape([m * n, 3]), maxeval=0)
print(ret.bounds, ret.n_eval, ret.eval_ratio)
ret = mrmsd.restart(ret, maxeval=0)
print(ret.bounds, ret.n_eval, ret.eval_ratio)
ret = mrmsd.restart(ret, maxeval=0)
print(ret.bounds, ret.n_eval, ret.eval_ratio)
