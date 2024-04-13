from . import coord_generator
from . import mobbrmsd
import numpy


def show_coordinates(n_apm, n_mol, x, y):
    print("  ------------------------------------------------------------")
    print("   Reference coordinate |  Target coordinate   | displacement")
    for xi, yi in zip(x.reshape([n_mol, n_apm, 3]), y.reshape([n_mol, n_apm, 3])):
        print("  ----------------------|----------------------|-------------")
        for xij, yij in zip(xi, yi):
            d = numpy.sum(numpy.power(xij - yij, 2))
            print(
                f"  {xij[0]:7.3f}{xij[1]:7.3f}{xij[2]:7.3f} |{yij[0]:7.3f}{yij[1]:7.3f}{yij[2]:7.3f} |{d:9.3f}",
            )
    d = numpy.sum(numpy.power(x - y, 2))
    print("  ---------------------------------------------|-------------")
    print(f"                squared deviation              |{d:9.3f}    ")
    d /= n_apm * n_mol
    print(f"             mean squared deviation            |{d:9.3f}    ")
    d = numpy.sqrt(d)
    print(f"           root mean squared deviation         |{d:9.3f}    ")
    print("  -----------------------------------------------------------\n")


cogen = coord_generator()
n_apm = 3
n_mol = 4
a = 0.8
b = 1.0

x = cogen.generate(n_apm, n_mol, a, b).reshape([n_apm * n_mol, 3])
y = cogen.generate(n_apm, n_mol, a, b).reshape([n_apm * n_mol, 3])
x -= numpy.mean(x, 0)
y -= numpy.mean(y, 0)

print("  ============================================================")
print("              --- demonstration of mobbrmsd ---")
print("  ============================================================\n")

show_coordinates(n_apm, n_mol, x, y)

mrmsd = mobbrmsd()
mrmsd.add_molecule(n_apm, n_mol, [[2, 3, 1], [3, 1, 2]])

ub, lb = numpy.inf, 0.0
ret = mrmsd.run(x, y, maxeval=0)

print("  ----------------------------------------------------------")
print("         Branch-and-bound for Molecular-oriented RMSD")
print("    N_eval  Coverage  Upperbound  Lowerbound      RMSD")
print("  ----------------------------------------------------------")
while not ret.is_finished:
    if ub > ret.bounds[0] or lb < ret.bounds[1]:
        print(
            f"  {ret.n_eval:8d} {ret.eval_ratio:9.3f}{ret.bounds[0]:12.6f}{ret.bounds[1]:12.6f}{ret.rmsd:12.6f} "
        )
    ub, lb = ret.bounds[0], ret.bounds[1]
    ret = mrmsd.restart(ret, maxeval=0)
ret = mrmsd.restart(ret, maxeval=0, Y=y)
print("  ----------------------------------------------------------")
print("      Final results")
print(
    f"  {ret.n_eval:8d} {ret.eval_ratio:9.3f}{ret.bounds[0]:12.6f}{ret.bounds[1]:12.6f}{ret.rmsd:12.6f}"
)
print("  ----------------------------------------------------------\n")

show_coordinates(n_apm, n_mol, x, y)

print("              End of computation")
print("  ============================================================\n")
