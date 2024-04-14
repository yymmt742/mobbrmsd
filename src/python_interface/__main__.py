from . import coord_generator
from . import mobbrmsd
import numpy
import time

start_wallclock_time = time.time()
start_cpu_time = time.process_time()

cogen = coord_generator()
n_apm = 3
n_mol = 8
sym = [[1, 2, 0], [2, 0, 1]]
a = 0.5
b = 1.0

x = cogen.generate(n_apm, n_mol, a, b).reshape([n_apm * n_mol, 3])
y = cogen.generate(n_apm, n_mol, a, b).reshape([n_apm * n_mol, 3])
x -= numpy.mean(x, 0)
y -= numpy.mean(y, 0)
z = y.copy()

print(
    "  ==============================================================================="
)
print("                      --- demonstration of mobbrmsd ---")
print(
    "  ===============================================================================\n"
)

print("      --System settings--")
print(f"    Atoms per molecule :{n_apm:2d}")
print(f"    Number of molecule :{n_mol:2d}")
print("    molecular symmetry : 0 ", [i for i in range(3)])
for i in range(len(sym)):
    print(f"                       :{i+1:2d} ", [si for si in sym[i]])
print()

mrmsd = mobbrmsd()
mrmsd.add_molecule(n_apm, n_mol, sym)

ub, lb = numpy.inf, 0.0
ret = mrmsd.run(x, y, maxeval=0)

print("    --------------------------------------------------------------------------")
print("                Molecular-oriented RMSD for Branch-and-bound")
print("    --------------------------------------------------------------------------")
print("        N_eval   Eval_ratio      Upperbound      Lowerbound          RMSD")
while not ret.is_finished:
    if ub > ret.bounds[0] or lb < ret.bounds[1]:
        print(
            f"  {ret.n_eval:12d} {ret.eval_ratio:12.8f}{ret.bounds[0]:16.6f}{ret.bounds[1]:16.6f}{ret.rmsd:16.6f} "
        )
    ub, lb = ret.bounds[0], ret.bounds[1]
    ret = mrmsd.restart(ret, maxeval=0)
ret = mrmsd.restart(ret, maxeval=0, Y=y)
print("    --------------------------------------------------------------------------")
print("      -- Final results --")
print("        N_eval   Eval_ratio      Upperbound      Lowerbound          RMSD")
print(
    f"  {ret.n_eval:12d} {ret.eval_ratio:12.8f}{ret.bounds[0]:16.6f}{ret.bounds[1]:16.6f}{ret.rmsd:16.6f} "
)
print(
    "    --------------------------------------------------------------------------\n"
)

print(
    "  ------------------------------------------------------------------------------"
)
print("        Reference     | Target (original) |  disp. |  Target (rotate)  |  disp.")
for xi, yi, zi in zip(
    x.reshape([n_mol, n_apm, 3]),
    y.reshape([n_mol, n_apm, 3]),
    z.reshape([n_mol, n_apm, 3]),
):
    print(
        "  --------------------|-------------------|--------|-------------------|--------"
    )
    for xij, yij, zij in zip(xi, yi, zi):
        d1 = numpy.sum(numpy.power(xij - zij, 2))
        d2 = numpy.sum(numpy.power(xij - yij, 2))
        print(
            f"  {xij[0]:6.2f}{xij[1]:6.2f}{xij[2]:6.2f} ",
            f"|{yij[0]:6.2f}{yij[1]:6.2f}{yij[2]:6.2f} |{d1:7.2f}",
            f"|{zij[0]:6.2f}{zij[1]:6.2f}{zij[2]:6.2f} |{d2:7.2f}",
        )
d1, d2 = numpy.sum(numpy.power(x - z, 2)), numpy.sum(numpy.power(x - y, 2))
print(
    "  ----------------------------------------|--------|-------------------|---------"
)
print(
    f"              squared deviation           |{d1:7.2f} |                   |{d2:7.2f}"
)
d1 /= n_apm * n_mol
d2 /= n_apm * n_mol
print(
    f"           mean squared deviation         |{d1:7.2f} |                   |{d2:7.2f}"
)
d1, d2 = numpy.sqrt(d1), numpy.sqrt(d2)
print(
    f"         root mean squared deviation      |{d1:7.2f} |                   |{d2:7.2f}"
)
print(
    "  -------------------------------------------------------------------------------"
)

end_wallclock_time = time.time()
end_cpu_time = time.process_time()
elapsed_cpu_time = end_wallclock_time - start_wallclock_time
elapsed_wallclock_time = end_wallclock_time - start_wallclock_time

print("    -- Elapsed times --")
print(
    f"    cpu time : {elapsed_cpu_time:.6f} sec   wallclock time : {elapsed_wallclock_time:.6f} sec"
)
print(
    "  ===============================================================================\n"
)
