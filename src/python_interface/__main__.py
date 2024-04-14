from . import __version__
from . import demo_bb
import sys
import numpy
import time

start_wallclock_time = time.time()
start_cpu_time = time.process_time()

print(
    "  ==============================================================================="
)
print(
    "                  --- demonstration of mobbrmsd ver.",
    __version__.__version__,
    " ---",
)
print(
    "  ===============================================================================\n"
)

demo_bb.main(n_apm=3, n_mol=10, sym=[[1, 2, 0], [2, 0, 1]], a=0.5, b=1.0)

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
