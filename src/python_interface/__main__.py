from . import __version__
from . import demo_bb
from . import demo_nn
import sys
import numpy
import time

bar1 = "  ------------------------------------------------------------------------------"
bar2 = " ================================================================================"

print(bar2)
print(
    "                  --- demonstration of mobbrmsd ver.",
    __version__.__version__,
    " ---",
)
print(bar2,"\n")

def run_demo(demo, after=None):

  start_wallclock_time = time.time()
  start_cpu_time = time.process_time()

  ret = demo()

  end_wallclock_time = time.time()
  end_cpu_time = time.process_time()
  elapsed_cpu_time = end_wallclock_time - start_wallclock_time
  elapsed_wallclock_time = end_wallclock_time - start_wallclock_time

  print("    -- Elapsed times --")
  print(
      f"       cpu time : {elapsed_cpu_time:16.3f} sec",
      f" wallclock time : {elapsed_wallclock_time:16.3f} sec"
  )
  print(bar2)

  if after is None: return

  after(ret)

while True:
      inp=input("select demo ['1', '2' or 'q'] >> ")
      if inp=='q':break
      if inp=='1': run_demo(demo_bb.main)
      if inp=='2': run_demo(demo_nn.main, after=demo_nn.show_graph)

