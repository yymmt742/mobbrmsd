from . import __version__
from . import demo_bb
from . import demo_mst
import sys
import numpy
import time
import pick

bar1 = (
    "  ------------------------------------------------------------------------------"
)
bar2 = (
    "  =============================================================================="
)
bar3 = "==========================================="
demo1 = demo_bb.title
demo2 = demo_mst.title
exit_ = "exit"
opts = (demo1, demo2, exit_)


def run_demo(demo, read_input, after=None):

    print(bar2)

    prms = read_input()

    start_wallclock_time = time.time()
    start_cpu_time = time.process_time()

    ret = demo(**prms)

    end_wallclock_time = time.time()
    end_cpu_time = time.process_time()
    elapsed_cpu_time = end_wallclock_time - start_wallclock_time
    elapsed_wallclock_time = end_wallclock_time - start_wallclock_time

    print("    -- Elapsed times --")
    print(
        f"       cpu time : {elapsed_cpu_time:16.3f} sec",
        f" wallclock time : {elapsed_wallclock_time:16.3f} sec",
    )
    print(bar2, "\n")

    if after is not None:
        after(ret)
        print(bar2, "\n")


title = (
    bar3
    + "\n"
    + "--- demonstration of mobbrmsd ver."
    + __version__.__version__
    + " ---\n"
    + bar3
    + "\n"
    + "    Select demo code :"
)

while True:
    option, index = pick.pick(opts, title, indicator="   >")
    if option == demo1:
        run_demo(demo_bb.main, demo_bb.read_input)
    elif option == demo2:
        run_demo(demo_mst.main, demo_mst.read_input, after=demo_mst.show_graph)
    elif option == exit_:
        break
    inp = input('  press any key to continue ("q" to exit) >> ')
    print(inp)
    if inp == "":
        continue
    if inp[0] == "q":
        break
