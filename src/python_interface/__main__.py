from . import __version__
from . import demo_cogen
from . import demo_bb
from . import demo_bb_sp
from . import demo_bb_multi
from . import demo_batch
from . import demo_batch_tri
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

demo_list = [
    demo_bb._demo_bb(),
    demo_cogen._demo_cogen(),
]
demo1 = demo_cogen.title
demo2 = demo_bb.title
demo3 = demo_bb_sp.title
demo4 = demo_bb_multi.title
demo5 = demo_batch.title
demo6 = demo_batch_tri.title
demo7 = demo_mst.title
exit_ = "exit"
# opts = (demo1, demo2, demo3, demo4, demo5, demo6, demo7, exit_)
opts = [l.title for l in demo_list] + ["exit"]


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
    if option == "exit":
        break
    for l in demo_list:
        if option == l.title:
            l.run_demo()
            break
    inp = input('  press any key to continue ("q" to exit) >> ')
    print(inp)
    if inp == "":
        continue
    if inp[0] == "q":
        break
