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


def main():
    demo_list = [
        demo_bb._demo_bb(),
        demo_bb_multi._demo_bb_multi(),
        demo_cogen._demo_cogen(),
        demo_batch._demo_batch(),
        demo_batch_tri._demo_batch_tri(),
        demo_mst._demo_mst(),
    ]
    opts = [l.title for l in demo_list] + ["exit"]

    ver = "--- demonstration of mobbrmsd ver." + __version__.__version__ + " ---"
    bar = "=" * len(ver)
    title = bar + "\n" + ver + "\n" + bar + "\n" + " - Select demo code :"

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


if __name__ == "__main__":
    main()
