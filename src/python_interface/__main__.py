from . import __version__
from . import demo_cogen
from . import demo_bb
from . import demo_bb_multi
from . import demo_batch
from . import demo_batch_tri
from . import demo_mst
import numpy
import pick


def main():
    demo_list = [
        demo_cogen._demo_cogen(),
        demo_cogen._demo_cogen(prec=numpy.float32),
        demo_bb._demo_bb(),
        demo_bb._demo_bb(prec=numpy.float32),
        demo_bb_multi._demo_bb_multi(),
        demo_bb_multi._demo_bb_multi(prec=numpy.float32),
        demo_batch._demo_batch(),
        demo_batch._demo_batch(prec=numpy.float32),
        demo_batch_tri._demo_batch_tri(),
        demo_batch_tri._demo_batch_tri(prec=numpy.float32),
        demo_mst._demo_mst(),
        demo_mst._demo_mst(prec=numpy.float32),
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
