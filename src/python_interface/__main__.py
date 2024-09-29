import argparse
from importlib.metadata import version

__version__ = version(__package__)


def command_demo(args):
    from . import demo_cogen
    from . import demo_bb
    from . import demo_bb_2d
    from . import demo_bb_multi
    from . import demo_batch
    from . import demo_batch_tri
    from . import demo_mst
    import numpy
    import pick

    no = -1 if (args.no is None) else args.no - 1
    prec = numpy.float32 if args.single else numpy.float64
    cli = args.cli

    demo_list = [
        demo_cogen.__demo__(cli=cli, prec=prec),
        demo_bb._demo_bb(cli=cli, prec=prec),
        demo_bb_2d._demo_bb(cli=cli, prec=prec),
        demo_bb_multi.__demo__(cli=cli, prec=prec),
        demo_batch._demo_batch(cli=cli, prec=prec),
        demo_batch_tri._demo_batch_tri(cli=cli, prec=prec),
        demo_mst.__demo__(cli=cli, prec=prec),
    ]

    if (no < 0) or len(demo_list) <= no:
        opts = [l.title for l in demo_list] + ["exit"]
        ver = (
            "--- demonstration of mobbrmsd v."
            + ".".join(__version__.split(".")[:3])
            + " ---"
        )
        bar = "=" * len(ver)
        title = bar + "\n" + ver + "\n" + bar + "\n" + " - Select demo code :"

        while True:
            option, index = pick.pick(opts, title, indicator="   >")
            if option == "exit":
                break
            for l in demo_list:
                if option == l.title:
                    l.run_demo(**args.argv)
                    break
            inp = input('  press any key to continue ("q" to exit) >> ')
            print(inp)
            if inp == "":
                continue
            if inp[0] == "q":
                break
    else:
        demo_list[no].run_demo(**args.argv)


def main():

    class ParamProc(argparse.Action):
        def __call__(self, parser, namespace, values, option_strings=None):
            prm_dict = getattr(namespace, self.dest, [])
            if prm_dict is None:
                prm_dict = {}
            for v in values:
                s = v.split("=", maxsplit=1)
                if len(s) == 2:
                    prm_dict[s[0]] = s[1]
            setattr(namespace, self.dest, prm_dict)

    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    parser_demo = sub.add_parser("demo", help="run demo codes")
    parser_demo.add_argument("--no", type=int, help="demo id.")
    parser_demo.add_argument(
        "-c", "--cli", action="store_true", help="CLI interaction mode."
    )
    parser_demo.add_argument(
        "-s", "--single", action="store_true", help="Use single precision."
    )
    parser_demo.add_argument("argv", nargs="*", action=ParamProc)
    parser_demo.set_defaults(handler=command_demo)
    args = parser.parse_args()
    if hasattr(args, "handler"):
        args.handler(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
