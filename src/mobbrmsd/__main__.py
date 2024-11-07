import argparse
from .dataclass import molecules, molecular_system, load
from ._mobbrmsd import mobbrmsd
from importlib.metadata import version

__version__ = version(__package__)


def command_run(args):
    import numpy
    import json
    import mdtraj
    from pathlib import Path

    prms = {**args.params}
    if args.inp is not None:
        try:
            with open(args.inp[0], "r") as f:
                try:
                    prms = {**json.load(f), **prms}
                except:
                    print(f"Warning : load json [{args.inp[0]}] is failed.")
                    pass
        except:
            print(f"Warning : open json [{args.inp[0]}] is failed.")
            pass

    def parse_coordinate(obj, top, dtype):

        if obj is None:
            return None

        if isinstance(obj, str):
            try:
                ref = mdtraj.load(obj, topology=top).xyz * 10.0
                if dtype != numpy.float32:
                    ref = numpy.array(ref, dtype=dtype)
            except:
                try:
                    ref = numpy.array(json.loads(obj), dtype=dtype)
                except:
                    raise IOError
        else:
            try:
                ref = numpy.array(obj, dtype=dtype)
            except:
                raise IOError
        return ref

    top = prms.get("top")
    dtype = numpy.float32 if prms.get("precision") == "single" else numpy.float64
    ref = parse_coordinate(prms.get("reference"), top, dtype)
    trg = parse_coordinate(prms.get("target"), top, dtype)
    if ref is None:
        raise IOError

    mols = load(prms.get("molecules"))

    mrmsd = mobbrmsd(mols=mols)
    if trg is None:
        ret = mrmsd.batch_run(ref, **prms)
    else:
        ret = mrmsd.batch_run(
            ref,
            trg,
            **prms,
        )
    print(ret)


def command_demo(args):
    import numpy
    import pick
    from .demo import cogen
    from .demo import bb
    from .demo import bb_2d
    from .demo import bb_multi
    from .demo import batch
    from .demo import batch_tri
    from .demo import mst

    no = -1 if (args.no is None) else args.no - 1
    prec = numpy.float32 if args.single else numpy.float64
    cli = args.cli

    demo_list = [
        cogen.__demo(cli=cli, prec=prec),
        bb.__demo(cli=cli, prec=prec),
        bb_2d.__demo(cli=cli, prec=prec),
        bb_multi.__demo(cli=cli, prec=prec),
        batch.__demo(cli=cli, prec=prec),
        batch_tri.__demo(cli=cli, prec=prec),
        mst.__demo(cli=cli, prec=prec),
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

    parser_run = sub.add_parser("run", help="run mobbrmsd")
    parser_run.add_argument(
        "-i",
        "--inp",
        nargs=1,
        default=None,
        help="Input file (json)",
    )
    parser_run.add_argument(
        "params",
        nargs="*",
        action=ParamProc,
        help="Keyword arguments.",
    )
    parser_run.set_defaults(handler=command_run)

    parser_demo = sub.add_parser("demo", help="compute demo codes")
    parser_demo.add_argument(
        "--no", type=int, help="demo id", choices=[1, 2, 3, 4, 5, 6, 7]
    )
    parser_demo.add_argument(
        "-c", "--cli", action="store_true", help="CLI interaction mode."
    )
    parser_demo.add_argument(
        "-s", "--single", action="store_true", help="Use single precision"
    )
    parser_demo.add_argument(
        "argv",
        nargs="*",
        action=ParamProc,
        help="demo-specific keyword arguments. Example: hoge=1 fuga=a,b,c",
    )
    parser_demo.set_defaults(handler=command_demo)

    args = parser.parse_args()
    if hasattr(args, "handler"):
        args.handler(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
