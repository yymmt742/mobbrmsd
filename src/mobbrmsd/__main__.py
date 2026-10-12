import argparse
import sys
from .dataclass import molecules, molecular_system, load
from importlib.metadata import version

__version__ = version(__package__)

_err = sys.stderr


def _print_species_table(ms):
    print(file=_err)
    print(f"  Discovered {len(ms)} species", file=_err)
    w_name = max((len(m.get("name") or "—") for m in ms), default=4)
    w_name = max(w_name, 4)
    hdr = f"  {'#':>3}  {'Name':<{w_name}}  {'n_apm':>5}  {'n_mol':>5}  sym"
    print(hdr, file=_err)
    print(f"  {'─'*3}  {'─'*w_name}  {'─'*5}  {'─'*5}  {'─'*3}", file=_err)
    for i, m in enumerate(ms, 1):
        name = m.get("name") or "—"
        sym = m.get("sym")
        sym_str = f"{len(sym)} perm" if sym else "—"
        print(f"  {i:>3}  {name:<{w_name}}  {m['n_apm']:>5}  {m['n_mol']:>5}  {sym_str}", file=_err)
    total = sum(m["n_apm"] * m["n_mol"] for m in ms)
    print(f"  Total atoms: {total}", file=_err)


def _write_json(mols, output):
    import json
    d = mols.to_dict()
    data = (json.dumps(d, indent=2) + "\n").encode("utf-8")
    if output:
        with open(output, "wb") as f:
            f.write(data)
        print(f"  Written to {output}", file=_err)
    else:
        sys.stdout.buffer.write(data)
        sys.stdout.buffer.flush()


def command_build(args) -> None:
    if getattr(args, "file", None):
        try:
            import mdtraj
        except ImportError:
            print("Error: mdtraj is required. Install with: pip install mdtraj", file=_err)
            return

        top_arg = getattr(args, "topology", None)
        print(f"  Parsing topology from {top_arg or args.file} …", file=_err)
        try:
            if top_arg:
                top = mdtraj.load(args.file, top=top_arg).topology
            else:
                top = mdtraj.load_topology(args.file)
        except Exception as e:
            print(f"  Error loading topology: {e}", file=_err)
            return

        from .dataclass import molecular_system as ms_cls
        try:
            mols = ms_cls.from_topology(top)
        except Exception as e:
            print(f"  Error: {e}", file=_err)
            print("  Launch topology builder to define species manually? [Y/n] ", end="", flush=True, file=_err)
            try:
                ans = input().strip().lower()
            except (EOFError, KeyboardInterrupt):
                print("\n  Cancelled.", file=_err)
                return
            if ans in ("n", "no"):
                return
            from .run.topology_builder import TopologyBuilder
            builder = TopologyBuilder()
            builder.run()
            try:
                mols = builder.build()
            except ValueError as e2:
                print(f"  Error: {e2}", file=_err)
                return
            _write_json(mols, args.output)
            return

        ms_dict = mols.to_dict().get("molecular_system", {})
        species = ms_dict.get("molecules", [])
        if isinstance(species, dict):
            species = [species]
        _print_species_table(species)

        print(file=_err)
        print("  Proceed? [Y]es / [n]o / [e]dit > ", end="", flush=True, file=_err)
        try:
            ans = input().strip().lower()
        except (EOFError, KeyboardInterrupt):
            print("\n  Cancelled.", file=_err)
            return

        if ans in ("n", "no"):
            print("  Cancelled.", file=_err)
            return
        elif ans in ("e", "edit"):
            from .run.parse import builder_from_topology
            print(f"  Launching topology builder ({len(species)} species pre-loaded) …", file=_err)
            print(file=_err)
            builder = builder_from_topology(top)
            builder.run()
            try:
                mols = builder.build()
            except ValueError as e:
                print(f"  Error: {e}", file=_err)
                return

    else:
        from .run.topology_builder import TopologyBuilder
        builder = TopologyBuilder()
        builder.run()
        try:
            mols = builder.build()
        except ValueError as e:
            print(f"Error: {e}", file=_err)
            return

    _write_json(mols, args.output)


def command_run(args) -> None:
    from .run import parse

    prms = {}
    for item in args.param or []:
        k, _, v = item.partition("=")
        if k:
            prms[k] = v

    if getattr(args, "reference", None):
        prms.setdefault("reference", args.reference)
    if getattr(args, "target", None):
        prms.setdefault("target", args.target)
    args.params = prms

    if getattr(args, "build", False):
        from .run.topology_builder import TopologyBuilder
        builder = TopologyBuilder()
        builder.run()
        try:
            mols = builder.build()
        except ValueError as e:
            print(f"Error: {e}", file=_err)
            return
        if "reference" not in prms:
            try:
                ref = input("  Reference coordinate file: ").strip()
            except (EOFError, KeyboardInterrupt):
                print("\nCancelled.", file=_err)
                return
            if ref:
                prms["reference"] = ref
        if "target" not in prms:
            try:
                tgt = input("  Target coordinate file (blank = skip): ").strip()
            except (EOFError, KeyboardInterrupt):
                tgt = ""
            if tgt:
                prms["target"] = tgt
        arg = parse.parser(args, mols=mols)
    else:
        arg = parse.parser(args)

    run_kwargs = {k: v for k, v in arg.prms.items() if k not in ("output", "trgmdt", "refmdt", "overwrite")}
    arg.run(
        arg.mols,
        arg.refxyz,
        arg.trgxyz,
        output=getattr(args, "output", None),
        trgmdt=getattr(arg, "trgmdt", None),
        refmdt=getattr(arg, "refmdt", None),
        overwrite=getattr(args, "overwrite", False),
        **run_kwargs,
    )


def command_demo(args) -> None:
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

    output_spec = getattr(args, "output", None)

    demo_list = [
        cogen.__demo(cli=cli, prec=prec),
        bb.__demo(cli=cli, prec=prec),
        bb_2d.__demo(cli=cli, prec=prec),
        bb_multi.__demo(cli=cli, prec=prec),
        batch.__demo(cli=cli, prec=prec),
        batch_tri.__demo(cli=cli, prec=prec),
        mst.__demo(cli=cli, prec=prec),
    ]

    def resolve_output(demo_idx):
        if output_spec is None:
            return None
        if output_spec == "":
            ext = demo_list[demo_idx].default_ext
            return f"mobbrmsd_demo_{demo_idx + 1:02d}{ext}"
        return output_spec

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
            if 0 <= index < len(demo_list):
                output = resolve_output(index)
                demo_list[index].run_demo(output=output, **args.argv)
            inp = input('  press any key to continue ("q" to exit) >> ')
            print(inp)
            if inp == "":
                continue
            if inp[0] == "q":
                break
    else:
        output = resolve_output(no)
        demo_list[no].run_demo(output=output, **args.argv)


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

    parser_run = sub.add_parser(
        "run",
        help="run mobbrmsd",
        usage="%(prog)s [-i JSON | -b] [-p KEY=VALUE ...] [reference [target]]",
    )
    run_input = parser_run.add_mutually_exclusive_group()
    run_input.add_argument(
        "-i",
        "--inp",
        nargs=1,
        default=None,
        metavar="JSON",
        help="Topology JSON file (molecular_system definition)",
    )
    run_input.add_argument(
        "-b",
        "--build",
        action="store_true",
        default=False,
        help="Launch topology builder interactively before running",
    )
    parser_run.add_argument(
        "reference",
        nargs="?",
        default=None,
        help="Reference coordinate file",
    )
    parser_run.add_argument(
        "target",
        nargs="?",
        default=None,
        help="Target coordinate file (optional)",
    )
    parser_run.add_argument(
        "-o",
        "--output",
        default=None,
        metavar="FILE",
        help="Save rotated target structure to FILE (MDTraj-supported format, e.g. .pdb, .xtc)",
    )
    parser_run.add_argument(
        "-O",
        "--overwrite",
        action="store_true",
        default=False,
        help="Overwrite output file(s) if they already exist (default: warn and skip)",
    )
    parser_run.add_argument(
        "-p",
        "--param",
        action="append",
        default=None,
        metavar="KEY=VALUE",
        help="Extra run parameter (repeatable), e.g. -p maxeval=1000 -p cutoff=2.0",
    )
    parser_run.set_defaults(handler=command_run)

    parser_build = sub.add_parser(
        "build", help="build molecular topology interactively and save to JSON"
    )
    parser_build.add_argument(
        "-o",
        "--output",
        default=None,
        metavar="FILE",
        help="Output JSON file (default: print to stdout)",
    )
    parser_build.add_argument(
        "-f",
        "--file",
        default=None,
        metavar="FILE",
        help="Parse topology from an MDTraj-supported file (PDB, GRO, …) instead of the interactive builder",
    )
    parser_build.add_argument(
        "--topology",
        default=None,
        metavar="TOP",
        help="Topology file for formats that require one (e.g. a PDB alongside a DCD)",
    )
    parser_build.set_defaults(handler=command_build)

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
        "-o", "--output",
        nargs="?",
        const="",
        default=None,
        metavar="FILE",
        help=(
            "Save plot to FILE instead of displaying. "
            "For animations: .gif/.webp saves animated file, .png/.jpg saves numbered frames. "
            "Omit FILE to use default name mobbrmsd_demo_NN.png."
        ),
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
