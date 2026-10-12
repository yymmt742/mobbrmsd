import time
import numpy
import pprint
import sys
from pathlib import Path

bar1 = (
    "  ------------------------------------------------------------------------------"
)
bar2 = (
    "  =============================================================================="
)
bar3 = "==========================================="


def generate_sym_indices(n_apm: int, n_sym: int = 1):

    import itertools

    per = itertools.permutations(range(n_apm))
    next(per)
    sym = []
    i = 1
    for iper in per:
        i += 1
        if i > n_sym:
            break
        sym += [iper]
    return sym


def print_system(molsys, title):
    print(bar1, file=sys.stderr)
    print(" " * ((len(bar1) - len(title)) // 2) + title, file=sys.stderr)
    print(bar1, file=sys.stderr)
    print("      --System settings--", file=sys.stderr)
    for mol in molsys.mols:
        print(
            f"    Atoms per molecule :{mol.n_apm:6d}",
            file=sys.stderr,
        )
        print(f"    Number of molecule :{mol.n_mol:6d}", file=sys.stderr)
        pp = pprint.pformat(
            tuple([i for i in range(mol.n_apm)]), width=48, compact=True
        ).split("\n")
        print("    Molecular symmetry :     0", pp[0], file=sys.stderr)
        for l in pp[1:]:
            print("                              ", l, file=sys.stderr)
        for i, s in enumerate(mol.sym):
            pp = pprint.pformat(s, width=48, compact=True).split("\n")
            print(f"                        {i+1:6d}", pp[0], file=sys.stderr)
            for l in pp[1:]:
                print("                              ", l, file=sys.stderr)
        print(file=sys.stderr)


def readinp(msg, default, check=None):
    sys.stderr.write(f"\n    {msg} (default : {default})   >> ")
    sys.stderr.flush()
    inp = sys.stdin.readline().rstrip("\n")

    def isfloat(string):
        try:
            float(string)
            return True
        except ValueError:
            return False

    while True:
        if inp == "":
            sys.stderr.write(
                f"\033[1A    {msg} (default : {default})   >> {default}\n"
            )
            sys.stderr.flush()
            return default
        if inp[0] == "q" or inp[0] == "Q":
            print(bar2, file=sys.stderr)
            exit()
        elif inp.isdigit():
            ret = int(inp)
        elif isfloat(inp):
            ret = float(inp)
        else:
            ret = inp

        if check is not None:
            if check(ret):
                return ret
        else:
            return ret
        sys.stderr.write(
            f"    | [{ret}] is invalid. \033[1F\033[0K    {msg} (default : {default})   >> "
        )
        sys.stderr.flush()
        inp = sys.stdin.readline().rstrip("\n")


def print_ret(
    ret, post="", end="\n", to_console: bool = False, header=False, footer=False
):
    sep1 = "  ------------------------------------------------------------------------"
    if footer or header:
        print(sep1)
        if footer:
            print("      -- Final results --")
        print(
            "    N_eval  Eval_ratio  Lowerbound  Upperbound RMSD(LB) RMSD(UB) RMSD_Gap"
        )
    ev = ret.n_eval()
    er = ret.eval_ratio()
    sb = ret.bounds()
    rb = ret.bounds_as_rmsd()
    gap = rb[1] - rb[0]
    if sys.stdout.isatty():
        print(
            f"\r{ev:10d}{er:12.8f}{sb[0]:12.4f}{sb[1]:12.4f}{rb[0]:9.3f}{rb[1]:9.3f}{gap:9.3f}",
            post,
            end=end,
        )
    else:
        if to_console:
            return
        print(
            f"{ev:10d}{er:12.8f}{sb[0]:12.4f}{sb[1]:12.4f}{rb[0]:9.3f}{rb[1]:9.3f}{gap:9.3f}"
        )
    if footer:
        print(sep1, "\n")


def _save_figure(output, fig, rot_func=None, rot_frames=None):
    """Save a matplotlib figure or animation to *output*.

    When *rot_func* and *rot_frames* are given the output is animated:
      - .gif / .webp  → single animated file (Pillow writer)
      - .png / .jpg   → numbered frame files  (foo.001.png, foo.002.png, …)
      - other         → static snapshot of the first frame

    When no animation function is given, ``fig.savefig(output)`` is used.
    """
    p = Path(output)
    suffix = p.suffix.lower()

    if rot_func is None or rot_frames is None:
        fig.savefig(output)
        return

    rot_frames = list(rot_frames)

    if suffix in (".gif", ".webp"):
        import matplotlib.animation as _anim
        ani = _anim.FuncAnimation(fig, rot_func, frames=rot_frames, interval=50)
        ani.save(output, writer="pillow")
        print(f"  Animation saved to {output}", file=sys.stderr)
    elif suffix in (".png", ".jpg", ".jpeg"):
        width = max(3, len(str(len(rot_frames))))
        for k, frame_val in enumerate(rot_frames):
            rot_func(frame_val)
            frame_path = str(p.parent / f"{p.stem}.{k + 1:0{width}d}{p.suffix}")
            fig.savefig(frame_path)
        print(
            f"  {len(rot_frames)} frames saved as "
            f"{p.stem}.001{p.suffix} … {p.stem}.{len(rot_frames):0{width}d}{p.suffix}",
            file=sys.stderr,
        )
    else:
        rot_func(rot_frames[0])
        fig.savefig(output)
        print(f"  Figure saved to {output}", file=sys.stderr)


class _demo:

    default_ext = ".png"

    def __init__(self, title="mobbrmsd demo", prec=numpy.float64, cli=False):
        self.prec = numpy.dtype(prec)
        self.cli = cli
        if self.prec == numpy.float64:
            self.title = title + " (double precision)"
        elif self.prec == numpy.float32:
            self.title = title + " (single precision)"
        else:
            raise ValueError

    def read_input(self):
        return

    def demo(self, **kwargs):
        return

    def after(self, **kwargs):
        return

    def run_demo(self, output=None, **kwargs):
        print(bar2, file=sys.stderr)
        if self.cli:
            prms = {**kwargs, **self.read_input()}
        else:
            prms = kwargs

        start_wallclock_time = time.time()
        start_cpu_time = time.process_time()

        ret = self.demo(**prms)

        end_wallclock_time = time.time()
        end_cpu_time = time.process_time()
        elapsed_cpu_time = end_wallclock_time - start_wallclock_time
        elapsed_wallclock_time = end_wallclock_time - start_wallclock_time

        print("    -- Elapsed times --", file=sys.stderr)
        print(
            f"       cpu time : {elapsed_cpu_time:16.3f} sec",
            f" wallclock time : {elapsed_wallclock_time:16.3f} sec",
            file=sys.stderr,
        )
        print(bar2, "\n", file=sys.stderr)

        if ret is None:
            self.after(output=output, **kwargs)
        else:
            self.after(output=output, **{**kwargs, **ret})

    def yes_or_no(self, msg):
        if not self.cli:
            return True
        sys.stderr.write(f"\n\033[1A\033[0K    {msg} [Y/n] > ")
        sys.stderr.flush()
        while True:
            inp = sys.stdin.readline().rstrip("\n")
            if inp == "":
                sys.stderr.write(f"\033[1A\033[0K    {msg} [Y/n] > ")
                sys.stderr.flush()
                continue
            if inp[0] == "q" or inp[0] == "Q":
                print(bar2, file=sys.stderr)
                exit()
            if inp[0] == "y" or inp[0] == "Y" or inp[0] == "n" or inp[0] == "N":
                return inp[0] == "y" or inp[0] == "Y"
            sys.stderr.write(f"\033[1A\033[0K    {msg} [Y/n] > ")
            sys.stderr.flush()
