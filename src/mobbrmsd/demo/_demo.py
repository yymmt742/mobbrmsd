from . import __version__
import time
import numpy
import pprint
import sys
import pick

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


def print_system(molecules, title):

    print(bar1)
    print(" " * ((len(bar1) - len(title)) // 2) + title)
    print(bar1)
    print("      --System settings--")
    for mol in molecules:

        print(
            f"    Atoms per molecule :{mol.n_apm:6d}",
        )
        print(f"    Number of molecule :{mol.n_mol:6d}")
        pp = pprint.pformat(
            tuple([i for i in range(mol.n_apm)]), width=48, compact=True
        ).split("\n")
        print("    Molecular symmetry :     0", pp[0])
        for l in pp[1:]:
            print("                              ", l)
        for i, s in enumerate(mol.sym):
            pp = pprint.pformat(s, width=48, compact=True).split("\n")
            print(f"                        {i+1:6d}", pp[0])
            for l in pp[1:]:
                print("                              ", l)
        print()


def readinp(msg, default, check=None):
    print()
    inp = input(f"\033[1A    {msg} (default : {default})   >> ")

    def isfloat(string):
        try:
            float(string)
            return True
        except ValueError:
            return False

    while True:
        if inp == "":
            print(f"\033[1A    {msg} (default : {default})   >> {default}")
            return default
        if inp[0] == "q" or inp[0] == "Q":
            print(bar2)
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
        inp = input(
            f"    | [{ret}] is invalid. \033[1F\033[0K    {msg} (default : {default})   >> "
        )


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


class _demo:

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

    def run_demo(self, **kwargs):
        print(bar2)
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

        print("    -- Elapsed times --")
        print(
            f"       cpu time : {elapsed_cpu_time:16.3f} sec",
            f" wallclock time : {elapsed_wallclock_time:16.3f} sec",
        )
        print(bar2, "\n")

        if ret is None:
            self.after(**kwargs)
        else:
            self.after(**{**kwargs, **ret})

    def yes_or_no(self, msg):
        if not self.cli:
            return True
        print()
        while True:
            inp = input(f"\033[1A\033[0K    {msg} [Y/n] > ")
            if inp == "":
                continue
            if inp[0] == "q" or inp[0] == "Q":
                print(bar2)
                exit()
            if inp[0] == "y" or inp[0] == "Y" or inp[0] == "n" or inp[0] == "N":
                return inp[0] == "y" or inp[0] == "Y"
