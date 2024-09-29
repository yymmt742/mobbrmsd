from . import __version__
import time
import numpy
import pprint

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
    i = 0
    for iper in per:
        i += 1
        if i > n_sym:
            break
        sym += [iper]
    return sym


def print_system(n_apm, n_mol, sym):

    print(
        f"    Atoms per molecule :{n_apm:6d}",
    )
    print(f"    Number of molecule :{n_mol:6d}")
    pp = pprint.pformat(tuple([i for i in range(n_apm)]), width=48, compact=True).split(
        "\n"
    )
    print("    Molecular symmetry :     0", pp[0])
    for l in pp[1:]:
        print("                              ", l)
    for i, s in enumerate(sym):
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


class _demo:

    def __init__(self, title="mobbrmsd demo", prec=numpy.float64, cli=False):
        self.prec = numpy.dtype(prec)
        self.cli = cli
        if self.prec == numpy.float64:
            self.title = title + " (double)"
        elif self.prec == numpy.float32:
            self.title = title + " (single)"
        else:
            raise ValueError

    def read_input(self):
        return

    def demo(self, **kwarg):
        return

    def after(self, **kwarg):
        return

    def run_demo(self, **kwarg):
        print(bar2)
        if self.cli:
            prms = {**self.read_input(), **kwarg}
        else:
            prms = kwarg

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

        if ret is not None:
            self.after(**ret)

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
