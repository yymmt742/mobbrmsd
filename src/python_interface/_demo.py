from . import __version__
import time
import numpy

bar1 = (
    "  ------------------------------------------------------------------------------"
)
bar2 = (
    "  =============================================================================="
)
bar3 = "==========================================="


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


def yes_or_no(msg):
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


class _demo:

    def __init__(self, title="mobbrmsd demo", prec=numpy.float64):
        self.title = title
        self.prec = prec

    def read_input(self):
        return

    def demo(self, **kwarg):
        return

    def after(self, **kwarg):
        return

    def run_demo(self):
        print(bar2)
        prms = self.read_input()

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
            print(bar2, "\n")
