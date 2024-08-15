from . import __version__
import time
import numpy


class _demo:
    bar1 = "  ------------------------------------------------------------------------------"
    bar2 = "  =============================================================================="
    bar3 = "==========================================="

    def __init__(self, prec=numpy.float64):
        self.title = "mobbrmsd demo"
        self.prec = prec

    def read_input(self):
        return None

    def demo(self):
        return None

    def after(self):
        return None

    def run_demo(self):
        print(bar2)
        prms = self.read_input()

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

        if ret is not None:
            self.after(ret)
            print(bar2, "\n")
