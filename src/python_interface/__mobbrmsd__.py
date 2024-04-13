import numpy


class mobbrmsd_result:
    def __init__(
        self,
        driver,
        header: numpy.ndarray,
        istate: numpy.ndarray,
        rstate: numpy.ndarray,
    ):
        self.header = header
        self.state = (istate, rstate)
        self.rmsd = driver.rmsd(istate, rstate)
        self.bounds = driver.bounds(istate, rstate)
        self.n_eval = driver.n_eval(istate, rstate)
        self.log_eval_ratio = driver.log_eval_ratio(istate, rstate)
        self.eval_ratio = numpy.exp(self.log_eval_ratio)


class mobbrmsd:
    def __init__(self, d: int = 3):

        if d == 2:
            from .mobbrmsd_2d import driver

            self.driver = driver
        elif d == 3:
            from .mobbrmsd_3d import driver

            self.driver = driver
        else:
            from .mobbrmsd import driver

            self.driver = driver
            self.driver.setup_dimension(d)

    def add_molecule(self, m=1, n=1, swp=None):

        if swp is None:
            self.driver.add_molecule(m, n, 1)
        else:
            swp_ = numpy.array(swp).reshape((-1, m))
            s = swp_.shape[0] + 1
            self.driver.add_molecule(m, n, s, swp_)

    def run(
        self,
        x: numpy.ndarray,
        y: numpy.ndarray,
        cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        rotate_y: bool = False,
    ) -> mobbrmsd_result | list:

        ndim = self.driver.n_dims()
        natom = self.driver.n_atoms()

        if x.ndim == 2:
            x_ = x.transpose()
        else:
            raise ValueError

        if y.ndim == 2:
            y_ = y.transpose().reshape((-1, y.shape[1], y.shape[0]))
        elif y.ndim == 3:
            y_ = y.transpose([2, 1, 0])
        else:
            raise ValueError

        if (
            x_.shape[0] != ndim
            or x_.shape[1] != natom
            or y_.shape[0] != ndim
            or y_.shape[1] != natom
        ):
            raise ValueError

        ntarget = y_.shape[2]
        n_header, n_int, n_float = self.driver.state_vector_lengthes()
        hret, iret, rret = self.driver.batch_run(
            n_header, n_int, n_float, x_, y_, cutoff, difflim, maxeval, rotate_y
        )

        if ntarget > 1:
            return [
                mobbrmsd_result(self.driver, hret, ir, rr)
                for ir, rr in zip(iret.T, rret.T)
            ]
        else:
            return mobbrmsd_result(self.driver, hret, iret, rret)
