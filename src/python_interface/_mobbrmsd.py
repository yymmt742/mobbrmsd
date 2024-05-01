import numpy
import networkx


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
        self.rmsd = driver.rmsd(rstate)
        self.bounds = driver.bounds(rstate)
        self.n_eval = driver.n_eval(rstate)
        self.log_eval_ratio = driver.log_eval_ratio(rstate)
        self.eval_ratio = numpy.exp(self.log_eval_ratio)
        self.is_finished = driver.is_finished(header, istate, rstate)


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

        self.ndim, self.natom = self.driver.n_atoms()
        self.memsize, self.njob = self.driver.workmemory_lengthes()
        self.n_header, self.n_int, self.n_float = self.driver.state_vector_lengthes()

    def add_molecule(self, n_apm, n_mol=1, swp=None):

        if swp is None:
            self.driver.add_molecule(n_apm, n_mol, 1)
        else:
            swp_ = numpy.array(swp).reshape((-1, n_apm)) + 1
            s = swp_.shape[0] + 1
            self.driver.add_molecule(n_apm, n_mol, s, swp_.flatten())

        self.ndim, self.natom = self.driver.n_atoms()
        self.memsize, self.njob = self.driver.workmemory_lengthes()
        self.n_header, self.n_int, self.n_float = self.driver.state_vector_lengthes()

    def run(
        self,
        x: numpy.ndarray,
        y: numpy.ndarray,
        cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        remove_com: bool = True,
        sort_by_g: bool = True,
        rotate_y: bool = False,
    ) -> mobbrmsd_result:

        x_ = self.varidation_coordinates_1(x)
        y_ = self.varidation_coordinates_1(y)

        if hasattr(self, "w"):
            if self.w.size != self.memsize:
                self.w = numpy.empty(self.memsize)
        else:
            self.w = numpy.empty(self.memsize)

        hret, iret, rret = self.driver.run(
            self.n_header,
            self.n_int,
            self.n_float,
            x_,
            y_,
            self.w,
            cutoff,
            difflim,
            maxeval,
            remove_com,
            sort_by_g,
            rotate_y,
        )

        return mobbrmsd_result(self.driver, hret, iret, rret)

    def restart(
        self,
        ret: mobbrmsd_result,
        cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        Y: numpy.ndarray = None,
    ) -> mobbrmsd_result:

        if not hasattr(self, "w"):
            raise ValueError

        self.driver.restart(
            ret.header,
            ret.state[0],
            ret.state[1],
            self.w,
            cutoff,
            difflim,
            maxeval,
        )
        if Y is not None:

            if Y.ndim != 2:
                raise ValueError

            self.driver.rotate_y(ret.header, ret.state[0], ret.state[1], Y.T)

        return mobbrmsd_result(self.driver, ret.header, ret.state[0], ret.state[1])

    def batch_run(
        self,
        x: numpy.ndarray,
        y: numpy.ndarray,
        cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        remove_com: bool = True,
        sort_by_g: bool = True,
        rotate_y: bool = False,
    ) -> list:

        x_ = self.varidation_coordinates_1(x)
        y_ = self.varidation_coordinates_2(y)

        if hasattr(self, "ww"):
            if self.ww.shape[1] != self.memsize or self.ww.shape[0] != self.njob:
                self.ww = numpy.empty((self.njob, self.memsize)).T
        else:
            self.ww = numpy.empty((self.njob, self.memsize)).T

        hret, iret, rret = self.driver.batch_run(
            self.n_header,
            self.n_int,
            self.n_float,
            x_,
            y_,
            self.ww,
            cutoff,
            difflim,
            maxeval,
            remove_com,
            sort_by_g,
            rotate_y,
        )

        return [
            mobbrmsd_result(self.driver, hret, ir, rr) for ir, rr in zip(iret.T, rret.T)
        ]

    def min_span_tree(
        self,
        x: numpy.ndarray,
        cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        remove_com: bool = True,
        sort_by_g: bool = True,
        verbose: bool = False,
    ) -> tuple:

        x_ = self.varidation_coordinates_2(x)

        if hasattr(self, "ww"):
            if self.ww.shape[1] != self.memsize or self.ww.shape[0] != self.njob:
                self.ww = numpy.empty((self.njob, self.memsize)).T
        else:
            self.ww = numpy.empty((self.njob, self.memsize)).T

        edges, weights, hret, iret, rret = self.driver.min_span_tree(
            self.n_header,
            self.n_int,
            self.n_float,
            x_,
            self.ww,
            cutoff,
            difflim,
            maxeval,
            remove_com,
            sort_by_g,
            verbose,
        )

        g = networkx.Graph()
        for e, w in zip(edges.T, weights):
            g.add_edge(e[0] - 1, e[1] - 1, weight=w)

        return (
            g,
            [
                [
                    mobbrmsd_result(self.driver, hret, irij, rrij)
                    for irij, rrij in zip(iri, rri)
                ]
                for iri, rri in zip(iret.T, rret.T)
            ],
        )

    def varidation_coordinates_1(self, x: numpy.ndarray) -> numpy.ndarray:

        if x.ndim == 2:
            x_ = x.transpose()
        else:
            raise ValueError

        if x_.shape[0] != self.ndim or x_.shape[1] != self.natom:
            raise ValueError
        return x_

    def varidation_coordinates_2(self, x: numpy.ndarray) -> numpy.ndarray:

        if x.ndim == 2:
            x_ = x.transpose().reshape((x.shape[1], x.shape[0], 1))
        elif x.ndim == 3:
            x_ = x.transpose([2, 1, 0])
        else:
            raise ValueError
        if x_.shape[0] != self.ndim or x_.shape[1] != self.natom:
            raise ValueError
        return x_

    def clear(self) -> None:
        self.driver.clear_molecule()
        self.ndim, self.natom = self.driver.n_atoms()
        self.memsize, self.njob = self.driver.workmemory_lengthes()
        self.n_header, self.n_int, self.n_float = self.driver.state_vector_lengthes()
