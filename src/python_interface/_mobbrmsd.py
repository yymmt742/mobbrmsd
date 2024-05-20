import numpy
import networkx


class mobbrmsd_result:
    def __init__(
        self,
        driver,
        header: numpy.ndarray,
        istate: numpy.ndarray,
        rstate: numpy.ndarray,
    ) -> None:
        self.header = header
        self.state = (istate, rstate)
        self.rmsd = driver.rmsd(rstate)
        self.autocorr = driver.autocorr(rstate)
        self.bounds = driver.bounds(rstate)
        self.n_eval = driver.n_eval(rstate)
        self.log_eval_ratio = driver.log_eval_ratio(rstate)
        self.eval_ratio = numpy.exp(self.log_eval_ratio)
        self.is_finished = driver.is_finished(header, istate, rstate)

    def __repr__(self):
        kws = [f"{key}={value!r}" for key, value in self.__dict__.items()]
        return "{}({})".format(type(self).__name__, ", ".join(kws))


class mobbrmsd:
    def __init__(
        self,
        molecules: dict | list[dict] = {"n_apm": 1, "n_mol": 1, "sym": None},
        d: int = 3,
    ) -> None:

        if d == 2:
            from .mobbrmsd_2d import driver

            self.driver = driver
        elif d == 3:
            from .mobbrmsd_3d import driver

            self.driver = driver
        elif d == 1 or d > 3:
            from .mobbrmsd import driver

            self.driver = driver
            self.driver.setup_dimension(d)
        else:
            raise ValueError

        def add_molecule(self, mol: dict) -> None:
            n_apm = int(mol.get("n_apm"))
            n_mol = int(mol.get("n_mol"))
            if n_apm < 1 or n_mol < 1:
                raise ValueError
            sym = mol.get("sym")
            if sym is None:
                self.molecules.append({"n_apm": n_apm, "n_mol": n_mol, "sym": sym})
                self.driver.add_molecule(n_apm, n_mol, 1)
            else:
                try:
                    sym_ = numpy.array(sym).reshape((-1, n_apm))
                except:
                    raise ValueError
                try:
                    self.driver.add_molecule(
                        n_apm,
                        n_mol,
                        sym_.shape[0] + 1,
                        sym_.flatten() + 1,
                    )
                except:
                    raise ValueError
                self.molecules.append(
                    {"n_apm": n_apm, "n_mol": n_mol, "sym": sym_.tolist()}
                )

        self.molecules = []

        if type(molecules) is list:
            for mol in molecules:
                add_molecule(self, mol)
        else:
            add_molecule(self, molecules)

        self.d, self.natom = self.driver.n_atoms()
        self.memsize, self.njob = self.driver.workmemory_lengthes()
        self.n_header, self.n_int, self.n_float = self.driver.state_vector_lengthes()

    def __del__(self):
        self.clear()
        del self.driver

    def __str__(self):
        kws = [
            f"molecules={self.molecules}",
            f"d={self.d}",
        ]
        return "{}({})".format(type(self).__name__, ", ".join(kws))

    def __repr__(self):
        kws = [f"{key}={value!r}" for key, value in self.__dict__.items()]
        return "{}({})".format(type(self).__name__, ", ".join(kws))

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
        Y: None | numpy.ndarray = None,
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
        y: None | numpy.ndarray = None,
        cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        remove_com: bool = True,
        sort_by_g: bool = True,
        rotate_y: bool = False,
        verbose: bool = False,
    ) -> list:

        if hasattr(self, "ww"):
            if self.ww.shape[1] != self.memsize or self.ww.shape[0] != self.njob:
                self.ww = numpy.empty((self.njob, self.memsize)).T
        else:
            self.ww = numpy.empty((self.njob, self.memsize)).T

        x_ = self.varidation_coordinates_2(x)
        if y is None:
            n_target = x_.shape[2]
            hret, iret, rret = self.driver.batch_run_tri(
                n_target,
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
            )

            def res(i, j):
                if i > j:
                    k = int(i * (i - 1) / 2) + j
                    return mobbrmsd_result(self.driver, hret, iret.T[k], rret.T[k])
                elif i < j:
                    k = int(j * (j - 1) / 2) + i
                    return mobbrmsd_result(self.driver, hret, iret.T[k], rret.T[k])
                else:
                    return None

            return [[res(i, j) for j in range(n_target)] for i in range(n_target)]

        else:
            y_ = self.varidation_coordinates_2(y)

            n_reference = x_.shape[2]
            n_target = y_.shape[2]

            hret, iret, rret = self.driver.batch_run(
                n_reference,
                n_target,
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
                [mobbrmsd_result(self.driver, hret, ir, rr) for ir, rr in zip(iry, rry)]
                for iry, rry in zip(
                    iret.transpose([2, 1, 0]), rret.transpose([2, 1, 0])
                )
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
        n_target = x.shape[0]

        if hasattr(self, "ww"):
            if self.ww.shape[1] != self.memsize or self.ww.shape[0] != self.njob:
                self.ww = numpy.empty((self.njob, self.memsize)).T
        else:
            self.ww = numpy.empty((self.njob, self.memsize)).T

        edges, weights, hret, iret, rret = self.driver.min_span_tree(
            n_target,
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

        if x.ndim != 2:
            raise ValueError

        if x.shape[1] != self.d or x.shape[0] != self.natom:
            raise ValueError
        return x.flatten()

    def varidation_coordinates_2(self, x: numpy.ndarray) -> numpy.ndarray:

        if x.ndim == 2:
            x_ = x.transpose().reshape((x.shape[1], x.shape[0], 1))
        elif x.ndim == 3:
            x_ = x.transpose([2, 1, 0])
        else:
            raise ValueError
        if x_.shape[0] != self.d or x_.shape[1] != self.natom:
            raise ValueError
        return x_

    def clear(self) -> None:
        self.driver.clear_molecule()
        self.d, self.natom = 0, 0
        self.memsize, self.njob = 0, 0
        self.n_header, self.n_int, self.n_float = 0, 0, 0
