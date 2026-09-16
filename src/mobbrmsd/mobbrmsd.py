# -*- coding: utf-8 -*-
"""mobbrmsd driver routine"""

import numpy
import numpy.typing as npt
import networkx
from typing import Union, Optional, List
from tqdm import trange
from .dataclass import molecules, molecular_system


class mobbrmsd_result:
    """Manages calculation results and memory for recalculations.

    Args:
        driver (): fortran driver
        d (int): spatial dimension
        header (npt.NDArray): header array
        istate (Optional[npt.NDArray]): state array (integer part)
        rstate (Optional[npt.NDArray]): state array (real part)
        rot (Optional[npt.NDArray]): rotation matrix, given by flatten [d*d]
        w (Optional[npt.NDArray]): working memory
    """

    def __init__(
        self,
        driver,
        d: int,
        header: npt.NDArray,
        istate: Optional[npt.NDArray] = None,
        rstate: Optional[npt.NDArray] = None,
        rot: Optional[npt.NDArray] = None,
        w: Optional[npt.NDArray] = None,
    ) -> None:

        self.header = header.copy()
        self.d = d

        if (istate is None) or (rstate is None):
            self.istate = numpy.array(["""DEFAULT_ISTATE"""])
            self.rstate = numpy.array(["""DEFAULT_RSTATE"""])
        else:
            self.istate = istate.copy()
            self.rstate = rstate.copy()

        if rot is None:
            self.rot = numpy.zeros([self.d * self.d])
        else:
            self.rot = rot.copy()

        if w is not None:
            self.w = w.copy()

    def autocorr(self) -> float:
        r"""Autocorrelation term, given by
           $G = \text{tr} [XX^{\top}] + \text{tr}[YY^\top]$

        Returns:
            ret (float): Autocorrelation, $G$.
        """
        return float(self.rstate["""INDEX_TO_AUTOCORR"""])

    def lowerbound(self) -> float:
        r"""Lowerbound of objective function,
            $L(X,Y) = \max_{R,S} \text{tr}[RYSX^\top]$

        Returns:
            ret (float): Lowerbound of objective function, $L$.
        """
        return float(self.rstate["""INDEX_TO_LOWERBOUND"""])

    def upperbound(self) -> float:
        r"""Upperbound of objective function,
            $L(X,Y) = \max_{R,S} \text{tr}[RYSX^\top]$

        Returns:
            ret (float): upperbound of objective function, $L$.
        """
        return float(self.rstate["""INDEX_TO_UPPERBOUND"""])

    def lowerbound_as_rmsd(self) -> float:
        r"""Lowerbound of objective function on RMSD scale (multiplied by 1/n and rooted).

        Returns:
            ret (float): Scaled lowerbound of objective function, $\sqrt{\frac1n \underline L}$.
        """
        rn = self.rstate["""RECIPROCAL_OF_N"""]
        return numpy.sqrt(
            numpy.max([0.0, rn * (2 * self.lowerbound() + self.autocorr())])
        )

    def upperbound_as_rmsd(self) -> float:
        r"""Upperbound of objective function on RMSD scale (multiplied by 1/n and rooted).

        Returns:
            ret (float): Scaled lowerbound of objective function, $\sqrt{\frac1n \bar L}$.
        """
        rn = self.rstate["""RECIPROCAL_OF_N"""]
        return numpy.sqrt(
            numpy.max([0.0, rn * (2 * self.upperbound() + self.autocorr())])
        )

    def sd(self) -> float:
        r"""Squared displacement.
            If the calculation terminates early,
            this is a preliminary value and may be larger than the true solution.

        Returns:
            ret (float): Squared displacement.
        """
        ret = float(2 * self.upperbound() + self.autocorr())
        if ret < 0.0:
            return 0.0
        else:
            return ret

    def msd(self) -> float:
        r"""Mean squared displacement.
            If the calculation terminates early,
            this is a preliminary value and may be larger than the true solution.

        Returns:
            ret (float): Mean squared displacement.
        """

        rn = self.rstate["""RECIPROCAL_OF_N"""]
        return float(rn * self.sd())

    def rmsd(self) -> float:
        r"""Root mean squared displacement.
            If the calculation terminates early,
            this is a preliminary value and may be larger than the true solution.

        Returns:
            ret (float): Root mean squared displacement.
        """

        return float(numpy.sqrt(self.msd()))

    def bounds(self) -> npt.NDArray:
        r"""Interval containing an exaxt solution.

        Returns:
            ret (npt.NDArray): Lowerbound and upperbound, $\left[\underline L, \bar L\right]$.
        """
        return numpy.array([self.lowerbound(), self.upperbound()])

    def bounds_as_rmsd(self) -> npt.NDArray:
        r"""Interval containing an exaxt solution on RMSD scale (multiplied by 1/n and rooted).

        Returns:
            ret (npt.NDArray): Lowerbound and upperbound, $\left[\sqrt{\frac1n \underline L}, \sqrt{\frac1n \bar L}\right]$.
        """
        rn = self.rstate["""RECIPROCAL_OF_N"""]
        return rn * (2 * self.bounds() + self.autocorr())

    def n_eval(self) -> int:
        r"""Number of Nodes Evaluated.

        Returns:
            ret(int): Number of evaluations.
        """
        return int(self.rstate["""INDEX_TO_N_EVAL"""])

    def log_eval_ratio(self) -> float:
        r"""
        Returns:
            ret(float): Percentage of computations relative to exhaustive search (logarithm).
        """
        return self.rstate["""INDEX_TO_LOG_RATIO"""]

    def eval_ratio(self) -> float:
        r"""
        Returns:
            ret(float): Percentage of computations relative to exhaustive search.
        """
        return numpy.exp(self.log_eval_ratio())

    def is_finished(self) -> bool:
        r"""
        Returns:
            ret(bool): True if calculation is complete.
        """
        return self.istate[-1] == ("""IS_FINISHED_FLAG""")

    def restart(
        self,
        cutoff: float = float("inf"),
        ub_cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        difflim_absolute: bool = False,
        get_rotation: bool = False,
    ) -> None:
        r"""Resume a calculation that was interrupted.
            Overwrite the instance with the calculation results.

        Args:
            cutoff (float): Terminate the calculation as soon as the lowerbound reaches `cutoff` or greater (RMSD-based).
            ub_cutoff (float): Terminate the calculation as soon as the upperbound reaches `ub_cutoff` or greater (RMSD-based).
            difflim (float): Terminate the calculation as soon as the difference between the upper and lowerbounds becomes `difflim` or less.
            maxeval (int): Terminate the calculation when the number of evaluated nodes exceeds `maxeval`.
                           At least one pruning cycle will be excuted.
                           When maxeval < 0, the calculation continues until an exact solution is found.
            difflim_absolute (bool): Use `difflim` in terms of RMSD (multiply the evaluation function by 1/n and take the square root).
                                     Scheduled for Discontinuation
            get_rotation (bool): Calculate the rotation matrix.
        """
        if not hasattr(self, "w"):
            return

        driver = _select_driver(self.d, dtype=self.w.dtype)
        ropts = numpy.array([cutoff, ub_cutoff, difflim], dtype=self.w.dtype)
        iopts = numpy.array([maxeval], dtype=self.header.dtype)
        driver.restart(
            self.header,
            self.istate,
            self.rstate,
            self.w,
            self.rot,
            ropts,
            iopts,
            difflim_absolute,
            get_rotation,
        )

        del driver

    def superpose(
        self,
        y: npt.NDArray,
    ) -> npt.NDArray:
        r"""Superpose $Y$.
            Returns swap and rotation target coordinate.

        Args(npt.NDArray):
           y: target coordinates. shape[n, d]

        Returns(npt.NDArray):
           Superposed $Y$, shape[n, d]
        """
        y_ = y.flatten().copy()
        driver = _select_driver(self.d, dtype=y_.dtype)
        driver.rotate_y(self.header, self.istate, self.rstate, self.rot, y_)
        del driver
        return y_.reshape(y.shape)

    def permutation_indices(
        self,
    ) -> npt.NDArray:
        """Permutation_indices.

        Returns:
           Superposed indices, shape[n, d]
        """
        driver = _select_driver(self.d)
        ret = numpy.empty(
            round(1 / self.rstate["""RECIPROCAL_OF_N"""]), dtype=self.header.dtype
        )
        driver.permutation_indices(self.header, self.istate, self.rstate, ret)
        del driver
        return ret

    def __repr__(self):
        kws = [f"{key}={value!r}" for key, value in self.__dict__.items()]
        return "{}({})".format(type(self).__name__, ", ".join(kws))

    def __str__(self):
        ev = self.n_eval()
        er = self.eval_ratio()
        sb = self.bounds()
        rb = self.bounds_as_rmsd()

        return (
            f"{ev:12d} {er:12.8f}{sb[0]:16.6f}{sb[1]:16.6f}{rb[0]:12.6f}{rb[1]:12.6f}"
        )


def _select_driver(d: int, dtype=None):
    dt = numpy.float64 if dtype is None else numpy.dtype(dtype)
    errmsg = lambda d, dt: f"Dimension {d} with {dt} is not supported."
    error = not ((dt == numpy.float64) or (dt == numpy.float32)) or (d < 1)
    if d == 2:
        if dt == numpy.float64:
            try:
                from .mobbrmsd_2ddp import driver
            except ModuleNotFoundError:
                error = True
        elif dt == numpy.float32:
            try:
                from .mobbrmsd_2dsp import driver
            except ModuleNotFoundError:
                error = True
    elif d == 3:
        if dt == numpy.float64:
            try:
                from .mobbrmsd_3ddp import driver
            except ModuleNotFoundError:
                error = True
        elif dt == numpy.float32:
            try:
                from .mobbrmsd_3dsp import driver
            except ModuleNotFoundError:
                error = True
    elif d == 1 or d > 3:
        if dt == numpy.float64:
            try:
                from .mobbrmsd_xddp import driver
            except ModuleNotFoundError:
                error = True
        elif dt == numpy.float32:
            try:
                from .mobbrmsd_xdsp import driver
            except ModuleNotFoundError:
                error = True
        driver.setup_dimension_(d)
    if error:
        raise ValueError(errmsg(d, dt))
    return driver


def _select_dtype(x: npt.NDArray, y: npt.NDArray):
    xdt = x.dtype
    ydt = y.dtype
    if x.dtype == numpy.float64:
        if y.dtype == numpy.float64:
            return numpy.float64
        elif y.dtype == numpy.float32:
            return numpy.float64
        else:
            ValueError
    elif x.dtype == numpy.float32:
        if y.dtype == numpy.float64:
            return numpy.float64
        elif y.dtype == numpy.float32:
            return numpy.float32
        else:
            ValueError
    else:
        ValueError


class mobbrmsd:
    """mobbrmsd driver class.

    Args:
        mols: molecules/molecular_system specifier.
        d (int): Spatial dimension. default=3.
    """

    def __init__(
        self,
        mols: Union[molecules, molecular_system],
        # mols: Union[molecules, molecular_system] = molecules(1, 1),
        d: int = 3,
    ) -> None:

        self.mols = mols

        def add_molecule(m: molecules):
            ret = [m.n_apm, m.n_mol]
            if m.sym is None:
                ret += [1]
            else:
                sym = numpy.array(m.sym, dtype=numpy.int32) + 1
                if sym.ndim == 1:
                    if sym.shape[0] != 0 and sym.shape[0] != m.n_apm:
                        raise ValueError
                elif sym.ndim == 2:
                    if sym.shape[1] != m.n_apm:
                        raise ValueError
                ret += [sym.shape[0] + 1] + sym.flatten().tolist()
            return ret

        ms = []
        if type(mols) is molecular_system:
            for m in mols.mols:
                ms += add_molecule(m)
        elif type(mols) is molecules:
            ms += add_molecule(mols)
        else:
            raise ValueError

        driver = _select_driver(d, dtype=None)
        att = driver.decode_attributes(ms)
        self.d = att[0]
        self.natom = att[1]
        self.n_header = att[2]
        self.n_int = att[3]
        self.n_float = att[4]
        self.n_rot = att[5]
        self.memsize = att[6]
        self.njob = att[7]
        self.header = driver.decode_header(ms, self.n_header)
        del driver

    def rmsd(
        self,
        x: npt.NDArray,
        y: npt.NDArray,
    ) -> float:
        r"""
        Calculate mobbRMSD of a structural pair. (Simplified interface)

        Args:
           x (npt.NDArray): Reference coordinates $X\in\mathbb R^{d\times n}$, shape [n,d].
           y (npt.NDArray): Target coordinates $Y\in\mathbb R^{d\times n}$, shape [n,d].

        Returns:
           ret (float): rmsd value.
        """
        dt = _select_dtype(x, y)
        driver = _select_driver(self.d, dtype=dt)
        w = numpy.empty(self.memsize, dtype=dt)

        _, rret, _ = driver.run(
            self.n_int,  # n_int
            self.n_float,  # n_float
            self.n_rot,  # n_rot
            self.header,  # header
            self.to_rank2_coordinates(x, dtype=dt),  # X
            self.to_rank2_coordinates(y, dtype=dt),  # Y
            w,  # W
            numpy.array([float("inf"), float("inf"), 0.0], dtype=dt),  # ropts
            numpy.array([-1], dtype=numpy.int32),  # iopts
            True,  # remove_com
            True,  # sort_by_g
            False,  # difflim_absolute
            False,  # rotate_y
            False,  # get_rotation
        )  # returns (int_states, float_states, rotation)

        ret = rret["""INDEX_TO_AUTOCORR"""] + 2 * rret["""INDEX_TO_UPPERBOUND"""]
        if ret < 0.0:
            ret = 0.0
        else:
            ret = numpy.sqrt(rret["""RECIPROCAL_OF_N"""] * ret)
        del driver, w
        return ret

    def run(
        self,
        x: npt.NDArray,
        y: npt.NDArray,
        cutoff: float = float("inf"),
        ub_cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        remove_com: bool = True,
        sort_by_g: bool = True,
        difflim_absolute: bool = False,
        rotate_y: bool = False,
        get_rotation: bool = False,
        *args,
        **kwargs,
    ) -> mobbrmsd_result:
        r"""Detailed interface.
            In addition to the calculated RMSD of a structural pair,
            intermediate calculation results are returned.
            If the calculation is terminated prematurely,
            the results retain the data needed to restart the calculation.

        Args:
            x (npt.NDArray): Reference coordinates $X\in\mathbb R^{d\times n}$, shape [n,d].
            y (npt.NDArray): Target coordinates $Y\in\mathbb R^{d\times n}$, shape [n,d].
                             if rotate_y=True, the best structure discovered is assigned.
            cutoff (float): Terminate the calculation as soon as the lowerbound reaches `cutoff` or greater (RMSD-based).
            ub_cutoff (float): Terminate the calculation as soon as the upperbound reaches `ub_cutoff` or greater (RMSD-based).
            difflim (float): Terminate the calculation as soon as the difference between the upper and lowerbounds becomes `difflim` or less.
            maxeval (int): Terminate the calculation when the number of evaluated nodes exceeds `maxeval`.
                           At least one pruning cycle will be excuted.
                           When maxeval < 0, the calculation continues until an exact solution is found.
            remove_com (bool): Remove the center of mass from the coodinates.
            sort_by_g (bool): Sort the reference structures in descending order of self-dispersion and perform the calculation.
                              This affects the computation time required to find an exact solution.
            difflim_absolute (bool): Use `difflim` in terms of RMSD (multiply the evaluation function by 1/n and take the square root).
                                     Scheduled for Discontinuation
            rotate_y (bool): Perform permutations and rotations on the target structure $Y$.
            get_rotation (bool): Calculate the rotation matrix.

        Returns:
          ret (mobbrmsd_result): Calculation results.
                                 Data for restart the calculation is included only in the event of early termination.
        """

        dt = _select_dtype(x, y)
        x_ = self.to_rank2_coordinates(x, dtype=dt)
        y_ = self.to_rank2_coordinates(y, dtype=dt)
        driver = _select_driver(self.d, dtype=dt)
        w = numpy.empty(self.memsize, dtype=dt)
        ropts = numpy.array([cutoff, ub_cutoff, difflim], dtype=dt)
        iopts = numpy.array([maxeval], dtype=numpy.int32)

        iret, rret, rot = driver.run(
            self.n_int,
            self.n_float,
            self.n_rot,
            self.header,
            x_,
            y_,
            w,
            ropts,
            iopts,
            remove_com,
            sort_by_g,
            difflim_absolute,
            rotate_y,
            get_rotation,
        )

        ret = mobbrmsd_result(driver, self.d, self.header, iret, rret, rot, w=w)

        del driver, w, ropts, iopts
        if rotate_y:
            if not numpy.may_share_memory(y, y_):
                y[...] = y_.reshape(y.shape)

        return ret

    def batch_run(
        self,
        x: npt.NDArray,
        y: Optional[npt.NDArray] = None,
        cutoff: float = float("inf"),
        ub_cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        remove_com: bool = True,
        sort_by_g: bool = True,
        difflim_absolute: bool = False,
        verbose: bool = True,
        n_chunk: int = 0,
        *args,
        **kwargs,
    ) -> npt.NDArray:
        r"""
           Batch RMSD runnner for multiple structures.
           Calculates the mobbRMSD matrix for a coordinate system.
           When only structure $X$ is provided,
           it returns the RMSD symmetric matrix $D(X,X)$;
           when both $X$ and $Y$ are provided,
           it returns the RMSD matrix $D(X,Y).

        Args:
            x (npt.NDArray): Reference coordinates $X\in\mathbb R^{d\times n}$, shape [n,d] or [m_X, n, d].
            y (npt.NDArray): Target coordinates $Y\in\mathbb R^{d\times n}$, shape [n,d] or [m_Y, n, d].
            cutoff (float): Terminate the calculation as soon as the lowerbound reaches `cutoff` or greater (RMSD-based).
            ub_cutoff (float): Terminate the calculation as soon as the upperbound reaches `ub_cutoff` or greater (RMSD-based).
            difflim (float): Terminate the calculation as soon as the difference between the upper and lowerbounds becomes `difflim` or less.
            maxeval (int): Terminate the calculation when the number of evaluated nodes exceeds `maxeval`.
                           At least one pruning cycle will be excuted.
                           When maxeval < 0, the calculation continues until an exact solution is found.
            remove_com (bool): Remove the center of mass from the coodinates.
            sort_by_g (bool): Sort the reference structures in descending order of self-dispersion and perform the calculation.
                              This affects the computation time required to find an exact solution.
            difflim_absolute (bool): Use `difflim` in terms of RMSD (multiply the evaluation function by 1/n and take the square root).
                                     Scheduled for Discontinuation
            get_rotation (bool): Calculate the rotation matrix.
            n_chunk (int): The maximum batch size for calculations performed in a single batch.
                           If set to <1, calculations are performed all at once.


        Returns:
            ret (npt.NDArray): A mobbRMSD matrix, shape[$n_X$, $n_X$] if $Y$ is None, else shape[$n_X$, $n_Y$].
        """

        n_eval = 0

        if y is None:
            dt = x.dtype
            x_ = self.to_rank3_coordinates(x)
            driver = _select_driver(self.d, dtype=dt)
            n_target = 1 if x.ndim == 2 else x.shape[0]
            n_tri = (n_target * (n_target - 1)) // 2
            n_chunk_ = n_tri if n_chunk < 1 else self.njob * n_chunk
            ww = numpy.empty((self.njob * self.memsize), dtype=dt)

            ropts = numpy.array([cutoff, ub_cutoff, difflim], dtype=dt)
            iopts = numpy.array([maxeval], dtype=numpy.int32)

            if n_tri == n_chunk_ or not verbose:
                r_tri, log_n_eval = driver.batch_run_tri(
                    n_target,
                    n_tri,
                    1,
                    self.header,
                    x_,
                    ww,
                    ropts,
                    iopts,
                    remove_com,
                    sort_by_g,
                    difflim_absolute,
                )
                n_eval += int(numpy.exp(log_n_eval))
            else:
                n_lower = 1
                nrep = (n_tri + n_chunk_ - 1) // n_chunk_
                r_tri = numpy.empty(n_tri, dtype=dt)
                for i in trange(nrep, *args, **kwargs):
                    l = n_lower - 1
                    u = min([l + n_chunk_, n_tri])
                    r_tri[l:u], log_n_eval = driver.batch_run_tri(
                        n_target,
                        min(n_chunk_, n_tri - n_lower + 1),
                        n_lower,
                        self.header,
                        x_,
                        ww,
                        ropts,
                        iopts,
                        remove_com,
                        sort_by_g,
                        difflim_absolute,
                    )
                    n_lower += n_chunk_
                    n_eval += int(numpy.exp(log_n_eval))
            ret = numpy.zeros([n_target, n_target], dtype=dt)
            k = 0
            for j in range(n_target):
                for i in range(j):
                    ret[i, j] = r_tri[k]
                    k += 1
            ret += ret.T
        else:
            dt = _select_dtype(x, y)
            x_ = self.to_rank3_coordinates(x, dtype=dt)
            y_ = self.to_rank3_coordinates(y, dtype=dt)
            driver = _select_driver(self.d, dtype=dt)
            n_reference = 1 if x.ndim == 2 else x.shape[0]
            n_target = 1 if y.ndim == 2 else y.shape[0]
            n_tri = n_reference * n_target
            n_chunk_ = n_tri if n_chunk < 1 else self.njob * n_chunk
            ww = numpy.empty((self.njob * self.memsize), dtype=dt)

            ropts = numpy.array([cutoff, ub_cutoff, difflim], dtype=dt)
            iopts = numpy.array([maxeval], dtype=numpy.int32)

            if n_tri == n_chunk_ or not verbose:
                ret, log_n_eval = driver.batch_run(
                    n_reference,
                    n_target,
                    n_tri,
                    1,
                    self.header,
                    x_,
                    y_,
                    ww,
                    ropts,
                    iopts,
                    remove_com,
                    sort_by_g,
                    difflim_absolute,
                )
                n_eval += int(numpy.exp(log_n_eval))
            else:
                n_lower = 1
                nrep = (n_tri + n_chunk_ - 1) // n_chunk_
                ret = numpy.empty(n_tri, dtype=dt)
                for i in trange(nrep, *args, **kwargs):
                    l = n_lower - 1
                    u = min([l + n_chunk_, n_tri])
                    ret[l:u], log_n_eval = driver.batch_run(
                        n_reference,
                        n_target,
                        min(n_chunk_, n_tri - n_lower + 1),
                        n_lower,
                        self.header,
                        x_,
                        y_,
                        ww,
                        ropts,
                        iopts,
                        remove_com,
                        sort_by_g,
                        difflim_absolute,
                    )
                    n_lower += n_chunk_
                    n_eval += int(numpy.exp(log_n_eval))
            ret = ret.reshape([n_target, n_reference])
        del driver

        self.n_eval = n_eval

        return ret

    def min_span_tree(
        self,
        x: numpy.ndarray,
        remove_com: bool = True,
        sort_by_g: bool = True,
        verbose: bool = False,
        n_work: int = 0,
        *args,
        **kwargs,
    ) -> networkx.Graph:
        r"""
           Min_span_tree batch calculator.
           Calculate the minimum spanning tree (MST)
           for a sequence of coordinates.

        Args:
            x (npt.NDArray): Reference coordinates $X\in\mathbb R^{d\times n}$, shape [m, n, d].
            verbose(bool): If the calculation takes a long time, display a progress bar.
            remove_com (bool): Remove the center of mass from the coodinates.
            sort_by_g (bool): Sort the reference structures in descending order of self-dispersion and perform the calculation.
                              This affects the computation time required to find an exact solution.
            n_work(int): Maximum working memory size. If <1, then w*n*(n-1)/2.

        Returns:
            ret (networkx.Graph): A minimum spanning tree.
        """

        dt = x.dtype
        x_ = self.to_rank3_coordinates(x)
        n_target = 1 if x.ndim == 2 else x.shape[0]

        driver = _select_driver(self.d, dtype=dt)
        ropts = numpy.array([0.0, 0.0, 0.0], dtype=dt)
        iopts = numpy.array([n_work], dtype=numpy.int32)

        edges, weights, log_n_eval = driver.min_span_tree(
            n_target,
            self.header,
            x_,
            ropts,
            iopts,
            remove_com,
            sort_by_g,
            verbose,
        )
        del driver

        g = networkx.Graph()
        for e, w in zip(edges.T, weights):
            g.add_edge(e[0] - 1, e[1] - 1, weight=w)

        self.n_eval = int(numpy.exp(log_n_eval))

        return g

    def to_rank2_coordinates(self, x: npt.NDArray, dtype=None) -> npt.NDArray:
        """座標が Rank2 であるかのバリデーションを行います.
           Rank2 座標は [self.n_atom, self.d] の次元を持ち, 単一の構造を意味します.


        Args:
           x: reference coordinates for test.
           dtype: Any object that can be interpreted as a numpy data type.
            x: reference coordinates for test.
            dtype: Any object that can be interpreted as a numpy data type.
                   See `numpy.org <https://numpy.org/doc/2.1/reference/arrays.dtypes.html>`__ for detail.
                   default = None.
                  See `numpy.org <https://numpy.org/doc/2.1/reference/arrays.dtypes.html>`__ for detail.
                  default = None.

        Returns:
            Rank2 座標.

        Raises:
            ValueError: Shape が合わない場合。
        """

        if x.ndim != 2:
            raise ValueError

        if x.shape[1] != self.d or x.shape[0] != self.natom:
            raise ValueError

        return numpy.asfortranarray(x, dtype=dtype).flatten()

    def to_rank3_coordinates(self, x: npt.NDArray, dtype=None) -> npt.NDArray:
        """座標が Rank3 であるかのバリデーションを行います.
           Rank3 座標は [nframe, self.n_atom, self.d] の次元を持ち, 構造の系列を意味します.

        Args:
           x: reference coordinates for test.
           dtype: Any object that can be interpreted as a numpy data type.
                  See `numpy.org <https://numpy.org/doc/2.1/reference/arrays.dtypes.html>`__ for detail.
                  default = None.

        Returns:
            Rank3 座標.

        Raises:
            ValueError: Shape が合わない場合。
        """

        if x.ndim == 2:
            pass
            # x_ = x.transpose().reshape((x.shape[1], x.shape[0], 1))
        elif x.ndim == 3:
            pass
            # x_ = x.transpose([2, 1, 0])
        else:
            raise ValueError
        if x.shape[-1] != self.d or x.shape[-2] != self.natom:
            raise ValueError

        return numpy.asfortranarray(x, dtype=dtype).flatten()

    def __del__(self):
        if hasattr(self, "mols"):
            del self.mols
        if hasattr(self, "att"):
            del self.att

    def __str__(self):
        kws = [
            f"mols={self.mols}",
            f"d={self.d}",
        ]
        return "{}({})".format(type(self).__name__, ", ".join(kws))

    def __repr__(self):
        kws = [f"{key}={value!r}" for key, value in self.__dict__.items()]
        return "{}({})".format(type(self).__name__, ", ".join(kws))
