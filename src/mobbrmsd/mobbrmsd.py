# -*- coding: utf-8 -*-
"""mobbrmsd driver routine

"""

import numpy
import numpy.typing as npt
import networkx
from typing import Union
from tqdm import trange
from .dataclass import molecules, molecular_system


class mobbrmsd_result:
    """mobbrmsd の計算結果と再計算用のメモリを管理するクラス.
       ほとんどの場合, ユーザーはインスタンスを生成する必要はありません.

    :param driver: fortran driver
    :param d: spatial dimension
    :type d: int
    :param header: header array
    :type header: npt.NDArray
    :param istate: state array 1
    :type istate: None | npt.NDArray
    :param rstate: state array 2
    :type rstate: None | npt.NDArray (float32/float64)
    :param rot: rotation matrix
    :type rot: None | npt.NDArray, shape[d,d] (float32/float64)
    :param w: working memory
    :type w: None | npt.NDArray (float32/float64)
    """

    def __init__(
        self,
        driver,
        d: int,
        header: npt.NDArray,
        istate: Union[None, npt.NDArray],
        rstate: Union[None, npt.NDArray],
        rot: Union[None, npt.NDArray],
        w: Union[None, npt.NDArray] = None,
        # istate: None | npt.NDArray,
        # rstate: None | npt.NDArray,
        # rot: None | npt.NDArray,
        # w: None | npt.NDArray = None,
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
        """自己相関項.

        :return: 自己相関項.
        :rtype: float
        """
        return float(self.rstate["""INDEX_TO_AUTOCORR"""])

    def lowerbound(self) -> float:
        """目的関数の変分下限.

        :return: 変分下限.
        :rtype: float
        """
        return float(self.rstate["""INDEX_TO_LOWERBOUND"""])

    def upperbound(self) -> float:
        """目的関数の上限.

        :return: 変分上限.
        :rtype: float
        """
        return float(self.rstate["""INDEX_TO_UPPERBOUND"""])

    def lowerbound_as_rmsd(self) -> float:
        """RMSD 換算された変分下限.

        :return: RMSD 換算された変分下限.
        :rtype: float
        """
        rn = self.rstate["""RECIPROCAL_OF_N"""]
        return numpy.sqrt(
            numpy.max([0.0, rn * (2 * self.lowerbound() + self.autocorr())])
        )

    def upperbound_as_rmsd(self) -> float:
        """RMSD 換算された上限.

        :return: RMSD 換算された上限.
        :rtype: float
        """
        rn = self.rstate["""RECIPROCAL_OF_N"""]
        return numpy.sqrt(
            numpy.max([0.0, rn * (2 * self.upperbound() + self.autocorr())])
        )

    def sd(self) -> float:
        """二乗変位.
           計算が早期終了されている場合, これは暫定値で, 真の解よりも大きい可能性があります。

        :return: 二乗変位.
        :rtype: float
        """
        ret = float(2 * self.upperbound() + self.autocorr())
        if ret < 0.0:
            return 0.0
        else:
            return ret

    def msd(self) -> float:
        """平均二乗変位.
           計算が早期終了されている場合, これは暫定値で, 真の解よりも大きい可能性があります。

        :return: 平均二乗変位.
        :rtype: float
        """
        rn = self.rstate["""RECIPROCAL_OF_N"""]
        return float(rn * self.sd())

    def rmsd(self) -> float:
        """RMSD.
           計算が早期終了されている場合, これは暫定値で, 真の解よりも大きい可能性があります。

        :return: RMSD.
        :rtype: float
        """
        return float(numpy.sqrt(self.msd()))

    def bounds(self) -> npt.NDArray:
        """目的関数の厳密解が含まれる区間.

        :return: lowerbound, upperbound
        :rtype: npt.NDArray
        """
        return numpy.array([self.lowerbound(), self.upperbound()])

    def bounds_as_rmsd(self) -> npt.NDArray:
        """RMSD 換算の厳密解が含まれる区間.

        :return: lowerbound, upperbound
        :rtype: npt.NDArray
        """
        rn = self.rstate["""RECIPROCAL_OF_N"""]
        return rn * (2 * self.bounds() + self.autocorr())

    def n_eval(self) -> int:
        """計算回数.

        :return: 計算回数.
        :rtype: int
        """
        return int(self.rstate["""INDEX_TO_N_EVAL"""])

    def log_eval_ratio(self) -> float:
        """全探索に対する、計算回数の割合の対数.

        :return: 計算回数の割合の対数.
        :rtype: float
        """
        return self.rstate["""INDEX_TO_LOG_RATIO"""]

    def eval_ratio(self) -> float:
        """全探索に対する、計算回数の割合.

        :return: 計算回数の割合.
        :rtype: float
        """
        return numpy.exp(self.log_eval_ratio())

    def is_finished(self) -> bool:
        """探索が最後まで完了しているか.

        :return: 探索が完了しているならTrue.
        :rtype: bool
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
        """途中終了した計算を再開します.
           計算結果はインスタンスに上書きされます.

        :param cutoff: BBの下限がRMSD換算でcutoff以上になったとき, 計算を終了する. default=float(inf).
        :type cutoff: float
        :param ub_cutoff: BBの上限がRMSD換算でub_cutoff以上であれば, 計算を終了する. default=float(inf).
        :type ub_cutoff: float
        :param difflim: BBの上限と下限の差が difflim 以下になったとき,計算を終了する. default=0.0.
        :type difflim: float
        :param maxeval: BBのノード評価数がmaxevalを超えたとき,計算を終了する.
                        ただし,最低でも expand と closure の1サイクルは実行される.
                        maxeval < 0 のとき,無制限.
                        default=-1.
        :type maxeval: int
        :param difflim_absolute: True なら difflim を RMSD 換算で用いる. default=False.
        :type difflim_absolute: bool
        :param get_rotation: 回転行列を計算する. default=False.
        :type get_rotation: bool
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

    ##
    # @brief swap and rotate Y
    # @details
    # @param
    #   y (numpy.ndarray): 対象構造
    # @return None
    def rotate_y(
        self,
        y: npt.NDArray,
    ) -> None:
        """swap and rotate Y.
           swap and rotation target coordinate.

        :param y: target coordinates.
        :type y: npt.NDArray
        """
        y_ = y.flatten()
        driver = _select_driver(self.d, dtype=y.dtype)
        driver.rotate_y(self.header, self.istate, self.rstate, self.rot, y_)
        del driver
        return y_.reshape(y.shape)

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
    """mobbrmsd driver

    :param mols: molecules/molecular_system specifier
    :type mols: molecules | molecular_system
    :param d: Spatial dimension. default=3.
    :type d: int
    """

    def __init__(
        self,
        mols: Union[molecules, molecular_system] = molecules(1, 1),
        # mols: molecules | molecular_system = molecules(1, 1),
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
        """calculate RMSD of a structural pair. (Simplified interface)
           RMSDを計算します.

        :param x: reference coordinates, rank 2.
        :type x: npt.NDArray
        :param y: target coordinates, rank 2.
        :type y: npt.NDArray
        :return: rmsd value.
        :rtype: float
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
        """calculate RMSD of a structural pair.
           RMSD, 自己相関, 上下限, 分子の置換インデックス, 回転行列 (optional)を計算します.
           計算が早期終了された場合, mobbrmsd_result はインスタンスは計算再開用のデータを保持します.

        :param x: reference coordinates, rank 2.
        :type x: npt.NDArray
        :param y: target coordinates, rank 2.
        :type y: npt.NDArray
        :param cutoff: BBの下限がRMSD換算でcutoff以上になったとき, 計算を終了する. default=float(inf).
        :type cutoff: float
        :param ub_cutoff: BBの上限がRMSD換算でub_cutoff以上であれば, 計算を終了する. default=float(inf).
        :type ub_cutoff: float
        :param difflim: BBの上限と下限の差が difflim 以下になったとき,計算を終了する. default=0.0.
        :type difflim: float
        :param maxeval: BBのノード評価数がmaxevalを超えたとき,計算を終了する.
                        ただし,最低でも expand と closure の1サイクルは実行される.
                        maxeval < 0 のとき,無制限.
                        default=-1.
        :type maxeval: int
        :param remove_com: 参照構造と対照構造から重心を除去する. default=True.
        :type remove_com: bool
        :param sort_by_g: 参照構造を自己分散の大きい順に並び替えて計算を実行する. default=True.
        :type sort_by_g: bool
        :param difflim_absolute: True なら difflim を RMSD 換算で用いる. default=False.
        :type difflim_absolute: bool
        :param rotate_y: 対象構造に対して置換と回転を実行する. default=False.
        :type rotate_y: bool
        :param get_rotation: 回転行列を計算する. default=False.
        :type get_rotation: bool
        :return: RMSD, 自己相関, 上下限, 分子の置換インデックス, 回転行列 (optional), 計算再開用データ(早期終了時).
        :rtype: mobbrmsd_result
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
        y: Union[None, npt.NDArray] = None,
        # y: None | npt.NDArray = None,
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
        """Batch RMSD runnner for multiple structures.
           座標系列に対して mobbRMSD 行列を計算します.
           構造 x のみが与えられたとき, RMSD 対称行列 D(x, x) を返します.
           構造 x, y が与えられたとき, RMSD 行列 D(x, y) を返します.

        :param x: reference coordinates, rank 2/3.
        :type x: npt.NDArray
        :param y: target coordinates, rank 2/3.
        :type y: npt.NDArray
        :param cutoff: BBの下限がRMSD換算でcutoff以上になったとき, 計算を終了する. default=float(inf).
        :type cutoff: float
        :param ub_cutoff: BBの上限がRMSD換算でub_cutoff以上であれば, 計算を終了する. default=float(inf).
        :type ub_cutoff: float
        :param difflim: BBの上限と下限の差が difflim 以下になったとき,計算を終了する. default=0.0.
        :type difflim: float
        :param maxeval: BBのノード評価数がmaxevalを超えたとき,計算を終了する.
                        ただし,最低でも expand と closure の1サイクルは実行される.
                        maxeval < 0 のとき,無制限.
                        default=-1.
        :type maxeval: int
        :param remove_com: 参照構造と対照構造から重心を除去する. default=True.
        :type remove_com: bool
        :param sort_by_g: 参照構造を自己分散の大きい順に並び替えて計算を実行する. default=True.
        :type sort_by_g: bool
        :param difflim_absolute: True なら difflim を RMSD 換算で用いる. default=False.
        :type difflim_absolute: bool
        :param rotate_y: 対象構造に対して置換と回転を実行する. default=False.
        :type rotate_y: bool
        :param verbose: 計算が長くなる場合, 進捗バーを表示する. default=True.
        :type verbose: bool
        :param n_chunk: 一度にまとめて計算されるバッチサイズ上限. <1 の場合, 一括計算. default=1000.
        :type n_chunk: int
        :return: RMSD 行列
        :rtype: npt.NDArray
        """

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
                r_tri = driver.batch_run_tri(
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
            else:
                n_lower = 1
                nrep = (n_tri + n_chunk_ - 1) // n_chunk_
                r_tri = numpy.empty(n_tri, dtype=dt)
                for i in trange(nrep, *args, **kwargs):
                    l = n_lower - 1
                    u = min([l + n_chunk_, n_tri])
                    r_tri[l:u] = driver.batch_run_tri(
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
                ret = driver.batch_run(
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
            else:
                n_lower = 1
                nrep = (n_tri + n_chunk_ - 1) // n_chunk_
                ret = numpy.empty(n_tri, dtype=dt)
                for i in trange(nrep, *args, **kwargs):
                    l = n_lower - 1
                    u = min([l + n_chunk_, n_tri])
                    ret[l:u] = driver.batch_run(
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
            ret = ret.reshape([n_target, n_reference])
        del driver
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
        """Min_span_tree batch runner.
           座標の系列について mobbRMSD の最小全域木を計算します.

        :param x: reference coordinates, rank 3.
        :type x: npt.NDArray
        :param remove_com: 参照構造と対照構造から重心を除去する. default=True.
        :type remove_com: bool
        :param sort_by_g: 参照構造を自己分散の大きい順に並び替えて計算を実行する. default=True.
        :type sort_by_g: bool
        :param verbose: 計算が長くなる場合, 進捗バーを表示する. default=True.
        :type verbose: bool
        :param n_work: メモリサイズ上限. <1 の場合、n*(n-1)/2. default=None.
        :type n_work: int
        :return: 最小全域木
        :rtype: networkx.Graph
        """

        dt = x.dtype
        x_ = self.to_rank3_coordinates(x)
        n_target = 1 if x.ndim == 2 else x.shape[0]

        driver = _select_driver(self.d, dtype=dt)
        ropts = numpy.array([0.0, 0.0, 0.0], dtype=dt)
        iopts = numpy.array([n_work], dtype=numpy.int32)

        edges, weights = driver.min_span_tree(
            n_target,
            self.header,
            x_,
            ropts,
            iopts,
            remove_com,
            sort_by_g,
        )
        del driver

        g = networkx.Graph()
        for e, w in zip(edges.T, weights):
            g.add_edge(e[0] - 1, e[1] - 1, weight=w)

        return g

    def to_rank2_coordinates(self, x: npt.NDArray, dtype=None) -> npt.NDArray:
        """座標が Rank2 であるかのバリデーションを行います.
           Rank2 座標は [self.n_atom, self.d] の次元を持ち, 単一の構造を意味します.

        :param x: reference coordinates for test.
        :type x: npt.NDArray
        :param dtype: Any object that can be interpreted as a numpy data type.
                      See `numpy.org <https://numpy.org/doc/2.1/reference/arrays.dtypes.html>`__ for detail.
                      default = None.
        :return: Rank2 座標.
        :rtype: npt.NDArray
        :raise: ValueError
        """

        if x.ndim != 2:
            raise ValueError

        if x.shape[1] != self.d or x.shape[0] != self.natom:
            raise ValueError

        return numpy.asfortranarray(x, dtype=dtype).flatten()

    def to_rank3_coordinates(self, x: npt.NDArray, dtype=None) -> npt.NDArray:
        """座標が Rank3 であるかのバリデーションを行います.
           Rank3 座標は [nframe, self.n_atom, self.d] の次元を持ち, 構造の系列を意味します.

        :param x: reference coordinates for test.
        :type x: npt.NDArray
        :param dtype: Any object that can be interpreted as a numpy data type.
                      See `numpy.org <https://numpy.org/doc/2.1/reference/arrays.dtypes.html>`__ for detail.
                      default = None.
        :return: Rank3 座標.
        :rtype: npt.NDArray
        :raise: ValueError
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
