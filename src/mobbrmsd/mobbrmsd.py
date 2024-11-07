# -*- coding: utf-8 -*-
"""mobbrmsd ドライバールーチン

   class mobbrmsd_result:

"""

import numpy
import numpy.typing as npt
import networkx
from tqdm import trange
from .dataclass import molecules, molecular_system


class mobbrmsd_result:
    """mobbrmsd の計算結果と再計算用のメモリを管理するクラス

    mobbrmsd の計算結果、再計算用のメモリを管理するクラス

    """

    def __init__(
        self,
        driver,
        d: int,
        header: npt.NDArray,
        istate: None | npt.NDArray,
        rstate: None | npt.NDArray,
        rot: None | npt.NDArray,
        w: None | npt.NDArray = None,
    ) -> None:
        """initializer

        initializer

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
        :return None
        """

        self.header = header.copy()
        self.d = d

        if (istate is None) or (rstate is None):
            self.istate = numpy.array([0])
            self.rstate = numpy.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
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
        return float(self.rstate[1])

    def lowerbound(self) -> float:
        return float(self.rstate[3])

    def upperbound(self) -> float:
        return float(self.rstate[2])

    def lowerbound_as_rmsd(self) -> float:
        rn = self.rstate[0]
        return numpy.sqrt(
            numpy.max([0.0, rn * (2 * self.lowerbound() + self.autocorr())])
        )

    def upperbound_as_rmsd(self) -> float:
        rn = self.rstate[0]
        return numpy.sqrt(
            numpy.max([0.0, rn * (2 * self.upperbound() + self.autocorr())])
        )

    def sd(self) -> float:
        return float(numpy.max([0.0, (2 * self.upperbound() + self.autocorr())]))

    def msd(self) -> float:
        rn = self.rstate[0]
        return float(rn * self.sd())

    def rmsd(self) -> float:
        return float(numpy.sqrt(self.msd()))

    def bounds(self) -> npt.NDArray:
        return numpy.array([self.lowerbound(), self.upperbound()])

    def bounds_as_rmsd(self) -> npt.NDArray:
        rn = self.rstate[0]
        return rn * (2 * self.bounds() + self.autocorr())

    def n_eval(self) -> int:
        return int(self.rstate[4])

    def log_eval_ratio(self) -> float:
        return self.rstate[5]

    def eval_ratio(self) -> float:
        return numpy.exp(self.log_eval_ratio())

    def is_finished(self) -> bool:
        return self.istate[-1] == (0)

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

    ##
    # @brief restart mobbrmsd
    # @details 途中終了したmobbRMSD を計算する。
    #   RMSD値の他に自己相関部分、BBの上下限、分子の置換インデックス、回転行列などを含む。
    #   所与の終了条件により計算が中断された場合、インスタンスは計算再開用のメモリを保持する。
    # @param
    #   x (numpy.ndarray): 参照構造
    #   y (numpy.ndarray): 対照構造
    #   cutoff (float): BBの下限がcutoff以上になったとき、計算を終了する。 default float(inf)
    #   ub_cutoff (float): BBの上限がcutoff以上になったとき、計算を終了する。 default float(inf)
    #   difflim (float): BBの上限と下限の差がdifflim以下になったとき、計算を終了する。 default 0
    #   maxeval (int): BBのノード評価数がmaxevalを超えたとき、計算を終了する。
    #                  ただし、最低でも expand と closure の1サイクルは実行される。
    #                  maxeval < 0 のとき、無制限。
    #                  default -1
    #   difflim_absolute (bool): True なら difflim を RMSD 換算の絶対値で用いる。 default False
    # @return None
    def restart(
        self,
        cutoff: float = float("inf"),
        ub_cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        difflim_absolute: bool = False,
        get_rotation: bool = False,
    ) -> None:

        if not hasattr(self, "w"):
            return

        driver = select_driver(self.d, dtype=self.w.dtype)
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
    # @details swap and rotation target coordinate
    # @param
    #   y (numpy.ndarray): 対象構造
    # @return None
    def rotate_y(
        self,
        y: npt.NDArray,
    ) -> None:
        y_ = y.flatten()
        driver = select_driver(self.d, dtype=y.dtype)
        driver.rotate_y(self.header, self.istate, self.rstate, self.rot, y_)
        del driver
        return y_.reshape(y.shape)


def select_driver(d: int, dtype=None):
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


def select_dtype(x: npt.NDArray, y: npt.NDArray):
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


def varidation_coordinates_1(x: npt.NDArray, d, natom, dtype=None) -> npt.NDArray:

    if x.ndim != 2:
        raise ValueError

    if x.shape[1] != d or x.shape[0] != natom:
        raise ValueError

    return numpy.asfortranarray(x, dtype=dtype).flatten()


def varidation_coordinates_2(x: npt.NDArray, d, natom, dtype=None) -> npt.NDArray:

    if x.ndim == 2:
        pass
        # x_ = x.transpose().reshape((x.shape[1], x.shape[0], 1))
    elif x.ndim == 3:
        pass
        # x_ = x.transpose([2, 1, 0])
    else:
        raise ValueError
    if x.shape[-1] != d or x.shape[-2] != natom:
        raise ValueError

    return numpy.asfortranarray(x, dtype=dtype).flatten()


##
# @class mobbrmsd
# @brief mobbrmsdユーザーインターフェース
# @details molecular-oriented RMSD with branch-and-bound.


class mobbrmsd:
    ##
    # @brief initializer
    # @details initializer
    # @param mols (molecules | molecular_system): 系の分子
    #        d(int): Spatial dimension. default=3.
    # @return None
    def __init__(
        self,
        mols: molecules | molecular_system = molecules(1, 1),
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

        driver = select_driver(d, dtype=None)
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

    ##
    # @brief rmsd
    # @details RMSD を計算する。
    # @param
    #   x (numpy.ndarray): 参照構造
    #   y (numpy.ndarray): 対照構造
    # @return float
    def rmsd(
        self,
        x: npt.NDArray,
        y: npt.NDArray,
    ) -> float:
        dt = select_dtype(x, y)
        x_ = varidation_coordinates_1(x, self.d, self.natom, dtype=dt)
        y_ = varidation_coordinates_1(y, self.d, self.natom, dtype=dt)
        driver = select_driver(self.d, dtype=dt)
        w = numpy.empty(self.memsize, dtype=dt)
        ropts = numpy.array([float("inf"), 0.0], dtype=dt)
        iopts = numpy.array([-1], dtype=numpy.int32)

        _, rret, _ = driver.run(
            self.n_int,  # n_int
            self.n_float,  # n_float
            self.n_rot,  # n_rot
            self.header,  # header
            x_,  # X
            y_,  # Y
            w,  # W
            ropts,  # ropts
            iopts,  # iopts
            True,  # remove_com
            True,  # sort_by_g
            False,  # difflim_absolute
            False,  # rotate_y
            False,  # get_rotation
        )  # returns (int_states, float_states, rotation)

        ret = driver.rmsd(rret)
        del driver, w, ropts, iopts
        return ret

    ##
    # @brief general runner
    # @details RMSD を計算する。
    #   RMSD の他に自己相関、BBの上下限、分子の置換インデックス、回転行列、再計算用情報を返す。
    # @param
    #   x (numpy.ndarray): 参照構造
    #   y (numpy.ndarray): 対照構造
    #   cutoff (float): BBの下限がcutoff以上になったとき、計算を終了する。 default float(inf)
    #   ub_cutoff (float): BBの上限がcutoff以上になったとき、計算を終了する。 default float(inf)
    #   difflim (float): BBの上限と下限の差がdifflim以下になったとき、計算を終了する。 default 0
    #   maxeval (int): BBのノード評価数がmaxevalを超えたとき、計算を終了する。
    #                  ただし、最低でも expand と closure の1サイクルは実行される。
    #                  maxeval < 0 のとき、無制限。
    #                  default -1
    #   remove_com (bool): 参照構造と対照構造から重心を除去する。 default True
    #   sort_by_g (bool): 参照構造を自己分散の大きい順に並び替えて計算を実行する。 default True
    #   rotate_y (bool): 対象構造に対して置換と回転を実行する。 default False
    #   get_rotation (bool): 回転行列を計算する。 default False
    # @return mobbrmsd_result
    #   所与の終了条件により計算が中断された場合、mobbrmsd_result はインスタンスは計算再開用のデータを保持する。
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

        dt = select_dtype(x, y)
        x_ = varidation_coordinates_1(x, self.d, self.natom, dtype=dt)
        y_ = varidation_coordinates_1(y, self.d, self.natom, dtype=dt)
        driver = select_driver(self.d, dtype=dt)
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

    ##
    # @brief batch runner
    # @details 一連の座標について mobbRMSD を計算する。
    #   構造 x のみが与えられたとき、x-x の間の RMSD 対称行列を返す。
    #   構造 x, y が与えられたとき、x-y の間の RMSD 行列を返す。
    # @param
    #   x (numpy.ndarray): 構造
    #   y (numpy.ndarray): 対照構造, optional
    #   cutoff (float): BBの下限がcutoff以上になったとき、計算を終了する。 default float(inf)
    #   ub_cutoff (float): BBの上限がcutoff以上になったとき、計算を終了する。 default float(inf)
    #   difflim (float): BBの上限と下限の差が (G/2) * difflim 以下になったとき、計算を終了する。 default 0.0
    #   maxeval (int): BBのノード評価数がmaxevalを超えたとき、計算を終了する。
    #                  ただし、最低でも expand と closure の1サイクルは実行される。
    #                  maxeval < 0 のとき、無制限。
    #                  default -1
    #   remove_com (bool): 参照構造と対照構造から重心を除去する。 default True
    #   sort_by_g (bool): 参照構造を自己分散の大きい順に並び替えて計算を実行する。 default True
    #   difflim_absolute (bool): True なら difflim を RMSD 換算の絶対値で用いる。 default False
    #   verbose (bool): 進捗を表示する。 default True
    #   n_chunk (int): 一度に計算される回数。 <1 の場合一括計算。 default 1000
    # @return npt.NDArray
    #   RMSD 行列
    def batch_run(
        self,
        x: npt.NDArray,
        y: None | npt.NDArray = None,
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

        if y is None:
            dt = x.dtype
            x_ = varidation_coordinates_2(x, self.d, self.natom)
            driver = select_driver(self.d, dtype=dt)
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
            dt = select_dtype(x, y)
            x_ = varidation_coordinates_2(x, self.d, self.natom, dtype=dt)
            y_ = varidation_coordinates_2(y, self.d, self.natom, dtype=dt)
            driver = select_driver(self.d, dtype=dt)
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

    ##
    # @brief min_span_tree runner
    # @details 一連の座標について mobbRMSD の最小全域木を計算する。
    # @param
    #   x (numpy.ndarray): 参照構造
    #   cutoff (float): BBの下限がcutoff以上になったとき、計算を終了する。 default float(inf)
    #   ub_cutoff (float): BBの上限がcutoff以上になったとき、計算を終了する。 default float(inf)
    #   difflim (float): BBの上限と下限の差が difflim 以下になったとき、計算を終了する。 default 0
    #   maxeval (int): BBのノード評価数がmaxevalを超えたとき、計算を終了する。
    #                  ただし、最低でも expand と closure の1サイクルは実行される。
    #                  maxeval < 0 のとき、無制限。
    #                  default -1
    #   remove_com (bool): 参照構造と対照構造から重心を除去する。 default True
    #   sort_by_g (bool): 参照構造を自己分散の大きい順に並び替えて計算を実行する。 default True
    #   difflim_absolute (bool): True なら difflim を RMSD 換算の絶対値で用いる。 default False
    #   verbose (bool): 進捗を表示する。 default True
    #   n_chunk (int): 一度に計算される回数。 <1 の場合一括計算。 default 1000
    # @return networkx.Graph
    def min_span_tree(
        self,
        x: numpy.ndarray,
        cutoff: float = float("inf"),
        ub_cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        remove_com: bool = True,
        sort_by_g: bool = True,
        difflim_absolute: bool = False,
        verbose: bool = False,
        n_chunk: int = 0,
        *args,
        **kwargs,
    ) -> networkx.Graph:

        dt = x.dtype
        x_ = varidation_coordinates_2(x, self.d, self.natom)
        n_target = 1 if x.ndim == 2 else x.shape[0]
        ww = numpy.empty((n_target * (n_target - 1) // 2 * self.memsize), dtype=dt)

        driver = select_driver(self.d, dtype=dt)
        ropts = numpy.array([cutoff, ub_cutoff, difflim], dtype=dt)
        iopts = numpy.array([maxeval], dtype=numpy.int32)

        edges, weights = driver.min_span_tree(
            n_target,
            self.header,
            x_,
            ww,
            ropts,
            iopts,
            remove_com,
            sort_by_g,
            difflim_absolute,
        )
        del driver

        g = networkx.Graph()
        for e, w in zip(edges.T, weights):
            g.add_edge(e[0] - 1, e[1] - 1, weight=w)

        return g

    def __del__(self):
        if hasattr(self, "molecules"):
            del self.molecules
        if hasattr(self, "att"):
            del self.att

    def __str__(self):
        kws = [
            f"molecules={self.molecules}",
            f"d={self.d}",
        ]
        return "{}({})".format(type(self).__name__, ", ".join(kws))

    def __repr__(self):
        kws = [f"{key}={value!r}" for key, value in self.__dict__.items()]
        return "{}({})".format(type(self).__name__, ", ".join(kws))
