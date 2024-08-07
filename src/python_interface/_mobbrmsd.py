import numpy
import numpy.typing as npt
import networkx
from tqdm import trange
import dataclasses


@dataclasses.dataclass(frozen=True)
class DataclassMolecule:
    n_apm: int
    n_mol: int
    sym: any = None
    name: str = ""


class mobbrmsd_result:
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
        self.header = header.copy()
        if (istate is None) or (rstate is None):
            self.state = None
            self.rmsd = 0.0
            self.autocorr = 0.0
            self.bounds = (0.0, 0.0)
            self.n_eval = 0
            self.log_eval_ratio = None
            self.eval_ratio = 0.0
            self.is_finished = True
            return
        self.d = d
        self.istate = istate.copy()
        self.rstate = rstate.copy()
        self.rmsd = driver.rmsd(rstate)
        self.autocorr = driver.autocorr(rstate)
        self.bounds = driver.bounds(rstate)
        self.n_eval = driver.n_eval(rstate)
        self.log_eval_ratio = driver.log_eval_ratio(rstate)
        self.eval_ratio = numpy.exp(self.log_eval_ratio)
        self.is_finished = driver.is_finished(header, istate, rstate)
        if rot is None:
            self.rot = numpy.zeros([self.d * self.d])
        else:
            self.rot = rot.copy()
        if w is not None:
            self.w = w.copy()

    def __repr__(self):
        kws = [f"{key}={value!r}" for key, value in self.__dict__.items()]
        return "{}({})".format(type(self).__name__, ", ".join(kws))

    ##
    # @brief restart mobbrmsd
    # @details 途中終了したmobbRMSD を計算する。
    #   RMSD値の他に自己相関部分、BBの上下限、分子の置換インデックス、回転行列などを含む。
    #   所与の終了条件により計算が中断された場合、インスタンスは計算再開用のメモリを保持する。
    # @param
    #   x (numpy.ndarray): 参照構造
    #   y (numpy.ndarray): 対照構造
    #   cutoff (float): BBの上限がcutoff以下になったとき、計算を終了する。 default float(inf)
    #   difflim (float): BBの上限と下限の差がdifflim以下になったとき、計算を終了する。 default 0
    #   maxeval (int): BBのノード評価数がmaxevalを超えたとき、計算を終了する。
    #                  ただし、最低でも expand と closure の1サイクルは実行される。
    #                  maxeval < 0 のとき、無制限。
    #                  default -1
    # @return mobbrmsd_result
    def restart(
        self,
        cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        get_rotation: bool = False,
    ) -> None:

        if not hasattr(self, "w"):
            return

        driver = select_driver(self.d, dtype=self.w.dtype)
        driver.restart(
            self.header,
            self.istate,
            self.rstate,
            self.w,
            self.rot,
            cutoff,
            difflim,
            maxeval,
            get_rotation,
        )

        self.rmsd = driver.rmsd(self.rstate)
        self.autocorr = driver.autocorr(self.rstate)
        self.bounds = driver.bounds(self.rstate)
        self.n_eval = driver.n_eval(self.rstate)
        self.log_eval_ratio = driver.log_eval_ratio(self.rstate)
        self.eval_ratio = numpy.exp(self.log_eval_ratio)
        self.is_finished = driver.is_finished(self.header, self.istate, self.rstate)

        del driver

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
    if d == 2:
        if dt == numpy.float64:
            from .mobbrmsd_2ddp import driver
        elif dt == numpy.float32:
            from .mobbrmsd_2dsp import driver
        else:
            raise ValueError
    elif d == 3:
        if dt == numpy.float64:
            from .mobbrmsd_3ddp import driver
        elif dt == numpy.float32:
            from .mobbrmsd_3dsp import driver
        else:
            raise ValueError
    elif d == 1 or d > 3:
        if dt == numpy.float64:
            from .mobbrmsd_xddp import driver
        elif dt == numpy.float32:
            from .mobbrmsd_xdsp import driver
        else:
            raise ValueError

        driver.setup_dimension_(d)
    else:
        raise ValueError
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

    if dtype is None:
        return x.flatten()
    elif x.dtype == dtype:
        return x.flatten()
    else:
        return numpy.astype(x, dtype=dtype).flatten()


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

    if dtype is None:
        return x.flatten()
    elif x.dtype == dtype:
        return x.flatten()
    else:
        return numpy.astype(x, dtype=dtype).flatten()


##
# @class mobbrmsd
# @brief mobbrmsdユーザーインターフェース
# @details molecular-oriented RMSD with branch-and-bound.


class mobbrmsd:
    ##
    # @brief initializer
    # @details initializer
    # @param molecules (DataclassMolecule か DataclassMolecule のリスト): 系の分子
    #        d(int): Spatial dimension. default=3.
    # @return None
    def __init__(
        self,
        molecules: DataclassMolecule | list[DataclassMolecule] = DataclassMolecule(
            1, 1
        ),
        d: int = 3,
    ) -> None:

        self.molecules = molecules
        molecular_sequence = []

        def add_molecule(mol: DataclassMolecule):
            ret = [mol.n_apm, mol.n_mol]
            if mol.sym is None:
                ret += [1]
            else:
                sym = numpy.array(mol.sym, dtype=numpy.int32) + 1
                if sym.ndim == 1:
                    if sym.shape[0] != 0 and sym.shape[0] != mol.n_apm:
                        raise ValueError
                elif sym.ndim == 2:
                    if sym.shape[1] != mol.n_apm:
                        raise ValueError
                ret += [sym.shape[0] + 1] + sym.flatten().tolist()
            return ret

        if type(molecules) is list:
            for mol in molecules:
                molecular_sequence += add_molecule(mol)
        elif type(molecules) is DataclassMolecule:
            molecular_sequence += add_molecule(molecules)
        else:
            raise ValueError

        driver = select_driver(d, dtype=None)

        (
            self.d,
            self.natom,
            self.memsize,
            self.njob,
            self.n_header,
            self.n_int,
            self.n_float,
            self.n_rot,
        ) = driver.decode_attributes(molecular_sequence)
        self.header = driver.decode_header(molecular_sequence, self.n_header)
        del driver

    ##
    # @brief mobbRMSD
    # @details mobbRMSD を計算する。
    # @param
    #   x (numpy.ndarray): 参照構造
    #   y (numpy.ndarray): 対照構造
    # @return float
    def mobbrmsd(
        self,
        x: npt.NDArray,
        y: npt.NDArray,
    ) -> float:
        dt = select_dtype(x, y)
        x_ = varidation_coordinates_1(x, self.d, self.natom, dtype=dt)
        y_ = varidation_coordinates_1(y, self.d, self.natom, dtype=dt)
        driver = select_driver(self.d, dtype=dt)
        w = numpy.empty(self.memsize, dtype=dt)

        _, rret, _ = driver.run(
            self.n_int,
            self.n_float,
            self.n_rot,
            self.header,
            x_,
            y_,
            w,
            float("inf"),
            0.0,
            -1,
            True,
            True,
            False,
            False,
        )

        ret = driver.rmsd(rret)
        del driver
        del w
        return ret

    ##
    # @brief general runner
    # @details mobbRMSD を計算する。
    #   RMSD の他に自己相関、BBの上下限、分子の置換インデックス、回転行列、再計算用情報を返す。
    # @param
    #   x (numpy.ndarray): 参照構造
    #   y (numpy.ndarray): 対照構造
    #   cutoff (float): BBの上限がcutoff以下になったとき、計算を終了する。 default float(inf)
    #   difflim (float): BBの上限と下限の差がdifflim以下になったとき、計算を終了する。 default 0
    #   maxeval (int): BBのノード評価数がmaxevalを超えたとき、計算を終了する。
    #                  ただし、最低でも expand と closure の1サイクルは実行される。
    #                  maxeval < 0 のとき、無制限。
    #                  default -1
    #   remove_com (bool): 参照構造と対照構造から重心を除去する。 default True
    #   sort_by_g (bool): 参照構造を自己分散の大きい順に並び替えて計算を実行する。 default True
    #   rotate_y (bool): 対象構造に対して置換と回転を実行する。 default False
    #   rotate_y (bool): 対象構造に対して置換と回転を実行する。 default False
    # @return mobbrmsd_result
    #   所与の終了条件により計算が中断された場合、mobbrmsd_result はインスタンスは計算再開用のデータを保持する。
    def run(
        self,
        x: npt.NDArray,
        y: npt.NDArray,
        cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        remove_com: bool = True,
        sort_by_g: bool = True,
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

        iret, rret, rot = driver.run(
            self.n_int,
            self.n_float,
            self.n_rot,
            self.header,
            x_,
            y_,
            w,
            cutoff,
            difflim,
            maxeval,
            remove_com,
            sort_by_g,
            rotate_y,
            get_rotation,
        )

        ret = mobbrmsd_result(driver, self.d, self.header, iret, rret, rot, w=w)
        del driver
        del w
        return ret

    ##
    # @brief batch runner
    # @details 一連の座標について mobbRMSD を計算する。
    #   構造 x のみが与えられたとき、x-x の間の RMSD 対称行列を返す。
    #   構造 x, y が与えられたとき、x-y の間の RMSD 行列を返す。
    # @param
    #   x (numpy.ndarray): 構造
    #   y (numpy.ndarray): 対照構造, optional
    #   cutoff (float): BBの上限がcutoff以下になったとき、計算を終了する。 default float(inf)
    #   difflim (float): BBの上限と下限の差がdifflim以下になったとき、計算を終了する。 default 0
    #   maxeval (int): BBのノード評価数がmaxevalを超えたとき、計算を終了する。
    #                  ただし、最低でも expand と closure の1サイクルは実行される。
    #                  maxeval < 0 のとき、無制限。
    #                  default -1
    #   remove_com (bool): 参照構造と対照構造から重心を除去する。 default True
    #   sort_by_g (bool): 参照構造を自己分散の大きい順に並び替えて計算を実行する。 default True
    #   verbose (bool): 進捗を表示する。 default True
    #   n_chunk (int): 一度に計算される回数。 <1 の場合一括計算。 default 1000
    # @return npt.NDArray
    #   RMSD 行列
    def batch_run(
        self,
        x: npt.NDArray,
        y: None | npt.NDArray = None,
        cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        remove_com: bool = True,
        sort_by_g: bool = True,
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

            if n_tri == n_chunk_ or not verbose:
                r_tri = driver.batch_run_tri(
                    n_target,
                    n_tri,
                    1,
                    self.header,
                    x_,
                    ww,
                    cutoff,
                    difflim,
                    maxeval,
                    remove_com,
                    sort_by_g,
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
                        cutoff,
                        difflim,
                        maxeval,
                        remove_com,
                        sort_by_g,
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
                    cutoff,
                    difflim,
                    maxeval,
                    remove_com,
                    sort_by_g,
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
                        cutoff,
                        difflim,
                        maxeval,
                        remove_com,
                        sort_by_g,
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
    #   cutoff (float): BBの上限がcutoff以下になったとき、計算を終了する。 default float(inf)
    #   difflim (float): BBの上限と下限の差がdifflim以下になったとき、計算を終了する。 default 0
    #   maxeval (int): BBのノード評価数がmaxevalを超えたとき、計算を終了する。
    #                  ただし、最低でも expand と closure の1サイクルは実行される。
    #                  maxeval < 0 のとき、無制限。
    #                  default -1
    #   remove_com (bool): 参照構造と対照構造から重心を除去する。 default True
    #   sort_by_g (bool): 参照構造を自己分散の大きい順に並び替えて計算を実行する。 default True
    #   verbose (bool): 進捗を表示する。 default True
    #   n_chunk (int): 一度に計算される回数。 <1 の場合一括計算。 default 1000
    # @return networkx.Graph
    def min_span_tree(
        self,
        x: numpy.ndarray,
        cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        remove_com: bool = True,
        sort_by_g: bool = True,
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
        edges, weights = driver.min_span_tree(
            n_target,
            self.n_int,
            self.n_float,
            self.header,
            x_,
            ww,
            cutoff,
            difflim,
            maxeval,
            remove_com,
            sort_by_g,
        )
        del driver

        g = networkx.Graph()
        for e, w in zip(edges.T, weights):
            g.add_edge(e[0] - 1, e[1] - 1, weight=w)

        return g

    def __del__(self):
        del self.molecules
        del self.d
        del self.natom
        del self.memsize
        del self.njob
        del self.n_header
        del self.n_int
        del self.n_float
        del self.header

    def __str__(self):
        kws = [
            f"molecules={self.molecules}",
            f"d={self.d}",
        ]
        return "{}({})".format(type(self).__name__, ", ".join(kws))

    def __repr__(self):
        kws = [f"{key}={value!r}" for key, value in self.__dict__.items()]
        return "{}({})".format(type(self).__name__, ", ".join(kws))
