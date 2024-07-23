import numpy
import numpy.typing as npt
import networkx
from tqdm import tqdm
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
    ) -> None:

        if not hasattr(self, "w"):
            return

        driver = select_driver(self.d, dtype=self.w.dtype)
        driver.restart(
            self.header,
            self.istate,
            self.rstate,
            self.w,
            cutoff,
            difflim,
            maxeval,
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
        driver.rotate_y(self.header, self.istate, self.rstate, y_)
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
            from .mobbrmsd_gddp import driver
        elif dt == numpy.float32:
            from .mobbrmsd_gdsp import driver
        else:
            raise ValueError

        driver.setup_dimension(d)
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


def varidation_coordinates_2(x: npt.NDArray) -> npt.NDArray:

    if x.ndim == 2:
        x_ = x.transpose().reshape((x.shape[1], x.shape[0], 1))
    elif x.ndim == 3:
        x_ = x.transpose([2, 1, 0])
    else:
        raise ValueError
    if x_.shape[0] != d or x_.shape[1] != natom:
        raise ValueError

    if dtype is None:
        return x_.flatten()
    elif x.dtype == dtype:
        return x_.flatten()
    else:
        return numpy.astype(x_, dtype=dtype).flatten()
    return x_


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
            if type(mol.sym) is None:
                ret += [1]
            else:
                sym = numpy.array(mol.sym, dtype=numpy.int32) + 1
                if sym.ndim == 1:
                    if sym.shape[0] != mol.n_apm:
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
        ) = driver.decode_attributes(molecular_sequence)
        self.header = driver.decode_header(molecular_sequence, self.n_header)
        del driver

    ##
    # @brief general runner
    # @details mobbRMSD を計算する。
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
    #   remove_com (bool): 参照構造と対照構造から重心を除去する。 default True
    #   sort_by_g (bool): 参照構造を自己分散の大きい順に並び替えて計算を実行する。 default True
    #   rotate_y (bool): 対象構造に対して置換と回転を実行する。 default False
    # @return mobbrmsd_result
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
    ) -> mobbrmsd_result:

        dt = select_dtype(x, y)
        x_ = varidation_coordinates_1(x, self.d, self.natom, dtype=dt)
        y_ = varidation_coordinates_1(y, self.d, self.natom, dtype=dt)
        driver = select_driver(self.d, dtype=dt)
        w = numpy.empty(self.memsize, dtype=dt)

        iret, rret = driver.run(
            self.n_int,
            self.n_float,
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
        )

        ret = mobbrmsd_result(driver, self.d, self.header, iret, rret, w=w)
        del driver
        return ret

    """


    def batch_run(
        self,
        x: numpy.ndarray,
        y: None | numpy.ndarray = None,
        cutoff: float = float("inf"),
        difflim: float = 0.0,
        maxeval: int = -1,
        remove_com: bool = True,
        sort_by_g: bool = True,
        verbose: bool = False,
        n_chunk: int = 1,
        full_info: bool = False,
    ) -> list:

        if hasattr(self, "ww"):
            if self.ww.shape[1] != self.memsize or self.ww.shape[0] != self.njob:
                self.ww = numpy.empty((self.njob, self.memsize), dtype=self.dtype).T
        else:
            self.ww = numpy.empty((self.njob, self.memsize), dtype=self.dtype).T

        x_ = self.varidation_coordinates_2(x)
        if y is None:

            def res(i, j, hret, iret, rret):
                if i > j:
                    k = (i - 1) * (i - 2) // 2 + i + j - 1
                    if full_info:
                        return mobbrmsd_result(self.driver, hret, iret[k], rret[k])
                    else:
                        return self.driver.rmsd(rret[k])
                elif i < j:
                    k = (j - 1) * (j - 2) // 2 + j + i - 1
                    if full_info:
                        return mobbrmsd_result(self.driver, hret, iret[k], rret[k])
                    else:
                        return self.driver.rmsd(rret[k])
                else:
                    if full_info:
                        return mobbrmsd_result(self.driver, hret, None, None)
                    else:
                        return 0.0

            n_target = x_.shape[2]
            n_tri = (n_target * (n_target - 1)) // 2
            n_chunk_ = n_tri if n_chunk < 1 else self.njob * n_chunk
            if n_tri == n_chunk_:
                hret, iret, rret = self.driver.batch_run_tri(
                    n_target,
                    self.n_header,
                    self.n_int,
                    self.n_float,
                    n_tri,
                    1,
                    x_,
                    self.ww,
                    cutoff,
                    difflim,
                    maxeval,
                    remove_com,
                    sort_by_g,
                )
                return [
                    [res(i, j, hret, iret.T, rret.T) for j in range(n_target)]
                    for i in range(n_target)
                ]
            else:
                n_lower = 1
                nrep = (n_tri + n_chunk_ - 1) // n_chunk_
                hret = numpy.empty([self.n_header])
                iret = numpy.empty([n_tri, self.n_int])
                rret = numpy.empty([n_tri, self.n_float])
                for i in tqdm(range(nrep)):
                    l = n_lower - 1
                    u = min([l + n_chunk_, n_tri])
                    hret, iret_, rret_ = self.driver.batch_run_tri(
                        n_target,
                        self.n_header,
                        self.n_int,
                        self.n_float,
                        min(n_chunk_, n_tri - n_lower + 1),
                        n_lower,
                        x_,
                        self.ww,
                        cutoff,
                        difflim,
                        maxeval,
                        remove_com,
                        sort_by_g,
                    )
                    iret[l:u] = iret_.T
                    rret[l:u] = rret_.T
                    n_lower += n_chunk_

                return [
                    [res(i, j, hret, iret, rret) for j in range(n_target)]
                    for i in range(n_target)
                ]

        else:

            def res(n_reference, i, j, hret, iret, rret):
                k = j * n_reference + i
                if full_info:
                    return mobbrmsd_result(self.driver, hret, iret[k], rret[k])
                else:
                    return self.driver.rmsd(rret[k])

            y_ = self.varidation_coordinates_2(y)

            n_reference = x_.shape[2]
            n_target = y_.shape[2]
            n_tri = n_reference * n_target
            n_chunk_ = n_tri if n_chunk < 1 else self.njob * n_chunk
            if n_tri == n_chunk_:
                hret, iret, rret = self.driver.batch_run(
                    n_reference,
                    n_target,
                    self.n_header,
                    self.n_int,
                    self.n_float,
                    n_tri,
                    1,
                    x_,
                    y_,
                    self.ww,
                    cutoff,
                    difflim,
                    maxeval,
                    remove_com,
                    sort_by_g,
                )
                return [
                    [
                        res(n_reference, i, j, hret, iret.T, rret.T)
                        for j in range(n_target)
                    ]
                    for i in range(n_reference)
                ]
            else:
                n_lower = 1
                nrep = (n_tri + n_chunk_ - 1) // n_chunk_
                hret = numpy.empty([self.n_header])
                iret = numpy.empty([n_tri, self.n_int])
                rret = numpy.empty([n_tri, self.n_float])
                for i in tqdm(range(nrep)):
                    l = n_lower - 1
                    u = min([l + n_chunk_, n_tri])
                    hret, iret_, rret_ = self.driver.batch_run(
                        n_reference,
                        n_target,
                        self.n_header,
                        self.n_int,
                        self.n_float,
                        min(n_chunk_, n_tri - n_lower + 1),
                        n_lower,
                        x_,
                        y_,
                        self.ww,
                        cutoff,
                        difflim,
                        maxeval,
                        remove_com,
                        sort_by_g,
                    )
                    iret[l:u] = iret_.T
                    rret[l:u] = rret_.T
                    n_lower += n_chunk_

                return [
                    [res(n_reference, i, j, hret, iret, rret) for j in range(n_target)]
                    for i in range(n_reference)
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
                self.ww = numpy.empty((self.njob, self.memsize), dtype=self.dtype).T
        else:
            self.ww = numpy.empty((self.njob, self.memsize), dtype=self.dtype).T

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


    def clear(self) -> None:
        self.d, self.natom = 0, 0
        self.memsize, self.njob = 0, 0
        self.n_header, self.n_int, self.n_float = 0, 0, 0
    """

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
        if hasattr(self, "ww"):
            del self.ww

    def __str__(self):
        kws = [
            f"molecules={self.molecules}",
            f"d={self.d}",
        ]
        return "{}({})".format(type(self).__name__, ", ".join(kws))

    def __repr__(self):
        kws = [f"{key}={value!r}" for key, value in self.__dict__.items()]
        return "{}({})".format(type(self).__name__, ", ".join(kws))
