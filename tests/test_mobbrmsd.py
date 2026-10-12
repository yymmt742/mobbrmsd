import numpy
import pytest
import networkx
from mobbrmsd.dataclass import molecules, molecular_system
from mobbrmsd.mobbrmsd import mobbrmsd
from mobbrmsd.demo._coord_generator import coord_generator


def make_cogen(d=3, seed=42):
    cog = coord_generator(d=d)
    cog.rng = numpy.random.default_rng(seed)
    cog.sog.rng = numpy.random.default_rng(seed)
    return cog


class TestRun3D:
    """Single-pair 3D RMSD (corresponds to demo/bb.py)."""

    def test_terminates(self):
        cog = make_cogen()
        x, y = cog.generate_pair(3, 4, alpha=0.8, beta=1.0)
        mols = molecules(n_apm=3, n_mol=4, sym=[[1, 0, 2]])
        ret = mobbrmsd(mols=mols).run(x, y)
        assert ret.is_finished()

    def test_rmsd_nonnegative(self):
        cog = make_cogen()
        x, y = cog.generate_pair(3, 4, alpha=0.8, beta=1.0)
        mols = molecules(n_apm=3, n_mol=4)
        ret = mobbrmsd(mols=mols).run(x, y)
        assert ret.rmsd() >= 0.0

    def test_same_coords_zero_rmsd(self):
        cog = make_cogen()
        x, _ = cog.generate_pair(3, 4, alpha=0.8, beta=1.0)
        mols = molecules(n_apm=3, n_mol=4)
        ret = mobbrmsd(mols=mols).run(x, x)
        assert ret.rmsd() == pytest.approx(0.0, abs=1e-6)

    def test_bounds_enclose_rmsd(self):
        cog = make_cogen()
        x, y = cog.generate_pair(3, 4, alpha=0.8, beta=1.0)
        mols = molecules(n_apm=3, n_mol=4, sym=[[1, 0, 2]])
        ret = mobbrmsd(mols=mols).run(x, y)
        assert ret.lowerbound_as_rmsd() <= ret.rmsd() + 1e-9
        assert ret.rmsd() <= ret.upperbound_as_rmsd() + 1e-9

    def test_rmsd_matches_simple_interface(self):
        cog = make_cogen()
        x, y = cog.generate_pair(3, 4, alpha=0.8, beta=1.0)
        mols = molecules(n_apm=3, n_mol=4)
        mrmsd = mobbrmsd(mols=mols)
        assert mrmsd.rmsd(x, y) == pytest.approx(mrmsd.run(x, y).rmsd(), rel=1e-6)

    def test_accepts_rank3_input(self):
        cog = make_cogen()
        x, y = cog.generate_pair(3, 4, alpha=0.8, beta=1.0)
        mols = molecules(n_apm=3, n_mol=4)
        mrmsd = mobbrmsd(mols=mols)
        ret_r2 = mrmsd.run(x, y)
        ret_r3 = mrmsd.run(x[numpy.newaxis], y[numpy.newaxis])
        assert ret_r2.rmsd() == pytest.approx(ret_r3.rmsd(), rel=1e-9)

    def test_multiframe_input_raises(self):
        cog = make_cogen()
        n_apm, n_mol = 3, 4
        xs = numpy.array([cog.generate(n_apm, n_mol, 0.5, 1.0) for _ in range(3)])
        mols = molecules(n_apm=n_apm, n_mol=n_mol)
        with pytest.raises(ValueError):
            mobbrmsd(mols=mols).run(xs, xs)


class TestRun2D:
    """Single-pair 2D RMSD (corresponds to demo/bb_2d.py)."""

    def test_terminates(self):
        cog = make_cogen(d=2)
        x, y = cog.generate_pair(3, 4, alpha=0.8, beta=1.0)
        mols = molecules(n_apm=3, n_mol=4, sym=[[1, 0, 2]])
        ret = mobbrmsd(d=2, mols=mols).run(x, y)
        assert ret.is_finished()

    def test_same_coords_zero_rmsd(self):
        cog = make_cogen(d=2)
        x, _ = cog.generate_pair(3, 4, alpha=0.8, beta=1.0)
        mols = molecules(n_apm=3, n_mol=4)
        ret = mobbrmsd(d=2, mols=mols).run(x, x)
        assert ret.rmsd() == pytest.approx(0.0, abs=1e-6)


class TestRunMultiComponent:
    """Multi-species molecular system (corresponds to demo/bb_multi.py)."""

    def test_terminates(self):
        cog = make_cogen()
        x, y = cog.generate_pair([2, 3], [3, 2], alpha=0.5, beta=1.0, remove_com=False)
        mols = molecular_system([
            molecules(n_apm=2, n_mol=3, sym=[[1, 0]]),
            molecules(n_apm=3, n_mol=2, sym=[[1, 0, 2]]),
        ])
        ret = mobbrmsd(mols=mols).run(x, y)
        assert ret.is_finished()

    def test_same_coords_zero_rmsd(self):
        cog = make_cogen()
        x, _ = cog.generate_pair([2, 3], [3, 2], alpha=0.5, beta=1.0, remove_com=False)
        mols = molecular_system([
            molecules(n_apm=2, n_mol=3),
            molecules(n_apm=3, n_mol=2),
        ])
        ret = mobbrmsd(mols=mols).run(x, x)
        assert ret.rmsd() == pytest.approx(0.0, abs=1e-6)


class TestBatchRun:
    """Batch RMSD computation (corresponds to demo/batch.py and batch_tri.py)."""

    def test_shape_ref_vs_target(self):
        cog = make_cogen()
        n_apm, n_mol, n_ref, n_tgt = 3, 4, 5, 8
        xs = numpy.array([cog.generate(n_apm, n_mol, 0.5, 1.0) for _ in range(n_ref)])
        ys = numpy.array([cog.generate(n_apm, n_mol, 0.5, 1.0) for _ in range(n_tgt)])
        mols = molecules(n_apm=n_apm, n_mol=n_mol)
        mat = mobbrmsd(mols=mols).batch_run(xs, ys)
        assert mat.shape == (n_tgt, n_ref)

    def test_ref_target_nonnegative(self):
        cog = make_cogen()
        n_apm, n_mol, n_ref, n_tgt = 3, 4, 5, 8
        xs = numpy.array([cog.generate(n_apm, n_mol, 0.5, 1.0) for _ in range(n_ref)])
        ys = numpy.array([cog.generate(n_apm, n_mol, 0.5, 1.0) for _ in range(n_tgt)])
        mols = molecules(n_apm=n_apm, n_mol=n_mol)
        mat = mobbrmsd(mols=mols).batch_run(xs, ys)
        assert numpy.all(mat >= 0)

    def test_triangular_shape(self):
        cog = make_cogen()
        n_apm, n_mol, n_str = 3, 4, 6
        xs = numpy.array([cog.generate(n_apm, n_mol, 0.5, 1.0) for _ in range(n_str)])
        mols = molecules(n_apm=n_apm, n_mol=n_mol)
        mat = mobbrmsd(mols=mols).batch_run(xs)
        assert mat.shape == (n_str, n_str)

    def test_triangular_symmetric(self):
        cog = make_cogen()
        n_apm, n_mol, n_str = 3, 4, 6
        xs = numpy.array([cog.generate(n_apm, n_mol, 0.5, 1.0) for _ in range(n_str)])
        mols = molecules(n_apm=n_apm, n_mol=n_mol)
        mat = mobbrmsd(mols=mols).batch_run(xs)
        numpy.testing.assert_allclose(mat, mat.T, atol=1e-9)

    def test_triangular_diagonal_zero(self):
        cog = make_cogen()
        n_apm, n_mol, n_str = 3, 4, 6
        xs = numpy.array([cog.generate(n_apm, n_mol, 0.5, 1.0) for _ in range(n_str)])
        mols = molecules(n_apm=n_apm, n_mol=n_mol)
        mat = mobbrmsd(mols=mols).batch_run(xs)
        numpy.testing.assert_allclose(numpy.diag(mat), 0.0, atol=1e-6)

    def test_triangular_nonnegative(self):
        cog = make_cogen()
        n_apm, n_mol, n_str = 3, 4, 6
        xs = numpy.array([cog.generate(n_apm, n_mol, 0.5, 1.0) for _ in range(n_str)])
        mols = molecules(n_apm=n_apm, n_mol=n_mol)
        mat = mobbrmsd(mols=mols).batch_run(xs)
        assert numpy.all(mat >= 0)


class TestMinSpanTree:
    """Minimum spanning tree (corresponds to demo/mst.py)."""

    def test_is_spanning_tree(self):
        cog = make_cogen()
        n_apm, n_mol, n_str = 3, 4, 8
        xs = numpy.array([cog.generate(n_apm, n_mol, 0.5, 1.0) for _ in range(n_str)])
        mols = molecules(n_apm=n_apm, n_mol=n_mol)
        g = mobbrmsd(mols=mols).min_span_tree(xs)
        assert isinstance(g, networkx.Graph)
        assert networkx.is_connected(g)
        assert g.number_of_nodes() == n_str
        assert g.number_of_edges() == n_str - 1

    def test_edge_weights_nonnegative(self):
        cog = make_cogen()
        n_apm, n_mol, n_str = 3, 4, 6
        xs = numpy.array([cog.generate(n_apm, n_mol, 0.5, 1.0) for _ in range(n_str)])
        mols = molecules(n_apm=n_apm, n_mol=n_mol)
        g = mobbrmsd(mols=mols).min_span_tree(xs)
        weights = networkx.get_edge_attributes(g, "weight").values()
        assert all(w >= 0 for w in weights)
