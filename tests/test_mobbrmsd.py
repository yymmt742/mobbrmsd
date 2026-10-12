import numpy
import pytest
import networkx
from mobbrmsd.dataclass import molecules, molecular_system, _sort_species_by_variance
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


class TestAtomOrder:
    """atom_order field and sorted_by_variance on molecular_system."""

    def _make_two_species_coords(self, seed=0):
        """Two species: 3-atom x 4-mol (species A) and 2-atom x 3-mol (species B).
        Returns coords [natom, 3] in A-then-B order and the split molecular_system."""
        rng = numpy.random.default_rng(seed)
        n_apm_a, n_mol_a = 3, 4
        n_apm_b, n_mol_b = 2, 3
        natom = n_apm_a * n_mol_a + n_apm_b * n_mol_b
        coords = rng.random((natom, 3))
        mols_a = molecules(n_apm=n_apm_a, n_mol=n_mol_a, name="A")
        mols_b = molecules(n_apm=n_apm_b, n_mol=n_mol_b, name="B")
        ms = molecular_system([mols_a, mols_b])
        return coords, ms, n_apm_a, n_mol_a, n_apm_b, n_mol_b

    def test_no_atom_order_runs(self):
        """molecular_system without atom_order works as before."""
        cog = make_cogen()
        x, y = cog.generate_pair([3, 2], [4, 3], alpha=0.5, beta=1.0, remove_com=False)
        ms = molecular_system([molecules(n_apm=3, n_mol=4), molecules(n_apm=2, n_mol=3)])
        ret = mobbrmsd(mols=ms).run(x, y)
        assert ret.is_finished()

    def test_identity_atom_order_same_result(self):
        """atom_order = identity permutation gives same RMSD as no atom_order."""
        cog = make_cogen()
        x, y = cog.generate_pair(3, 4, alpha=0.8, beta=1.0)
        mols = molecules(n_apm=3, n_mol=4)
        natom = 3 * 4
        identity = numpy.arange(natom)
        ms_plain = molecular_system([mols])
        ms_ordered = molecular_system([mols], atom_order=identity)
        ret_plain = mobbrmsd(mols=ms_plain).run(x, y)
        ret_ordered = mobbrmsd(mols=ms_ordered).run(x, y)
        assert ret_plain.rmsd() == pytest.approx(ret_ordered.rmsd(), rel=1e-9)

    def test_atom_order_permuted_coords_same_result(self):
        """Permuting coords with atom_order inverse gives same RMSD as plain."""
        rng = numpy.random.default_rng(7)
        n_apm, n_mol = 3, 4
        natom = n_apm * n_mol
        x = rng.random((natom, 3))
        y = rng.random((natom, 3))

        mols = molecules(n_apm=n_apm, n_mol=n_mol)
        rmsd_plain = mobbrmsd(mols=mols).run(x, y).rmsd()

        # Shuffle atom order and tell molecular_system about it
        perm = rng.permutation(natom)
        inv_perm = numpy.argsort(perm)
        ms = molecular_system([mols], atom_order=perm)
        # Supply coords in shuffled order; atom_order undoes the shuffle
        rmsd_ordered = mobbrmsd(mols=ms).run(x[inv_perm], y[inv_perm]).rmsd()
        assert rmsd_plain == pytest.approx(rmsd_ordered, rel=1e-6)

    def test_sorted_by_variance_is_valid_permutation(self):
        """sorted_by_variance returns valid atom_order covering all atoms."""
        coords, ms, n_apm_a, n_mol_a, n_apm_b, n_mol_b = self._make_two_species_coords()
        ms2 = ms.sorted_by_variance(coords)
        natom = n_apm_a * n_mol_a + n_apm_b * n_mol_b
        assert ms2.atom_order is not None
        assert sorted(ms2.atom_order.tolist()) == list(range(natom))

    def test_sorted_by_variance_same_rmsd(self):
        """sorted_by_variance does not change the RMSD result."""
        rng = numpy.random.default_rng(13)
        n_apm_a, n_mol_a = 3, 4
        n_apm_b, n_mol_b = 2, 3
        natom = n_apm_a * n_mol_a + n_apm_b * n_mol_b
        x = rng.random((natom, 3))
        y = rng.random((natom, 3))

        ms = molecular_system([
            molecules(n_apm=n_apm_a, n_mol=n_mol_a),
            molecules(n_apm=n_apm_b, n_mol=n_mol_b),
        ])
        ms_sorted = ms.sorted_by_variance(x)

        rmsd_plain = mobbrmsd(mols=ms).run(x, y).rmsd()
        rmsd_sorted = mobbrmsd(mols=ms_sorted).run(x[ms_sorted.atom_order], y[ms_sorted.atom_order]).rmsd()
        assert rmsd_plain == pytest.approx(rmsd_sorted, rel=1e-6)

    def test_sort_species_by_variance_ordering(self):
        """_sort_species_by_variance puts higher-variance species first."""
        rng = numpy.random.default_rng(99)
        # Species A: very spread out; species B: very compact
        n_apm, n_mol = 1, 6
        coords_a = rng.random((n_mol, 3)) * 100.0   # large spread
        coords_b = rng.random((n_mol, 3)) * 0.01    # tiny spread
        coords = numpy.vstack([coords_a, coords_b])

        mols_a = molecules(n_apm=n_apm, n_mol=n_mol, name="A")
        mols_b = molecules(n_apm=n_apm, n_mol=n_mol, name="B")
        atom_order = numpy.arange(2 * n_mol)

        new_mols, new_order = _sort_species_by_variance([mols_a, mols_b], atom_order, coords)
        # Species A (high variance) should come first
        assert new_mols[0].name == "A"
        assert new_mols[1].name == "B"


class TestFromMDTraj:
    """molecular_system.from_mdtraj topology auto-detection."""

    @pytest.fixture
    def simple_topology(self):
        """3 water (3 atoms each) + 2 methanol (6 atoms each)."""
        mdtraj = pytest.importorskip("mdtraj")
        top = mdtraj.Topology()
        chain = top.add_chain()
        for _ in range(3):
            res = top.add_residue("WAT", chain)
            for elem_sym in ("O", "H", "H"):
                top.add_atom(elem_sym, mdtraj.element.Element.getBySymbol(elem_sym), res)
        for _ in range(2):
            res = top.add_residue("MOL", chain)
            for elem_sym in ("C", "O", "H", "H", "H", "H"):
                top.add_atom(elem_sym, mdtraj.element.Element.getBySymbol(elem_sym), res)
        return top

    def test_species_count(self, simple_topology):
        ms = molecular_system.from_mdtraj(simple_topology)
        assert len(ms.mols) == 2

    def test_n_mol_correct(self, simple_topology):
        ms = molecular_system.from_mdtraj(simple_topology)
        by_name = {m.name: m for m in ms.mols}
        assert by_name["WAT"].n_mol == 3
        assert by_name["MOL"].n_mol == 2

    def test_n_apm_correct(self, simple_topology):
        ms = molecular_system.from_mdtraj(simple_topology)
        by_name = {m.name: m for m in ms.mols}
        assert by_name["WAT"].n_apm == 3
        assert by_name["MOL"].n_apm == 6

    def test_atom_order_is_valid_permutation(self, simple_topology):
        ms = molecular_system.from_mdtraj(simple_topology)
        natom = 3 * 3 + 2 * 6
        assert ms.atom_order is not None
        assert sorted(ms.atom_order.tolist()) == list(range(natom))

    def test_rmsd_with_atom_order_zero_on_self(self, simple_topology):
        ms = molecular_system.from_mdtraj(simple_topology)
        natom = 3 * 3 + 2 * 6
        rng = numpy.random.default_rng(0)
        coords = rng.random((natom, 3))
        # Pass coords in original (topology) order; atom_order remaps them
        ret = mobbrmsd(mols=ms).run(coords, coords)
        assert ret.rmsd() == pytest.approx(0.0, abs=1e-6)

    def test_sort_by_variance_with_coords(self, simple_topology):
        mdtraj = pytest.importorskip("mdtraj")
        natom = 3 * 3 + 2 * 6
        rng = numpy.random.default_rng(5)
        coords = rng.random((natom, 3))
        ms = molecular_system.from_mdtraj(simple_topology, coords=coords)
        # Species order may differ but atom_order must still be a valid permutation
        assert sorted(ms.atom_order.tolist()) == list(range(natom))


class TestBuildJsonRoundtrip:
    """TopologyBuilder → build() → molecular_system.to_dict() → parse.load() round-trip.

    Verifies that the JSON produced by ``command_build`` (via ``mols.to_dict()``)
    can be consumed by ``run -i`` (via ``parse.load()``).
    """

    def test_to_dict_has_molecular_system_key(self):
        from mobbrmsd.run.topology_builder import TopologyBuilder

        mols = TopologyBuilder().add(n_apm=2, n_mol=3).build()
        d = mols.to_dict()
        assert "molecular_system" in d
        assert "molecules" in d["molecular_system"]

    def test_simple_roundtrip(self):
        from mobbrmsd.run.topology_builder import TopologyBuilder
        from mobbrmsd.run.parse import load

        mols = TopologyBuilder().add(n_apm=3, n_mol=4, name="Water").build()
        loaded = load(mols.to_dict())
        assert loaded is not None
        assert len(loaded.mols) == 1
        assert loaded.mols[0].n_apm == 3
        assert loaded.mols[0].n_mol == 4
        assert loaded.mols[0].name == "Water"

    def test_sym_preserved(self):
        from mobbrmsd.run.topology_builder import TopologyBuilder
        from mobbrmsd.run.parse import load

        mols = TopologyBuilder().add(n_apm=3, n_mol=4, sym=[[1, 0, 2]]).build()
        loaded = load(mols.to_dict())
        assert loaded is not None
        assert loaded.mols[0].sym == [[1, 0, 2]]

    def test_multispecies_roundtrip(self):
        from mobbrmsd.run.topology_builder import TopologyBuilder
        from mobbrmsd.run.parse import load

        mols = (
            TopologyBuilder()
            .add(n_apm=3, n_mol=4, name="Water")
            .add(n_apm=6, n_mol=2, name="Methanol")
            .build()
        )
        loaded = load(mols.to_dict())
        assert loaded is not None
        assert len(loaded.mols) == 2
        by_name = {m.name: m for m in loaded.mols}
        assert by_name["Water"].n_apm == 3
        assert by_name["Water"].n_mol == 4
        assert by_name["Methanol"].n_apm == 6
        assert by_name["Methanol"].n_mol == 2

    def test_json_file_roundtrip(self, tmp_path):
        import json
        from mobbrmsd.run.topology_builder import TopologyBuilder
        from mobbrmsd.run.parse import load

        mols = TopologyBuilder().add(n_apm=3, n_mol=4, sym=[[1, 0, 2]], name="Water").build()
        out = tmp_path / "topology.json"
        out.write_text(json.dumps(mols.to_dict(), indent=2))
        loaded = load(json.loads(out.read_text()))
        assert loaded is not None
        assert loaded.mols[0].n_apm == 3
        assert loaded.mols[0].n_mol == 4
        assert loaded.mols[0].sym == [[1, 0, 2]]

    def test_atom_names_not_in_output(self):
        """atom_names are builder-only and must not appear in molecular_system JSON."""
        from mobbrmsd.run.topology_builder import TopologyBuilder

        mols = (
            TopologyBuilder()
            .add(n_apm=3, n_mol=4, atom_names=["OW", "HW1", "HW2"])
            .build()
        )
        d = mols.to_dict()
        for mol_entry in d["molecular_system"]["molecules"]:
            assert "atom_names" not in mol_entry

    def test_loaded_system_runs(self):
        """A molecular_system loaded from build JSON can run mobbrmsd without error."""
        from mobbrmsd.run.topology_builder import TopologyBuilder
        from mobbrmsd.run.parse import load

        mols = TopologyBuilder().add(n_apm=3, n_mol=4, sym=[[1, 0, 2]]).build()
        loaded = load(mols.to_dict())

        rng = numpy.random.default_rng(0)
        natom = 3 * 4
        x = rng.random((natom, 3))
        y = rng.random((natom, 3))
        ret = mobbrmsd(mols=loaded).run(x, y)
        assert ret.is_finished()
        assert ret.rmsd() >= 0.0
