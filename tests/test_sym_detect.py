"""Tests for mobbrmsd.sym_detect — VF2 automorphisms + local distance filter.

Ported and refactored from ~/code/mobbrmsd_testcodes/ (standalone scripts) into
the pytest suite.  Tests in TestGraphAutomorphisms and TestLocalMaxDistErr use
only NetworkX + NumPy and require no MDTraj.  TestDetectSym uses MDTraj and
is skipped automatically when it is unavailable.
"""
import numpy
import pytest
import networkx

from mobbrmsd.sym_detect import graph_automorphisms, local_max_dist_err, max_dist_err


# ---------------------------------------------------------------------------
# Synthetic graph helpers (no MDTraj)
# ---------------------------------------------------------------------------

def _make_triangle():
    """All-H equilateral triangle.  Automorphism group = S3 (order 6)."""
    G = networkx.Graph()
    for i in range(3):
        G.add_node(i, element="H")
    G.add_edge(0, 1)
    G.add_edge(1, 2)
    G.add_edge(0, 2)
    return G


def _equilateral_xyz(a=1.5):
    """Perfect equilateral triangle in the xy-plane."""
    h = a * numpy.sqrt(3) / 2
    return numpy.array([[0.0, 0.0, 0.0], [a, 0.0, 0.0], [a / 2, h, 0.0]])


def _make_star3():
    """CH3-like star: C at index 0, 3 H at indices 1–3.  No H-H edges."""
    G = networkx.Graph()
    G.add_node(0, element="C")
    for i in range(1, 4):
        G.add_node(i, element="H")
        G.add_edge(0, i)
    return G


def _star3_xyz(r=1.09):
    """C at origin, 3 H at 120° intervals in the xy-plane (Å by default)."""
    angles = numpy.array([0.0, 2 * numpy.pi / 3, 4 * numpy.pi / 3])
    xyz = numpy.zeros((4, 3))
    xyz[1:, 0] = r * numpy.cos(angles)
    xyz[1:, 1] = r * numpy.sin(angles)
    return xyz


# ---------------------------------------------------------------------------
# TestGraphAutomorphisms
# ---------------------------------------------------------------------------

class TestGraphAutomorphisms:
    """VF2 automorphism enumeration on synthetic graphs — no MDTraj required."""

    def test_triangle_count(self):
        """S3 has order 6; minus identity → 5 non-identity automorphisms."""
        autos = graph_automorphisms(_make_triangle())
        assert len(autos) == 5

    def test_star3_count(self):
        """3! permutations of 3 H atoms; minus identity → 5."""
        autos = graph_automorphisms(_make_star3())
        assert len(autos) == 5

    def test_identity_excluded(self):
        G = _make_triangle()
        identity = list(range(G.number_of_nodes()))
        assert identity not in graph_automorphisms(G)

    def test_asymmetric_path_no_autos(self):
        """C-N-O path: all elements distinct, no non-trivial automorphisms."""
        G = networkx.Graph()
        G.add_node(0, element="C")
        G.add_node(1, element="N")
        G.add_node(2, element="O")
        G.add_edge(0, 1)
        G.add_edge(1, 2)
        assert graph_automorphisms(G) == []

    def test_perm_is_list_of_ints(self):
        for perm in graph_automorphisms(_make_triangle()):
            assert isinstance(perm, list)
            assert all(isinstance(v, int) for v in perm)

    def test_perm_length_equals_node_count(self):
        G = _make_star3()
        for perm in graph_automorphisms(G):
            assert len(perm) == G.number_of_nodes()


# ---------------------------------------------------------------------------
# TestLocalMaxDistErr
# ---------------------------------------------------------------------------

class TestLocalMaxDistErr:
    """local_max_dist_err returns a dimensionless fraction of mean bond length.

    Key property: the result is unit-free, so the same numeric value is
    returned whether coordinates are in Å or nm.
    """

    def test_identity_is_zero(self):
        G = _make_triangle()
        xyz = _equilateral_xyz()
        assert local_max_dist_err(xyz, [0, 1, 2], G) == pytest.approx(0.0, abs=1e-14)

    def test_perfect_triangle_c3_is_zero(self):
        G = _make_triangle()
        xyz = _equilateral_xyz()
        assert local_max_dist_err(xyz, [1, 2, 0], G) == pytest.approx(0.0, abs=1e-14)

    def test_imperfect_triangle_c3_is_nonzero(self):
        G = _make_triangle()
        xyz = _equilateral_xyz()
        xyz[2, 0] += 0.05         # break equilateral symmetry
        assert local_max_dist_err(xyz, [1, 2, 0], G) > 0.0

    def test_perfect_methyl_c3_is_zero(self):
        G = _make_star3()
        xyz = _star3_xyz()
        # C3 rotation: C stays, H0→H1→H2→H0
        assert local_max_dist_err(xyz, [0, 2, 3, 1], G) == pytest.approx(0.0, abs=1e-14)

    def test_imperfect_methyl_c3_is_nonzero(self):
        G = _make_star3()
        xyz = _star3_xyz()
        xyz[1, 0] += 0.02         # perturb one C-H bond
        assert local_max_dist_err(xyz, [0, 2, 3, 1], G) > 0.0

    def test_unit_free_angstrom_vs_nm(self):
        """Same geometry in Å and nm yields identical dimensionless error."""
        G = _make_star3()
        xyz_ang = _star3_xyz(r=1.09)
        xyz_ang[1, 0] += 0.02     # introduce imperfection
        xyz_nm = xyz_ang * 0.1    # convert to nm
        c3 = [0, 2, 3, 1]
        err_ang = local_max_dist_err(xyz_ang, c3, G)
        err_nm  = local_max_dist_err(xyz_nm,  c3, G)
        assert err_ang == pytest.approx(err_nm, rel=1e-10)

    def test_c2_swap_zero_on_symmetric(self):
        """H0↔H1 transposition is valid for a perfect equilateral triangle."""
        G = _make_triangle()
        xyz = _equilateral_xyz()
        assert local_max_dist_err(xyz, [1, 0, 2], G) == pytest.approx(0.0, abs=1e-14)


# ---------------------------------------------------------------------------
# TestMaxDistErr  (documents the global-check limitation)
# ---------------------------------------------------------------------------

class TestMaxDistErr:
    """max_dist_err checks the *global* distance matrix.

    It correctly identifies whole-molecule symmetries (e.g. water C2, benzene
    D6h) but incorrectly rejects valid *local* symmetries such as methyl C3
    because the distant parts of the molecule break the global symmetry.
    This class documents that known limitation.
    """

    def test_perfect_triangle_c3_global_zero(self):
        """For a pure-C3 triangle (no other atoms), global check is exact."""
        xyz = _equilateral_xyz()
        assert max_dist_err(xyz, [1, 2, 0]) == pytest.approx(0.0, abs=1e-14)

    def test_local_symmetry_fails_global_check(self):
        """Methyl C3 with a remote atom: global check gives large error."""
        xyz = _star3_xyz()
        # Add a remote atom that breaks C3 globally
        xyz_ext = numpy.vstack([xyz, [[5.0, 0.0, 0.0]]])   # shape (5, 3)
        # C3 perm on star; remote atom stays
        c3_ext = [0, 2, 3, 1, 4]
        err_global = max_dist_err(xyz_ext, c3_ext)
        err_local  = local_max_dist_err(xyz, [0, 2, 3, 1], _make_star3())
        # Global check incorrectly rejects; local check correctly accepts
        assert err_global > 0.1
        assert err_local  < 0.01


# ---------------------------------------------------------------------------
# TestDetectSym  (requires MDTraj)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def water_top_xyz():
    """Single-residue water topology with O-H1-H2 bonds, coordinates in Å."""
    mdtraj = pytest.importorskip("mdtraj")
    top = mdtraj.Topology()
    chain = top.add_chain()
    res = top.add_residue("WAT", chain)
    O  = top.add_atom("O",  mdtraj.element.oxygen,   res)
    H1 = top.add_atom("H1", mdtraj.element.hydrogen, res)
    H2 = top.add_atom("H2", mdtraj.element.hydrogen, res)
    top.add_bond(O, H1)
    top.add_bond(O, H2)
    # O-H bond = 0.96 Å, H-O-H angle = 104.5°
    angle = 104.5 * numpy.pi / 180
    xyz = numpy.array([
        [0.000, 0.000, 0.000],
        [0.960, 0.000, 0.000],
        [0.960 * numpy.cos(angle), 0.960 * numpy.sin(angle), 0.000],
    ])
    res_indices = [a.index for a in res.atoms]
    return top, xyz, res_indices


class TestDetectSym:
    """Full detect_sym pipeline — requires MDTraj (auto-skipped if absent)."""

    def test_no_template_returns_all_vf2_autos(self, water_top_xyz):
        from mobbrmsd.sym_detect import detect_sym
        top, xyz, res_indices = water_top_xyz
        sym, report = detect_sym(top, res_indices, template_xyz=None)
        # Water C2: exactly 1 non-identity automorphism (H-swap)
        assert len(sym) == 1
        # Pure-VF2 mode: local_dist_err is set to 0.0 for all candidates
        assert all(r["local_dist_err"] == 0.0 for r in report)

    def test_water_c2_automorphism_is_hswap(self, water_top_xyz):
        from mobbrmsd.sym_detect import detect_sym
        top, xyz, res_indices = water_top_xyz
        sym, _ = detect_sym(top, res_indices, template_xyz=None)
        assert sym[0] == [0, 2, 1]   # O stays, H1 ↔ H2

    def test_water_c2_detected_with_template(self, water_top_xyz):
        from mobbrmsd.sym_detect import detect_sym
        top, xyz, res_indices = water_top_xyz
        sym, report = detect_sym(top, res_indices, template_xyz=xyz, eps=0.05)
        assert len(sym) == 1
        assert report[0]["local_dist_err"] < 0.05

    def test_distorted_geometry_rejected(self, water_top_xyz):
        """Moving H1 far off makes the H-swap geometrically invalid."""
        from mobbrmsd.sym_detect import detect_sym
        top, xyz, res_indices = water_top_xyz
        xyz_bad = xyz.copy()
        xyz_bad[1, 0] += 0.30        # shift H1 by 0.3 Å → local_dist_err ≫ 0.05
        sym, _ = detect_sym(top, res_indices, template_xyz=xyz_bad, eps=0.05)
        assert len(sym) == 0

    def test_report_sorted_ascending_by_dist_err(self, water_top_xyz):
        from mobbrmsd.sym_detect import detect_sym
        top, xyz, res_indices = water_top_xyz
        _, report = detect_sym(top, res_indices, template_xyz=xyz, eps=1.0)
        errs = [r["local_dist_err"] for r in report]
        assert errs == sorted(errs)

    def test_report_contains_all_vf2_candidates(self, water_top_xyz):
        """Number of report entries equals number of VF2 automorphisms."""
        from mobbrmsd.sym_detect import detect_sym
        top, xyz, res_indices = water_top_xyz
        _, report_no_filter   = detect_sym(top, res_indices, template_xyz=None)
        _, report_with_filter = detect_sym(top, res_indices, template_xyz=xyz, eps=1.0)
        assert len(report_no_filter) == len(report_with_filter)

    def test_unit_free_eps_angstrom_vs_nm(self, water_top_xyz):
        """Same eps accepts the same permutations regardless of coordinate unit."""
        from mobbrmsd.sym_detect import detect_sym
        top, xyz_ang, res_indices = water_top_xyz
        xyz_nm = xyz_ang * 0.1
        sym_ang, _ = detect_sym(top, res_indices, template_xyz=xyz_ang, eps=0.05)
        sym_nm,  _ = detect_sym(top, res_indices, template_xyz=xyz_nm,  eps=0.05)
        assert sym_ang == sym_nm
