"""Intramolecular symmetry detection via VF2 graph automorphisms + local distance filter.

A permutation σ is accepted as a valid molecular symmetry when it preserves
pairwise interatomic distances within the *local neighbourhood* of the permuted
atoms (1-hop by default)::

    d(i, j) ≈ d(σ(i), σ(j))   for pairs (i, j) near the permuted atoms

This correctly detects *local* symmetry (e.g. methyl C3 rotation) as well as
*global* symmetry (e.g. water C2), while naturally rejecting pseudo-symmetric
operations that break local geometry.

Typical usage::

    from mobbrmsd.sym_detect import detect_sym

    sym, report = detect_sym(top, res_indices, template_xyz, eps=0.05)
    # sym    : list of permutation lists accepted by molecules(sym=...)
    # report : list of {"perm": [...], "local_dist_err": float} for ALL
    #          non-identity graph-automorphic candidates, sorted by local_dist_err

Public API
----------
residue_graph(top, res_indices)                    -> nx.Graph  (element-labelled)
graph_automorphisms(G)                             -> list[list[int]]  (non-identity)
max_dist_err(xyz, perm)                            -> float  (global check)
local_max_dist_err(xyz, perm, G, hop=1)            -> float  (local neighbourhood check)
detect_sym(top, res_indices, template_xyz, eps)    -> (sym, report)
"""
from __future__ import annotations

import numpy
import numpy.typing as npt
from typing import List, Dict, Any, Tuple, Optional


# ---------------------------------------------------------------------------
# Step 1 – build residue graph
# ---------------------------------------------------------------------------

def residue_graph(top: Any, res_indices: List[int]):
    """Build an element-labelled NetworkX graph for *res_indices*.

    Parameters
    ----------
    top         : ``mdtraj.Topology``
    res_indices : global atom indices that belong to one residue (in local order)

    Returns
    -------
    G : ``networkx.Graph``
        Nodes are 0-based local indices with attribute ``element`` (str).
        Edges are bonds between atoms that are both in *res_indices*.
    """
    import networkx as nx

    local = {gi: li for li, gi in enumerate(res_indices)}

    G = nx.Graph()
    for li, gi in enumerate(res_indices):
        G.add_node(li, element=top.atom(gi).element.symbol)
    for bond in top.bonds:
        ai, bi = bond.atom1.index, bond.atom2.index
        if ai in local and bi in local:
            G.add_edge(local[ai], local[bi])
    return G


# ---------------------------------------------------------------------------
# Step 2 – enumerate graph automorphisms
# ---------------------------------------------------------------------------

def graph_automorphisms(G) -> List[List[int]]:
    """Return all non-identity automorphisms of *G* as permutation lists.

    Uses the VF2 algorithm (``networkx.GraphMatcher``).  The identity
    permutation ``[0, 1, 2, ...]`` is excluded.

    Each entry ``perm`` satisfies: ``perm[i]`` is the node that node *i* maps
    to under the automorphism, i.e. the permuted coordinate array is
    ``X_perm = X[perm]``.
    """
    from networkx.algorithms.isomorphism import GraphMatcher, categorical_node_match

    n = G.number_of_nodes()
    identity = list(range(n))
    node_match = categorical_node_match("element", "X")
    gm = GraphMatcher(G, G, node_match=node_match)
    perms = []
    for mapping in gm.isomorphisms_iter():
        perm = [mapping[i] for i in range(n)]
        if perm != identity:
            perms.append(perm)
    return perms


# ---------------------------------------------------------------------------
# Step 3 – distance-matrix preservation checks
# ---------------------------------------------------------------------------

def max_dist_err(xyz: npt.NDArray, perm: List[int]) -> float:
    """Maximum absolute change in any pairwise distance under *perm*.

    Checks the **global** distance matrix.  Returns ≈ 0 only for true
    whole-molecule symmetries (e.g. C2v of water, D6h of benzene).  For
    *local* symmetries such as methyl C3 rotation the value is large because
    distant atoms break the symmetry of the rest of the molecule.  Use
    :func:`local_max_dist_err` for local symmetry detection.

    Parameters
    ----------
    xyz  : ``(n, 3)`` coordinate array (any consistent unit)
    perm : permutation as a list of ints of length *n*

    Returns
    -------
    float : max |d(i,j) - d(σ(i), σ(j))| over all atom pairs (i, j)
    """
    xyz = numpy.asarray(xyz, dtype=float)
    p = numpy.asarray(perm, dtype=numpy.intp)
    diff = xyz[:, None, :] - xyz[None, :, :]        # (n, n, 3)
    D = numpy.sqrt((diff ** 2).sum(-1))              # (n, n)
    Dp = D[numpy.ix_(p, p)]                          # permuted distance matrix
    return float(numpy.abs(D - Dp).max())


def local_max_dist_err(xyz: npt.NDArray, perm: List[int], G, hop: int = 1) -> float:
    """Max pairwise-distance change within the local neighbourhood, as a fraction of mean bond length.

    Only atom pairs where **at least one atom** lies within *hop* bonds of a
    permuted atom are checked.  The raw max error is divided by the mean bond
    length of edges fully within the neighbourhood, making the result
    **unit-free** (same threshold works for coordinates in Å or nm).

    Typical values with hop=1 on MM-optimised geometries:

    * ``≈ 0``    for exact symmetry (water H-swap, ideal methyl C3)
    * ``≈ 0.01`` for slightly imperfect MM structures (MMFF94 methyl C3)
    * ``> 0.1``  for diastereotopic atoms (NH₂ in a chiral amino acid)

    Parameters
    ----------
    xyz  : ``(n, 3)`` coordinate array (any consistent unit)
    perm : permutation as a list of ints of length *n*
    G    : NetworkX graph (from :func:`residue_graph`) providing bond topology
    hop  : number of bond hops that define the local neighbourhood (default 1)

    Returns
    -------
    float : max |d(i,j) - d(σ(i), σ(j))| / mean_bond_length, dimensionless
    """
    n = len(perm)
    perm_arr = numpy.asarray(perm, dtype=numpy.intp)
    xyz = numpy.asarray(xyz, dtype=float)

    moved = {i for i in range(n) if perm[i] != i}
    if not moved:
        return 0.0
    neighborhood = set(moved)
    for _ in range(hop):
        ext: set = set()
        for a in neighborhood:
            ext.update(G.neighbors(a))
        neighborhood.update(ext)
    nb = sorted(neighborhood)
    nb_set = set(nb)

    # Normalise by mean bond length of edges fully within the neighbourhood.
    bond_lens = [
        float(numpy.linalg.norm(xyz[a] - xyz[b]))
        for (a, b) in G.edges()
        if a in nb_set and b in nb_set
    ]
    ref = float(numpy.mean(bond_lens)) if bond_lens else 1.0

    diff = xyz[:, None, :] - xyz[None, :, :]
    D  = numpy.sqrt((diff ** 2).sum(-1))
    Dp = D[numpy.ix_(perm_arr, perm_arr)]
    nb_arr = numpy.asarray(nb, dtype=numpy.intp)
    return float(numpy.abs(D - Dp)[numpy.ix_(nb_arr, nb_arr)].max()) / ref


# ---------------------------------------------------------------------------
# Step 4 – full pipeline
# ---------------------------------------------------------------------------

def detect_sym(
    top: Any,
    res_indices: List[int],
    template_xyz: Optional[npt.NDArray] = None,
    eps: float = 0.05,
    hop: int = 1,
) -> Tuple[List[List[int]], List[Dict[str, Any]]]:
    """Detect intramolecular symmetry permutations for one residue.

    Algorithm
    ---------
    1. Build element-labelled bond graph from *top* / *res_indices*.
    2. Enumerate **all** graph automorphisms via VF2 (non-identity only).
    3. If *template_xyz* is provided, filter candidates by
       ``local_max_dist_err(xyz, σ, G, hop=hop) < eps`` — checking only
       pairwise distances within *hop* bonds of the permuted atoms.
       If *template_xyz* is ``None``, **all** VF2 candidates are accepted
       (pure graph-topology symmetry, no geometric filter).
    4. Return accepted permutations and a full report for all candidates.

    Why local (not global) distance check?
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A global distance-matrix check rejects valid *local* symmetries such as
    methyl C3 rotation: while the three H atoms form an equilateral triangle
    locally, their distances to distant atoms (e.g. OXT) differ because the
    whole molecule is not C3-symmetric.  Checking only the 1-hop neighbourhood
    (the permuted atoms + their direct bond partners) captures the local
    geometry without being confused by the rest of the molecule.

    Parameters
    ----------
    top          : ``mdtraj.Topology``
    res_indices  : global atom indices for one representative residue,
                   in the local order matching ``sym`` indexing
                   (``res_indices[i]`` = global index of local atom *i*)
    template_xyz : ``(n_atoms_total, 3)`` or ``(n_frames, n_atoms_total, 3)``
                   MM-optimised coordinates (any consistent unit).
                   When 3-D, the first frame is used.
                   When ``None``, all VF2 automorphisms are accepted.
    eps          : acceptance threshold for ``local_max_dist_err`` — a
                   **unit-free fraction of the mean bond length** (default 0.05,
                   i.e. 5 % of the local bond length).  The same value works
                   regardless of whether *template_xyz* is in Å or nm.
    hop          : neighbourhood radius for the local check (default 1,
                   i.e. direct bond partners of permuted atoms)

    Returns
    -------
    sym    : list of accepted non-identity permutations, each a ``list[int]``
             of length ``n_apm`` — directly usable as ``molecules(sym=sym)``
    report : list of dicts for **every** non-identity VF2 candidate, sorted by
             ``local_dist_err`` (or 0 when no template is given).  Each dict:
             ``{"perm": list[int], "local_dist_err": float}`` where
             ``local_dist_err`` is dimensionless (fraction of bond length).
             Entries with ``local_dist_err < eps`` (or all when no template)
             are the accepted ones.
    """
    G = residue_graph(top, res_indices)
    candidates = graph_automorphisms(G)

    if template_xyz is None:
        # Pure VF2: accept all topological automorphisms
        report: List[Dict[str, Any]] = [
            {"perm": perm, "local_dist_err": 0.0} for perm in candidates
        ]
        return candidates, report

    xyz = numpy.asarray(template_xyz, dtype=float)
    if xyz.ndim == 3:
        xyz = xyz[0]
    xyz_res = xyz[res_indices]  # (n_apm, 3)

    report = []
    for perm in candidates:
        err = local_max_dist_err(xyz_res, perm, G, hop=hop)
        report.append({"perm": perm, "local_dist_err": err})

    report.sort(key=lambda d: d["local_dist_err"])
    sym = [d["perm"] for d in report if d["local_dist_err"] < eps]
    return sym, report
