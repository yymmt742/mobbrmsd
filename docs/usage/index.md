# Usage

`mobbrmsd` represents a molecular assembly as a collection of molecular
species.

Each molecular species is described by:

- the number of atoms per molecule,
- the number of molecules,
- optional intramolecular permutations,
- an optional species name.

## Molecular system

A molecular system can contain multiple molecular species.

For example:

```python
from mobbrmsd.dataclass import molecules, molecular_system

water = molecules(
    n_apm=3,
    n_mol=4,
    sym=[[0, 2, 1]],
    name="Water",
)

system = molecular_system(
    mols=[water],
)
```

See the [Data classes API](../api/dataclass.md) for details.

## Python API

### Single RMSD (`run`)

`run` returns a detailed `mobbrmsd_result` for a single pair of structures.

```python
import numpy
from mobbrmsd.dataclass import molecules, molecular_system
from mobbrmsd.mobbrmsd import mobbrmsd

mols = molecular_system([molecules(n_apm=3, n_mol=4, sym=[[0, 2, 1]])])
mrmsd = mobbrmsd(mols=mols)

# x and y: shape [n_atoms, 3] or [n_frames, n_atoms, 3]
result = mrmsd.run(x, y)
print(result.rmsd)
```

Additional options (all optional, keyword only):

| option | default | description |
|---|---|---|
| `cutoff` | `inf` | Stop when lowerbound ≥ cutoff |
| `ub_cutoff` | `inf` | Stop when upperbound ≥ ub_cutoff |
| `difflim` | `0.0` | Stop when upper − lower ≤ difflim |
| `maxeval` | `-1` | Max nodes to evaluate (< 0 for exact) |
| `remove_com` | `True` | Remove centre of mass |
| `sort_by_g` | `True` | Sort reference by self-dispersion for better pruning |
| `get_rotation` | `False` | Also compute rotation matrix |
| `rotate_y` | `False` | Apply optimal rotation to target |

### Batch RMSD (`batch_run`)

`batch_run` returns an RMSD distance matrix.

- `batch_run(x)` → symmetric matrix *D(X, X)*, shape `[n_X, n_X]`
- `batch_run(x, y)` → matrix *D(X, Y)*, shape `[n_X, n_Y]`

```python
mrmsd = mobbrmsd(mols=mols)
D = mrmsd.batch_run(x)          # symmetric distance matrix
D_xy = mrmsd.batch_run(x, y)   # cross distance matrix
```

`batch_run`-specific options:

| option | default | description |
|---|---|---|
| `verbose` | `True` | Show tqdm progress bar for large runs |
| `n_chunk` | `0` | Max pairs per batch (0 = all at once) |

All options accepted by `run` (except `get_rotation` and `rotate_y`) can also be
passed to `batch_run` via `**run_kwargs`:

```python
D = mrmsd.batch_run(x, cutoff=2.0, maxeval=1000, verbose=False)
```

### MDTraj auto-detection

When coordinates are loaded via [MDTraj](https://mdtraj.org/), the molecular
system can be constructed automatically from the topology.

#### `from_topology` (no coordinates needed)

Groups atoms by residue name and sorts species by ascending molecule count,
then descending atoms-per-molecule as a tiebreaker — a structural heuristic
that improves branch-and-bound pruning without requiring coordinates.
All residues of the same name must have identical atom counts; otherwise a
`ValueError` is raised (typical for all-atom peptides whose terminal residues
differ from the interior).

```python
import mdtraj
from mobbrmsd.dataclass import molecular_system
from mobbrmsd.mobbrmsd import mobbrmsd

traj = mdtraj.load("simulation.dcd", top="topology.pdb")
coords = traj.xyz * 10.0  # nm → Å

mols = molecular_system.from_topology(traj.topology)
mrmsd = mobbrmsd(mols=mols)
D = mrmsd.batch_run(coords)
```

#### `from_mdtraj` (variance-based ordering)

When reference coordinates are provided, species are instead sorted by
descending COM variance (most spatially dispersed species first):

```python
mols = molecular_system.from_mdtraj(traj.topology, coords)
```

Without coordinates, `from_mdtraj` falls back to the same structural ordering
as `from_topology`.

#### Generating JSON from topology

`to_dict()` serialises a `molecular_system` to a dict compatible with the
CLI JSON input format, including the correct `atom_order`:

```python
import json
mols = molecular_system.from_topology(traj.topology)
print(json.dumps(mols.to_dict(), indent=2))
```

Output:
```json
{
  "molecular_system": {
    "molecules": [{"n_apm": 3, "n_mol": 4, "name": "WAT"}],
    "atom_order": [0, 1, 2, 3, ...]
  }
}
```

### Atom ordering (`atom_order` / `sorted_by_variance`)

`molecular_system.atom_order` is a permutation index array applied to
coordinates before they reach the Fortran core.  It can be set directly,
generated from a topology, or computed from coordinates:

```python
# Variance-based sort — requires reference coordinates
mols_sorted = mols.sorted_by_variance(coords)
mrmsd = mobbrmsd(mols=mols_sorted)
```

## Command-line interface

`mobbrmsd` exposes three subcommands:

| Subcommand | Description |
|------------|-------------|
| `run [-i JSON] [ref [tgt]] [-p KEY=VALUE ...]` | Run RMSD; topology from JSON, coordinates as positional args |
| `run -b [ref [tgt]]` | Launch the topology builder, then run RMSD |
| `run … -o FILE [-O]` | Write rotated target structure(s) to FILE (MDTraj format) |
| `build [-o FILE]` | Launch the topology builder and save JSON |
| `build -f FILE [-o FILE]` | Parse topology from an MDTraj-supported file, confirm, then save |
| `demo` | Run bundled demonstration scripts |

See [Input format](input.md) for the JSON schema used by `run -i`, and
[Topology Builder](builder.md) for the interactive topology editor.

## Structure output (`-o` / `-O`)

When `-o FILE` is given to `run`, the target structure is rotated to best
superimpose on the reference using the optimal permutation and rotation found
by mobbRMSD, then saved via MDTraj.
Any MDTraj-supported format is accepted (`.pdb`, `.xtc`, `.dcd`, …).

### Output file naming

| Reference frames | Target frames | Output |
|-----------------|---------------|--------|
| 1 | 1 | `FILE` — single file, no numbering |
| 1 | *N* | `FILE` — single multi-frame file |
| *M* | *N* | `FILE.001`, `FILE.002`, …, `FILE.`*M* — one file per reference, each with *N* frames |
| *M* | — *(self)* | `FILE.001`, …, `FILE.`*M* — one file per reference, all frames fitted to it |

When numbering is applied, the counter is zero-padded to at least three digits
(e.g., `.001`) and widens automatically if more than 999 reference frames are
present.

### Overwrite protection

If an output file already exists, `mobbrmsd` prints a warning to stderr and
skips that file without failing:

```
Warning: fitted.pdb already exists; skipping. Use -O to overwrite.
```

Pass `-O` (capital letter) to overwrite existing files silently.
