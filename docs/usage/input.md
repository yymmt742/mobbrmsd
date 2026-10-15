# Input format

The command-line interface uses JSON input.

## Minimal example

```json
{
  "reference": "./path/to/file1.pdb",
  "target": "./path/to/file2.xyz",
  "molecular_system": {
    "molecules": [
      {"n_apm": 2, "n_mol": 1, "name": "HydrogenFluoride"},
      {"n_apm": 3, "n_mol": 4, "sym": [[0, 2, 1]], "name": "Water"}
    ]
  }
}
```

In this example, the system contains one hydrogen fluoride and four waters.

## `molecular_system`

Describes the molecular topology.

### `molecules`

A list of molecular-species entries (see [Molecular species](#molecular-species) below).
A single-species system may also use a bare object instead of a list.

### `name`

Optional name for the molecular system.

### `atom_order`

Optional permutation index array.  When present, coordinates are reindexed
by `coords[atom_order]` before being passed to the Fortran core.
This is useful when the coordinate file has a different atom ordering than
the species-grouped layout expected by mobbRMSD.

```json
"atom_order": [0, 3, 6, 9, 1, 4, 7, 10, 2, 5, 8, 11]
```

You can generate the correct `atom_order` from Python with
`molecular_system.to_dict()` (see [Generating JSON from topology](#generating-json-from-topology)).

## Molecular species

Each entry in `molecules` describes one molecular species.

### `n_apm`

Number of atoms per molecule.

### `n_mol`

Number of molecules of this species.

### `sym`

Optional intramolecular atomic permutations.

```json
"sym": [[0, 2, 1]]
```

Represents a permutation of the atoms within a molecule (0-based atom indices).
The identity permutation does not need to be specified.

### `name`

Optional name of the molecular species.

## Coordinate files

The coordinate files are loaded using MDTraj.

Only Cartesian coordinates are used. Atom names, residue names,
and other topology information are not used for the RMSD calculation.

## Generating JSON from topology

`molecular_system.to_dict()` serialises a `molecular_system` object to a
dict that matches the JSON format above.  Combined with `from_topology`, this
lets you generate a ready-to-use input JSON from any MDTraj-supported
topology file:

```python
import json
import mdtraj
from mobbrmsd.dataclass import molecular_system

top = mdtraj.load_topology("topology.pdb")   # PDB, prmtop, gro: self-contained
mols = molecular_system.from_topology(top)
print(json.dumps(mols.to_dict(), indent=2))
```

`from_topology` groups atoms by residue name (same residue name → same
species), validates that all residues of a species have the same number of
atoms, and orders species by ascending molecule count then descending
atoms-per-molecule (tiebreaker).

If residues of the same name have different atom counts (e.g. N/C-terminal
residues in an all-atom peptide), `from_topology` raises a `ValueError`.
In that case define the species manually with the [Topology Builder](builder.md)
or via the Python API.

The CLI command `build -f FILE` runs this pipeline interactively and prompts
before writing the JSON (see [Topology Builder → build subcommand](builder.md#build-subcommand)).

## Running

Pass the topology JSON with `-i` and the coordinate files as positional arguments:

```bash
# topology from JSON, coordinates as positional args
python -m mobbrmsd run -i topology.json ref.pdb tgt.pdb

# everything in one JSON file
python -m mobbrmsd run -i input.json

# coordinate files without topology JSON (single all-atom species assumed)
python -m mobbrmsd run ref.pdb tgt.pdb

# extra run parameters via -p (repeatable)
python -m mobbrmsd run -i topology.json ref.pdb tgt.pdb -p maxeval=1000 -p cutoff=2.0
```

### Run parameters (`-p KEY=VALUE`)

Extra keyword arguments are passed to the RMSD calculation via `-p`.

| Key | Example | Description |
|-----|---------|-------------|
| `precision` | `single` | Use single-precision floats (`float32`) |
| `mask` | `name CA` | MDTraj atom selection mask applied when loading coordinates |
| `topology` | `top.pdb` | Topology file for trajectory formats that do not embed it (DCD, XTC, …) |
| `maxeval` | `1000` | Maximum branch-and-bound nodes (< 0 = exact) |
| `cutoff` | `2.0` | Stop when lower bound ≥ cutoff |
| `ub_cutoff` | `2.0` | Stop when upper bound ≥ ub_cutoff |
| `difflim` | `0.001` | Stop when upper − lower ≤ difflim |

### Dispatch rules

| Arguments | Behavior |
|-----------|----------|
| `reference` only | Symmetric `batch_run(x)` — *n_frame* × *n_frame* distance matrix |
| `reference target`, both 1 frame | Single-pair `run(x, y)` with full output |
| `reference target`, multi-frame | Cross `batch_run(x, y)` — *n_ref* × *n_trg* distance matrix |

When `-o FILE` is given, the dispatch changes for multi-frame or self-comparison cases:
instead of the fast Fortran batch routine, `run(x, y, rotate_y=True)` is called
per pair in Python so that the rotated coordinates can be captured.
See [Structure output](index.md#structure-output-o-o) for the full naming rules.

