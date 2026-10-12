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
      {"n_apm": 3, "n_mol": 4, "sym": [[1, 3, 2]], "name": "Water"}
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
"sym": [[1, 3, 2]]
```

Represents a permutation of the atoms within a molecule.
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

top = mdtraj.load_topology("topology.pdb")
mols = molecular_system.from_topology(top)
print(json.dumps(mols.to_dict(), indent=2))
```

`from_topology` groups atoms by residue name (same residue name → same
species), validates that all residues of a species have the same number of
atoms, and orders species by descending atoms-per-molecule then ascending
molecule count.

## Running

```bash
python -m mobbrmsd run -i input.json
```

