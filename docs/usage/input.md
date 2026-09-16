# Input format

The command-line interface uses JSON input.

A minimal example is:

```json
{
  "reference": "./path/to/file1.pdb",
  "target": "./path/to/file2.xyz",
  "mols": [
    {
      "n_apm": 2,
      "n_mol": 1,
      "name": "HydrogenFluoride"
    },
    {
      "n_apm": 3,
      "n_mol": 4,
      "sym": [[1, 3, 2]],
      "name": "Water"
    }
  ]
}
```

## Molecular species

Each entry in `mols` describes one molecular species.

### `n_apm`

Number of atoms per molecule.

### `n_mol`

Number of molecules of this species.

### `sym`

Optional intramolecular atomic permutations.

For example:

```json
"sym": [[1, 3, 2]]
```

represents a permutation of the atoms within a molecule.

The identity permutation does not need to be specified.

### `name`

Optional name of the molecular species.

## Coordinate files

The coordinate files are loaded using MDTraj.

Only Cartesian coordinates are used. Atom names, residue names,
and other topology information are not used for the RMSD calculation.

## Running

```bash
python -m mobbrmsd run -i input.json
```

