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
    sym=[[1, 3, 2]],
    name="Water",
)

system = molecular_system(
    mols=[water],
)
```

See the [Data classes API](../api/dataclass.md) for details.

