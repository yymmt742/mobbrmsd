# Topology Builder

The Topology Builder is an interactive editor for defining molecular systems.
It provides a command-line wizard, a Python API, and can save the resulting
topology as a JSON file compatible with `run -i`.

## Quick start

### Build, save, then run

```bash
# 1. Define topology interactively and save to JSON
python -m mobbrmsd build -o topology.json

# 2. Run with the saved topology and coordinate files as positional args
python -m mobbrmsd run -i topology.json ref.pdb tgt.pdb
```

### Build and run in one step

```bash
python -m mobbrmsd run -b ref.pdb tgt.pdb
```

The `-b` / `--build` flag launches the builder, then runs the RMSD
calculation using the given coordinate files.  If `reference` is omitted,
the tool prompts for it interactively after the builder exits.

---

## `build` subcommand

```
python -m mobbrmsd build [-o FILE] [-f FILE [--topology TOP]]
```

| Option | Description |
|--------|-------------|
| `-o FILE` / `--output FILE` | Write JSON to FILE (default: print to stdout) |
| `-f FILE` / `--file FILE` | Parse topology from an MDTraj-supported file (PDB, GRO, …) instead of the interactive builder |
| `--topology TOP` | Topology file required by some formats (e.g. a PDB alongside a DCD) |

When `-f` is given, the parsed topology is shown as a species table and
the user is prompted before writing:

```
  Discovered 3 species
  #    Name          n_apm  n_mol  sym
  ───  ────────────  ─────  ─────  ───
    1  HOH               3    216  —
    2  NA                1     23  —
    3  CL                1     18  —
  Total atoms: 689

  Proceed? [Y]es / [n]o / [e]dit >
```

Choosing `e` launches the interactive REPL with the discovered species
pre-loaded so they can be inspected or edited before saving.
If `molecular_system.from_topology()` raises an error (e.g. inconsistent
atom counts), the user is offered to define species manually via the REPL.

The output uses the same `"molecular_system"` format as `molecular_system.to_dict()`
and is directly consumable by `run -i`.

---

## Interactive REPL

On startup the REPL prints a summary table of defined species and accepts the
following commands.  All commands accept their first letter as a shorthand.

| Command | Shorthand | Description |
|---------|-----------|-------------|
| `list` | `l` | Show the species summary table |
| `list N` | `l N` | Show species N in detail (atom names, sym) |
| `add` | `a` | Start the add-species wizard |
| `edit N` | `e N` | Re-run the wizard for species N |
| `remove N` | `r N` | Delete species N |
| `move N M` | `mv N M` | Move species N to position M |
| `sort` | `s` | Sort species (ascending n_mol, descending n_apm) |
| `name [NAME]` | `n` | Set or clear the system name |
| `order LIST` | `o` | Set atom_order (space-separated index list) |
| `done` | `d` | Build and exit |
| `quit` | `q` | Exit without building |

### Add-species wizard

`add` steps through five fields.  Press **Up** to go back one step; press
**Ctrl+C** to cancel.

| Step | Field | Description |
|------|-------|-------------|
| 1 | `n_apm` | Atoms per molecule |
| 2 | `n_mol` | Number of molecules |
| 3 | `name` | Species name (optional; Enter to skip) |
| 4 | `atom_names` | Per-atom labels (optional; Ctrl+C to skip) |
| 5 | `sym` | Intramolecular permutation symmetry |

### Atom name editor

Atom names are short labels (up to 4 ASCII characters) displayed in the sym
editor.  They are **not** included in the output JSON or `molecular_system`.

| Key | Action |
|-----|--------|
| **Left** / **Right** | Move cursor between atom positions |
| *name* + **Enter** | Set the name at the cursor and advance |
| **Enter** (blank) | Auto-generate a placeholder (`X001`, `X002`, …) and advance |
| **Enter** at last position / **Down** / **Ctrl+D** | Confirm all and exit |
| **Up** | Go back to the `name` step |

Validation rules:

- Must be ASCII only.
- Must be unique within the species.
- Names longer than 4 characters are accepted but truncated to 4 with a warning.
- Empty strings and non-ASCII input are rejected.

### Sym editor

Permutations are indexed from 0 (identity, always present).  User-defined
permutations start at index 1.

| Command | Description |
|---------|-------------|
| `add` / `a` | Open the permutation editor to define a new sym |
| `remove N` / `r N` / `N` | Remove user-defined permutation N (N ≥ 1) |
| `list` / `l` | Redisplay the sym list |
| `done` / `d` / **Enter** | Confirm and exit |
| `back` / `b` / **Up** | Return to the atom_names step |

#### Permutation editor

Displays atoms as `Name` / `Map` rows.  For molecules with many atoms the
display is paginated — one chunk of atoms per screen — with a `(n/N)` page
indicator at the end of each row.

| Key / input | Action |
|-------------|--------|
| **Left** / **Right** | Move the cursor to another atom position (view auto-follows) |
| **Up** / **Down** | Scroll the display page without moving the cursor |
| *N* + **Enter** | Swap the atom at the cursor with atom N |
| **Enter** (blank) | Confirm the permutation |
| **Ctrl+C** / **Ctrl+D** | Cancel (discard edits) |

---

## Python API

### Programmatic construction

```python
from mobbrmsd.run.topology_builder import TopologyBuilder

builder = (
    TopologyBuilder()
    .add(n_apm=3, n_mol=4, sym=[[1, 0, 2]], name="Water")
    .add(n_apm=6, n_mol=2, name="Methanol")
)

mols = builder.build()           # -> molecular_system
```

**`TopologyBuilder` methods** (all return `self` for chaining):

| Method | Description |
|--------|-------------|
| `.add(n_apm, n_mol, *, sym, name, atom_names)` | Append a species |
| `.edit(index, **kwargs)` | Update fields of species `index` |
| `.remove(index)` | Delete species at `index` |
| `.sort()` | Sort by ascending n_mol then descending n_apm |
| `.set_name(name)` | Set the system name |
| `.set_atom_order(order)` | Set the atom-order permutation |
| `.build()` | Return a frozen `molecular_system` |
| `.to_dict()` | Serialise to a JSON-compatible dict |
| `.from_dict(d)` | *(classmethod)* Restore from dict |
| `.run()` | Launch the interactive REPL |

### Saving and loading

```python
import json
from mobbrmsd.run.topology_builder import TopologyBuilder
from mobbrmsd.run.parse import load

# Build and save
mols = TopologyBuilder().add(n_apm=3, n_mol=4, sym=[[1, 0, 2]]).build()
with open("topology.json", "w") as f:
    json.dump(mols.to_dict(), f, indent=2)

# Load and pass to mobbrmsd
loaded = load(json.load(open("topology.json")))
```

---

## JSON output format

```json
{
  "molecular_system": {
    "molecules": [
      {"n_apm": 3, "n_mol": 4, "sym": [[1, 0, 2]], "name": "Water"},
      {"n_apm": 6, "n_mol": 2, "name": "Methanol"}
    ],
    "name": "MySystem",
    "atom_order": [0, 3, 6, 9, 1, 4, 7, 10, 2, 5, 8, 11]
  }
}
```

`name`, `sym`, and `atom_order` are omitted when `None`.

The format is identical to `molecular_system.to_dict()` and is compatible with
the `run -i` input format described in [Input format](input.md).
