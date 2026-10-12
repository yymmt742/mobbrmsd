# Quick Start

## Installation

```bash
pip install git+https://github.com/yymmt742/mobbrmsd.git
```

For full installation options (Windows, development install) see [Installation](installation.md).

---

## Try the demo

```bash
python -m mobbrmsd demo
```

An interactive menu lets you choose from several bundled demonstrations.

---

## Command-line interface

### Step 1 — Build a topology file

The topology JSON describes how atoms are grouped into molecular species.
Generate it from any MDTraj-supported coordinate file:

```bash
python -m mobbrmsd build -f system.pdb -o topology.json
```

`mobbrmsd` groups atoms by residue name and displays the detected species for
confirmation before writing the file.  Edit interactively if needed (`e` at
the prompt launches the topology builder).

??? tip "Generate from Python"
    ```python
    import json, mdtraj
    from mobbrmsd.dataclass import molecular_system

    top = mdtraj.load_topology("system.pdb")
    mols = molecular_system.from_topology(top)
    print(json.dumps(mols.to_dict(), indent=2))
    ```

### Step 2 — Run an RMSD calculation

```bash
python -m mobbrmsd run -i topology.json reference.pdb target.pdb
```

**Single-pair output** (`nref = 1, ntrg = 1`):

```
Reference frame       1 | Target    frame       1
  RMSD        :            1.234567890
  ----------
  AUTO CORR.  :            2.345678901
  UPPER BOUND :            1.234567890  (as RMSD) :            1.234567890
  LOWER BOUND :            0.000000000  (as RMSD) :            0.000000000
  N EVAL      :                  12345  (ratio)   :  0.012345678901234568
  ----------
```

**Batch output** (multi-frame files → space-separated `ref trg rmsd` per line):

```bash
python -m mobbrmsd run -i topology.json reference.dcd target.dcd \
       -p topology=topology.pdb
```

For the full list of `-p` parameters (`maxeval`, `cutoff`, `mask`, …) see
[Input format → Run parameters](input.md#run-parameters).

### Step 3 — Save the fitted structure

Add `-o FILE` to write the target structure rotated onto the reference:

```bash
# single pair → fitted.pdb
python -m mobbrmsd run -i topology.json reference.pdb target.pdb -o fitted.pdb

# multiple reference frames → fitted.001.pdb, fitted.002.pdb, ...
python -m mobbrmsd run -i topology.json reference.dcd target.pdb -o fitted.pdb

# no target (self-comparison) → fit every frame to each reference frame
python -m mobbrmsd run -i topology.json trajectory.dcd -o frame.pdb
```

Output file naming:

| Situation | Output files |
|-----------|-------------|
| `nref = 1` | `fitted.pdb` (no numbering) |
| `nref > 1` | `fitted.001.pdb`, `fitted.002.pdb`, … |
| `nref > 1`, `ntrg > 1` | one multi-frame file per reference frame |

If an output file already exists, `mobbrmsd` prints a warning and skips it.
Use `-O` to overwrite without prompting.

---

## Python API

### Minimal example

```python
import numpy
from mobbrmsd.dataclass import molecules, molecular_system
from mobbrmsd.mobbrmsd import mobbrmsd

# 4 water molecules; H atoms (positions 1 and 2) are interchangeable
water = molecules(n_apm=3, n_mol=4, sym=[[0, 2, 1]])
mols = molecular_system([water])
mrmsd = mobbrmsd(mols=mols)

x = numpy.random.randn(12, 3)   # reference, shape [n_atoms, 3]
y = numpy.random.randn(12, 3)   # target

result = mrmsd.run(x, y)
print(f"RMSD = {result.rmsd():.6f} Å")
```

### Batch RMSD matrix

```python
frames = numpy.random.randn(50, 12, 3)

D = mrmsd.batch_run(frames)          # symmetric [50, 50]
D_xy = mrmsd.batch_run(frames, y)    # cross [50, 1]
```

### Load from MDTraj

```python
import mdtraj
from mobbrmsd.dataclass import molecular_system
from mobbrmsd.mobbrmsd import mobbrmsd

traj = mdtraj.load("simulation.dcd", top="topology.pdb")
coords = traj.xyz * 10.0   # nm → Å

# Topology auto-detection
mols = molecular_system.from_mdtraj(traj.topology, coords)
mrmsd = mobbrmsd(mols=mols)

D = mrmsd.batch_run(coords)          # [n_frames, n_frames] RMSD matrix
```

### Get the rotated structure

```python
y_copy = y.copy()
result = mrmsd.run(x, y_copy, rotate_y=True)
# y_copy now holds the optimally rotated target (mobbrmsd atom order)
```

---

## Next steps

- [Input format](usage/input.md) — full JSON schema and `-p` parameter reference
- [Topology Builder](usage/builder.md) — interactive topology editor
- [API Reference](api/index.md) — Python API details
- [Theory](theory/index.md) — algorithm and mathematical background
