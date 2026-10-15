"""
Download MDTraj test files and verify molecular_system.from_topology parsing.
Run from the repo root:  .venv\Scripts\python.exe example\test_parse_topology.py
"""
import sys
import urllib.request
import pathlib

HERE = pathlib.Path(__file__).parent
DATA_DIR = HERE / "test_data"
DATA_DIR.mkdir(exist_ok=True)

FILES = {
    "4waters.pdb": "https://raw.githubusercontent.com/mdtraj/mdtraj/master/tests/data/4waters.pdb",
    "ala_ala_ala.pdb": "https://raw.githubusercontent.com/mdtraj/mdtraj/master/tests/data/ala_ala_ala.pdb",
    "alanine-dipeptide-explicit.pdb": "https://raw.githubusercontent.com/mdtraj/mdtraj/master/tests/data/alanine-dipeptide-explicit.pdb",
    "alanine-dipeptide-explicit.prmtop": "https://raw.githubusercontent.com/mdtraj/mdtraj/master/tests/data/alanine-dipeptide-explicit.prmtop",
    "frame0.gro": "https://raw.githubusercontent.com/mdtraj/mdtraj/master/tests/data/frame0.gro",
}

print("=== Downloading test files ===")
for fname, url in FILES.items():
    dest = DATA_DIR / fname
    if dest.exists():
        print(f"  {fname}: cached")
        continue
    try:
        urllib.request.urlretrieve(url, dest)
        print(f"  {fname}: downloaded ({dest.stat().st_size} bytes)")
    except Exception as e:
        print(f"  {fname}: FAILED — {e}")

print()

try:
    import mdtraj
    from mobbrmsd.dataclass import molecular_system
except ImportError as e:
    print(f"Import error: {e}")
    sys.exit(1)


def check(label, file, top=None):
    print(f"--- {label} ---")
    path = DATA_DIR / file
    if not path.exists():
        print("  SKIP: file not found")
        return
    try:
        if top:
            # Trajectory format — load as traj to get topology
            traj = mdtraj.load(str(path), top=str(DATA_DIR / top))
            topology = traj.topology
        else:
            topology = mdtraj.load_topology(str(path))
        print(f"  MDTraj: {topology.n_atoms} atoms, "
              f"{topology.n_residues} residues, "
              f"{topology.n_chains} chains")
        mols = molecular_system.from_topology(topology)
        d = mols.to_dict().get("molecular_system", {})
        species = d.get("molecules", [])
        if isinstance(species, dict):
            species = [species]
        for sp in species:
            name = sp.get("name") or "—"
            print(f"    {name:>12}  n_apm={sp['n_apm']:3}  n_mol={sp['n_mol']:3}"
                  f"  total={sp['n_apm']*sp['n_mol']:4}")
        ao = d.get("atom_order")
        print(f"  atom_order: {len(ao)} entries" if ao else "  atom_order: none")
        print("  OK")
    except Exception as e:
        import traceback
        print(f"  ERROR: {e}")
        traceback.print_exc()
    print()


# Built-in example file
print("=== Parsing test ===\n")

print("--- PDB (local) — HF_4H2O_A.pdb ---")
try:
    topology = mdtraj.load_topology(str(HERE / "HF_4H2O_A.pdb"))
    print(f"  MDTraj: {topology.n_atoms} atoms, "
          f"{topology.n_residues} residues, "
          f"{topology.n_chains} chains")
    mols = molecular_system.from_topology(topology)
    d = mols.to_dict().get("molecular_system", {})
    species = d.get("molecules", [])
    if isinstance(species, dict):
        species = [species]
    for sp in species:
        name = sp.get("name") or "—"
        print(f"    {name:>12}  n_apm={sp['n_apm']:3}  n_mol={sp['n_mol']:3}"
              f"  total={sp['n_apm']*sp['n_mol']:4}")
    ao = d.get("atom_order")
    print(f"  atom_order: {len(ao)} entries" if ao else "  atom_order: none")
    print("  OK")
except Exception as e:
    import traceback
    print(f"  ERROR: {e}")
    traceback.print_exc()
print()

check("PDB — 4 waters (MDTraj test)", "4waters.pdb")
check("PDB — ala-ala-ala (3 ALA)", "ala_ala_ala.pdb")
check("PDB — alanine-dipeptide explicit solvent", "alanine-dipeptide-explicit.pdb")
check("AMBER prmtop — alanine-dipeptide explicit", "alanine-dipeptide-explicit.prmtop")
check("GROMACS gro — frame0", "frame0.gro")
