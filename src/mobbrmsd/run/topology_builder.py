"""Interactive builder for molecular_system.

Since ``molecular_system`` is a frozen dataclass, this module provides a
mutable ``TopologyBuilder`` that accumulates species edits interactively
and emits a ``molecular_system`` when the user is done.

``TopologyBuilder`` itself has no dependency on :mod:`mobbrmsd.dataclass`.
Conversion to/from ``molecular_system`` is handled by helper functions in
:mod:`mobbrmsd.run.parse` (``to_molecular_system``, ``builder_from_topology``,
``builder_from_mdtraj``), keeping the builder a pure data-editing layer.
"""

from __future__ import annotations

import copy
import json
from typing import Dict, List, Optional


_SEP = "  " + "-" * 72


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _parse_sym(raw: str) -> Optional[List[List[int]]]:
    """Parse a sym string like '[[1,2,0],[2,0,1]]' into a list of lists."""
    if raw.strip().lower() in ("", "none", "null"):
        return None
    try:
        val = json.loads(raw)
    except json.JSONDecodeError:
        raise ValueError(f"Cannot parse sym — expected JSON list of lists, got: {raw!r}")
    if isinstance(val, list) and all(isinstance(v, list) for v in val):
        return val
    raise ValueError(f"sym must be a list of lists, got: {val!r}")


def _ask(prompt: str, default: str = "") -> str:
    hint = f" (default: {default!r})" if default != "" else ""
    raw = input(f"    {prompt}{hint} >> ").strip()
    return raw if raw != "" else default


def _ask_int(prompt: str, default: int) -> int:
    while True:
        raw = _ask(prompt, str(default))
        try:
            return int(raw)
        except ValueError:
            print(f"    ! expected integer, got {raw!r}")


def _ask_optional_str(prompt: str, default: Optional[str]) -> Optional[str]:
    raw = _ask(prompt, "" if default is None else default)
    return None if raw == "" else raw


# ---------------------------------------------------------------------------
# TopologyBuilder
# ---------------------------------------------------------------------------

class TopologyBuilder:
    """Mutable builder for :class:`~mobbrmsd.dataclass.molecular_system`.

    ``TopologyBuilder`` is a pure data-editing layer with no dependency on
    the ``mobbrmsd`` dataclasses.  Conversion to/from ``molecular_system``
    is done via :mod:`mobbrmsd.run.parse`:

    * :func:`~mobbrmsd.run.parse.to_molecular_system` — builder → molecular_system
    * :func:`~mobbrmsd.run.parse.builder_from_topology` — MDTraj topology → builder
    * :func:`~mobbrmsd.run.parse.builder_from_mdtraj` — MDTraj topology + coords → builder

    Examples
    --------
    Start from a topology file and edit interactively::

        from mobbrmsd.run.parse import builder_from_topology
        import mdtraj

        top = mdtraj.load_topology("topology.pdb")
        mols = builder_from_topology(top).run()

    Programmatic fluent API::

        mols = (
            TopologyBuilder()
            .add(n_apm=3, n_mol=4, name="WAT")
            .add(n_apm=2, n_mol=1, name="HF")
            .sort()
            .build()
        )
    """

    def __init__(self) -> None:
        self._species: List[Dict] = []
        self._name: Optional[str] = None
        self._atom_order: Optional[List[int]] = None

    # ------------------------------------------------------------------
    # constructor from duck-typed molecular_system
    # ------------------------------------------------------------------

    @classmethod
    def from_molecular_system(cls, mols) -> "TopologyBuilder":
        """Initialise from any object with ``.mols``, ``.name``, ``.atom_order``.

        Works with :class:`~mobbrmsd.dataclass.molecular_system` without
        importing it — pure duck typing.
        """
        b = cls()
        b._name = mols.name
        b._atom_order = (
            mols.atom_order.tolist() if mols.atom_order is not None else None
        )
        for m in mols.mols:
            b._species.append(
                {
                    "n_apm": m.n_apm,
                    "n_mol": m.n_mol,
                    "sym": copy.deepcopy(m.sym),
                    "name": m.name,
                }
            )
        return b

    # ------------------------------------------------------------------
    # programmatic API (fluent)
    # ------------------------------------------------------------------

    def add(
        self,
        n_apm: int,
        n_mol: int,
        sym: Optional[List[List[int]]] = None,
        name: Optional[str] = None,
    ) -> "TopologyBuilder":
        """Append a new species."""
        self._species.append({"n_apm": n_apm, "n_mol": n_mol, "sym": sym, "name": name})
        return self

    def remove(self, index: int) -> "TopologyBuilder":
        """Remove species at *index* (0-based)."""
        self._species.pop(index)
        return self

    def edit(self, index: int, **kwargs) -> "TopologyBuilder":
        """Update fields of species at *index*.

        Accepted keyword arguments: ``n_apm``, ``n_mol``, ``sym``, ``name``.
        """
        sp = self._species[index]
        for k in ("n_apm", "n_mol", "sym", "name"):
            if k in kwargs:
                sp[k] = kwargs[k]
        return self

    def reorder(self, indices: List[int]) -> "TopologyBuilder":
        """Reorder species according to *indices* (a permutation of 0..n-1)."""
        if sorted(indices) != list(range(len(self._species))):
            raise ValueError("indices must be a permutation of 0 .. n_species-1")
        self._species = [self._species[i] for i in indices]
        return self

    def sort(self) -> "TopologyBuilder":
        """Apply structural sort (ascending n_mol, descending n_apm)."""
        indices = sorted(
            range(len(self._species)),
            key=lambda i: (self._species[i]["n_mol"], -self._species[i]["n_apm"]),
        )
        self._species = [self._species[i] for i in indices]
        return self

    def set_name(self, name: Optional[str]) -> "TopologyBuilder":
        """Set the system name."""
        self._name = name
        return self

    def set_atom_order(self, order: Optional[List[int]]) -> "TopologyBuilder":
        """Set a custom atom_order permutation."""
        self._atom_order = list(order) if order is not None else None
        return self

    # ------------------------------------------------------------------
    # build — thin wrapper; conversion logic lives in parse.py
    # ------------------------------------------------------------------

    def build(self):
        """Return the final frozen ``molecular_system``.

        Delegates to :func:`mobbrmsd.run.parse.to_molecular_system`.
        """
        from .parse import to_molecular_system
        return to_molecular_system(self)

    # ------------------------------------------------------------------
    # display
    # ------------------------------------------------------------------

    def _show(self) -> None:
        print(_SEP)
        print(f"  System name : {self._name}" if self._name else "  System name : (none)")
        print(
            f"  atom_order  : {self._atom_order}"
            if self._atom_order is not None
            else "  atom_order  : (none)"
        )
        print(_SEP)
        if not self._species:
            print("  (no species defined)")
        else:
            print(f"  {'#':>3}  {'n_apm':>6}  {'n_mol':>6}  {'name':<16}  sym")
            for i, sp in enumerate(self._species):
                print(
                    f"  {i:>3}  {sp['n_apm']:>6}  {sp['n_mol']:>6}"
                    f"  {sp['name'] or '(none)':<16}  "
                    f"{repr(sp['sym']) if sp['sym'] else '(none)'}"
                )
        print(_SEP)

    # ------------------------------------------------------------------
    # interactive REPL
    # ------------------------------------------------------------------

    _HELP = """\
  Commands:
    list  (l)          — show current species
    add   (a)          — add a new species
    edit  (e) N        — edit species N
    remove(r) N        — remove species N
    move  (mv) N M     — move species N to position M
    sort  (s)          — apply structural sort (n_mol asc, n_apm desc)
    name  (n) [TEXT]   — set system name (empty to clear)
    order (o) [JSON]   — set atom_order from a JSON int list (empty to clear)
    done  (d)          — finish editing and return molecular_system
    quit  (q)          — abort without returning
    help  (h / ?)      — show this help
"""

    def run(self):
        """Start the interactive REPL.

        Returns the built ``molecular_system`` on ``done``, or ``None`` on
        ``quit``.
        """
        print(_SEP)
        print("  mobbrmsd topology editor  (type 'help' for commands)")
        self._show()

        while True:
            try:
                raw = input("  >> ").strip()
            except (EOFError, KeyboardInterrupt):
                print("\n  Aborted.")
                return None

            if not raw:
                continue

            parts = raw.split()
            cmd = parts[0].lower()
            args = parts[1:]

            if cmd in ("help", "h", "?"):
                print(self._HELP)

            elif cmd in ("list", "l"):
                self._show()

            elif cmd in ("add", "a"):
                print("  --- add species ---")
                n_apm = _ask_int("n_apm", 1)
                n_mol = _ask_int("n_mol", 1)
                name = _ask_optional_str("name (leave blank for none)", None)
                try:
                    sym = _parse_sym(
                        _ask("sym as JSON (e.g. [[1,0,2]]), leave blank for none", "")
                    )
                except ValueError as e:
                    print(f"  ! {e}")
                    continue
                self.add(n_apm=n_apm, n_mol=n_mol, sym=sym, name=name)
                self._show()

            elif cmd in ("edit", "e"):
                if len(args) != 1 or not args[0].isdigit():
                    print("  ! usage: edit N")
                    continue
                idx = int(args[0])
                if not (0 <= idx < len(self._species)):
                    print(f"  ! index {idx} out of range")
                    continue
                sp = self._species[idx]
                print(f"  --- edit species {idx} ---")
                n_apm = _ask_int("n_apm", sp["n_apm"])
                n_mol = _ask_int("n_mol", sp["n_mol"])
                name = _ask_optional_str(
                    "name (leave blank to keep, enter '-' to clear)",
                    sp["name"] or "",
                )
                if name == "-":
                    name = None
                try:
                    sym_raw = _ask(
                        "sym as JSON (leave blank to keep, enter 'none' to clear)",
                        json.dumps(sp["sym"]) if sp["sym"] else "",
                    )
                    sym = _parse_sym(sym_raw) if sym_raw != "" else sp["sym"]
                except ValueError as e:
                    print(f"  ! {e}")
                    continue
                self.edit(idx, n_apm=n_apm, n_mol=n_mol, sym=sym, name=name)
                self._show()

            elif cmd in ("remove", "r"):
                if len(args) != 1 or not args[0].isdigit():
                    print("  ! usage: remove N")
                    continue
                idx = int(args[0])
                if not (0 <= idx < len(self._species)):
                    print(f"  ! index {idx} out of range")
                    continue
                self.remove(idx)
                self._show()

            elif cmd in ("move", "mv"):
                if len(args) != 2 or not all(a.isdigit() for a in args):
                    print("  ! usage: move N M")
                    continue
                src, dst = int(args[0]), int(args[1])
                n = len(self._species)
                if not (0 <= src < n and 0 <= dst < n):
                    print(f"  ! indices out of range (0..{n-1})")
                    continue
                sp = self._species.pop(src)
                self._species.insert(dst, sp)
                self._show()

            elif cmd in ("sort", "s"):
                self.sort()
                self._show()

            elif cmd in ("name", "n"):
                self.set_name(" ".join(args) if args else None)
                print(f"  System name set to: {self._name!r}")

            elif cmd in ("order", "o"):
                raw_order = " ".join(args)
                if not raw_order:
                    self.set_atom_order(None)
                    print("  atom_order cleared.")
                else:
                    try:
                        order = json.loads(raw_order)
                        if not isinstance(order, list) or not all(
                            isinstance(x, int) for x in order
                        ):
                            raise ValueError("expected a JSON list of ints")
                        self.set_atom_order(order)
                        print(f"  atom_order set ({len(order)} elements).")
                    except (json.JSONDecodeError, ValueError) as e:
                        print(f"  ! {e}")

            elif cmd in ("done", "d"):
                try:
                    result = self.build()
                    self._show()
                    print("  molecular_system built successfully.")
                    print(_SEP)
                    return result
                except ValueError as e:
                    print(f"  ! {e}")

            elif cmd in ("quit", "q"):
                print("  Aborted.")
                return None

            else:
                print(f"  ! unknown command: {cmd!r}. Type 'help' for commands.")


# ---------------------------------------------------------------------------
# standalone entry point
# ---------------------------------------------------------------------------

def _main() -> None:
    import argparse
    import sys

    ap = argparse.ArgumentParser(
        prog="python -m mobbrmsd.run.topology_builder",
        description=(
            "Interactive molecular_system editor. "
            "On 'done', prints the resulting JSON to stdout."
        ),
    )
    ap.add_argument("-t", "--topology", metavar="FILE", default=None,
                    help="MDTraj topology/structure file to pre-load (PDB, GRO, …).")
    ap.add_argument("-c", "--coords", metavar="FILE", default=None,
                    help="Trajectory for variance-based ordering (requires --topology).")
    ap.add_argument("-n", "--name", metavar="NAME", default=None,
                    help="Initial system name.")
    ap.add_argument("-o", "--output", metavar="FILE", default=None,
                    help="Write JSON to FILE instead of stdout.")
    args = ap.parse_args()

    if args.topology is not None:
        try:
            import mdtraj
        except ImportError:
            print("ERROR: mdtraj is required to load topology files.", file=sys.stderr)
            sys.exit(1)
        try:
            from .parse import builder_from_topology, builder_from_mdtraj
            if args.coords is not None:
                traj = mdtraj.load(args.coords, top=args.topology)
                builder = builder_from_mdtraj(
                    traj.topology, coords=traj.xyz * 10.0, name=args.name
                )
            else:
                top = mdtraj.load_topology(args.topology)
                builder = builder_from_topology(top, name=args.name)
        except Exception as e:
            print(f"ERROR loading topology: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        builder = TopologyBuilder()
        if args.name:
            builder.set_name(args.name)

    result = builder.run()
    if result is None:
        sys.exit(1)

    output = json.dumps(result.to_dict(), indent=2)
    if args.output is not None:
        with open(args.output, "w") as f:
            f.write(output)
            f.write("\n")
        print(f"  Written to {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    _main()
