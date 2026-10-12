"""Interactive builder for molecular_system.

Since ``molecular_system`` is a frozen dataclass, this module provides a
mutable ``TopologyBuilder`` that accumulates species edits interactively
and emits a ``molecular_system`` when the user is done.

``TopologyBuilder`` itself has no dependency on :mod:`mobbrmsd.dataclass`.
Conversion to/from ``molecular_system`` is handled by helper functions in
:mod:`mobbrmsd.run.parse`.

Interactive features require ``prompt_toolkit`` (``pip install prompt-toolkit``).
Without it the REPL falls back to plain ``input()`` with ``!back``/``!cancel``
navigation.
"""

from __future__ import annotations

import copy
import json
from typing import Dict, List, Optional

_SEP = "  " + "-" * 72

# ── prompt_toolkit detection ──────────────────────────────────────────────────

try:
    from prompt_toolkit import prompt as _pt_prompt
    from prompt_toolkit.history import InMemoryHistory as _PTHistory
    _PT = True
except ImportError:
    _PT = False


class _Abort(Exception):
    """Cancel the current interactive operation without quitting the REPL."""


# ── prompt wrapper ────────────────────────────────────────────────────────────

class _PromptFn:
    """Thin wrapper around prompt() with optional history.

    Raises ``_Abort`` on Ctrl+C / EOF so callers don't need to catch
    ``KeyboardInterrupt`` directly.
    """

    def __init__(self, use_history: bool = False):
        self._hist = _PTHistory() if (_PT and use_history) else None

    def __call__(self, text: str, default: str = "") -> str:
        if _PT:
            try:
                return _pt_prompt(text, default=default, history=self._hist).strip()
            except KeyboardInterrupt:
                raise _Abort()
            except EOFError:
                raise _Abort()
        else:
            try:
                raw = input(text).strip()
                return raw if raw else default
            except (KeyboardInterrupt, EOFError):
                raise _Abort()


# ── permutation editor ────────────────────────────────────────────────────────

def _col_width(n_apm: int) -> int:
    return max(2, len(str(n_apm - 1)))


def _render_mapping(n_apm: int, mapping: List, cursor_pos: int) -> str:
    w = _col_width(n_apm)
    sep = " "
    pos_row = "  Pos: " + sep.join(f"{i:{w}d}" for i in range(n_apm))
    map_row = "  Map: " + sep.join(
        f"{v:{w}d}" if v is not None else " " * (w - 1) + "_" for v in mapping
    )
    # cursor arrow: each column is (w + 1) chars wide except the last
    cur_row = "  Cur: " + " " * (cursor_pos * (w + 1)) + "^"
    return "\n".join([pos_row, map_row, cur_row])


def _perm_editor_pt(n_apm: int, existing: Optional[List[int]] = None) -> Optional[List[int]]:
    """Cursor-based permutation editor using prompt_toolkit Application.

    Arrow keys move the cursor. Type the target index + Enter to assign.
    Backspace clears the current position. Enter with no input confirms
    when the permutation is complete and non-identity. Ctrl+C cancels.
    """
    from prompt_toolkit import Application
    from prompt_toolkit.buffer import Buffer
    from prompt_toolkit.document import Document
    from prompt_toolkit.key_binding import KeyBindings
    from prompt_toolkit.layout import Layout, HSplit, Window
    from prompt_toolkit.layout.controls import FormattedTextControl, BufferControl

    w = _col_width(n_apm)
    mapping: List = list(existing) if existing else [None] * n_apm
    cursor = [0]
    result: List[Optional[List[int]]] = [None]
    msg = [""]
    input_buf = Buffer(name="perm_input", multiline=False)
    kb = KeyBindings()

    def _clear_input():
        input_buf.set_document(Document(""))

    @kb.add("left")
    def _(event):
        cursor[0] = max(0, cursor[0] - 1)
        _clear_input()
        msg[0] = ""

    @kb.add("right")
    def _(event):
        cursor[0] = min(n_apm - 1, cursor[0] + 1)
        _clear_input()
        msg[0] = ""

    @kb.add("enter")
    def _(event):
        text = input_buf.text.strip()
        if not text:
            if None in mapping:
                msg[0] = "! incomplete - fill all positions first"
                return
            if mapping == list(range(n_apm)):
                msg[0] = "! identity - swap at least one pair"
                return
            result[0] = mapping.copy()
            event.app.exit()
            return
        try:
            val = int(text)
        except ValueError:
            msg[0] = f"! expected integer, got {text!r}"
            _clear_input()
            return
        if not (0 <= val < n_apm):
            msg[0] = f"! {val} out of range 0..{n_apm - 1}"
            _clear_input()
            return
        used = {mapping[i] for i in range(n_apm) if i != cursor[0] and mapping[i] is not None}
        if val in used:
            msg[0] = f"! {val} already assigned at another position"
            _clear_input()
            return
        mapping[cursor[0]] = val
        msg[0] = ""
        _clear_input()
        cursor[0] = min(n_apm - 1, cursor[0] + 1)

    @kb.add("backspace")
    def _(event):
        if input_buf.text:
            input_buf.delete_before_cursor()
        else:
            mapping[cursor[0]] = None
            msg[0] = ""

    @kb.add("c-c")
    @kb.add("c-d")
    def _(event):
        event.app.exit()

    def get_display() -> str:
        used_elsewhere = {
            mapping[i] for i in range(n_apm) if i != cursor[0] and mapping[i] is not None
        }
        avail = sorted(set(range(n_apm)) - used_elsewhere)
        lines = [
            "",
            _render_mapping(n_apm, mapping, cursor[0]),
            f"  Available: {avail}",
        ]
        if msg[0]:
            lines.append(f"  {msg[0]}")
        elif None not in mapping:
            if mapping == list(range(n_apm)):
                lines.append("  [identity - swap at least one pair]")
            else:
                lines.append("  [complete - press Enter (empty) to confirm]")
        else:
            lines.append(
                "  [<- -> move cursor   index + Enter to set   Backspace to clear   Ctrl+C cancel]"
            )
        lines.append(f"\n  Set pos {cursor[0]} to: ")
        return "\n".join(lines)

    layout = Layout(
        HSplit([
            Window(content=FormattedTextControl(get_display)),
            Window(content=BufferControl(buffer=input_buf), height=1),
        ])
    )
    app = Application(layout=layout, key_bindings=kb, full_screen=False, mouse_support=False)
    app.run()
    return result[0]


def _perm_editor_simple(n_apm: int, existing: Optional[List[int]] = None) -> Optional[List[int]]:
    """Sequential permutation editor using plain input() (fallback).

    Navigate with ``!back`` and ``!reset``; abort with ``!cancel``.
    """
    mapping: List = list(existing) if existing else [None] * n_apm
    pos = 0
    while 0 <= pos < n_apm:
        used = {mapping[i] for i in range(n_apm) if i != pos and mapping[i] is not None}
        avail = sorted(set(range(n_apm)) - used)
        print(_render_mapping(n_apm, mapping, pos))
        print(f"  Available: {avail}")
        try:
            raw = input(f"  pos {pos}  (!back  !reset  !cancel) >> ").strip()
        except (EOFError, KeyboardInterrupt):
            raise _Abort()
        low = raw.lower()
        if low == "!cancel":
            raise _Abort()
        if low == "!back":
            if pos > 0:
                mapping[pos - 1] = None
                pos -= 1
            continue
        if low == "!reset":
            mapping = [None] * n_apm
            pos = 0
            continue
        try:
            val = int(raw)
        except ValueError:
            print("  ! expected integer")
            continue
        if not (0 <= val < n_apm):
            print(f"  ! out of range 0..{n_apm - 1}")
            continue
        if val in {mapping[i] for i in range(n_apm) if i != pos and mapping[i] is not None}:
            print(f"  ! {val} already assigned at another position")
            continue
        mapping[pos] = val
        pos += 1

    if None in mapping or mapping == list(range(n_apm)):
        return None
    return mapping


def _perm_editor(n_apm: int, existing: Optional[List[int]] = None) -> Optional[List[int]]:
    """Open the permutation editor (prompt_toolkit if available, else simple)."""
    if _PT:
        return _perm_editor_pt(n_apm, existing)
    return _perm_editor_simple(n_apm, existing)


# ── sym session (list of permutations) ───────────────────────────────────────

def _sym_session(
    n_apm: int, existing: Optional[List[List[int]]] = None
) -> Optional[List[List[int]]]:
    """Interactively manage a list of permutations for sym.

    Returns the list (possibly empty -> None) or raises ``_Abort`` to go back.
    """
    sym: List[List[int]] = list(existing) if existing else []

    def show():
        print(_SEP)
        print(f"  Sym editor  (n_apm={n_apm})")
        print(f"  Identity {list(range(n_apm))} is implicit - do not add it.")
        if sym:
            for i, p in enumerate(sym):
                print(f"    {i}: {p}")
        else:
            print("  (no permutations - only identity applies)")
        print(_SEP)
        print("  add / remove N / list / done / !back")

    show()
    while True:
        try:
            raw = input("  sym >> ").strip()
        except (EOFError, KeyboardInterrupt):
            raise _Abort()
        if not raw:
            continue
        parts = raw.split()
        cmd = parts[0].lower()
        args = parts[1:]

        if cmd in ("done", "d"):
            return sym or None

        elif cmd in ("!back", "back", "b"):
            raise _Abort()

        elif cmd in ("add", "a"):
            print("  --- new permutation ---")
            perm = _perm_editor(n_apm)
            if perm is not None:
                if perm in sym:
                    print("  ! permutation already in list")
                else:
                    sym.append(perm)
                    show()
            else:
                print("  (cancelled or identity - not added)")

        elif cmd in ("remove", "r"):
            if len(args) != 1 or not args[0].isdigit():
                print("  ! usage: remove N")
                continue
            idx = int(args[0])
            if not (0 <= idx < len(sym)):
                print(f"  ! index {idx} out of range")
                continue
            sym.pop(idx)
            show()

        elif cmd in ("list", "l"):
            show()

        else:
            print(f"  ! unknown: {cmd!r}. Commands: add / remove N / list / done / !back")


# ── species wizard ────────────────────────────────────────────────────────────

def _wizard_species(
    prompt_fn: _PromptFn, existing: Optional[Dict] = None
) -> Dict:
    """Step-through wizard for one species.

    Steps: n_apm -> n_mol -> name -> sym.
    ``!back`` goes back one step; Ctrl+C / ``!cancel`` raises ``_Abort``.
    From the sym step, ``!back`` (or aborting the sym session) returns to name.
    """
    sp: Dict = copy.deepcopy(existing) if existing else {
        "n_apm": 1, "n_mol": 1, "name": None, "sym": None
    }

    STEPS = ["n_apm", "n_mol", "name", "sym"]
    step = 0

    while step < len(STEPS):
        field = STEPS[step]

        # ── sym step ──────────────────────────────────────────────────────────
        if field == "sym":
            print(f"  --- sym (n_apm={sp['n_apm']}) ---")
            try:
                sp["sym"] = _sym_session(sp["n_apm"], existing=sp.get("sym"))
                step += 1
            except _Abort:
                step -= 1  # back to name
            continue

        # ── text / int steps ─────────────────────────────────────────────────
        back_hint = "  [!back]" if step > 0 else ""
        cur_val = sp[field]
        default_str = str(cur_val) if cur_val is not None else ""

        try:
            raw = prompt_fn(
                f"    {field}{back_hint}  [!cancel] >> ", default=default_str
            )
        except _Abort:
            raise

        low = raw.lower()
        if low == "!cancel":
            raise _Abort()
        if low == "!back" and step > 0:
            step -= 1
            continue

        # ── validate by field ─────────────────────────────────────────────────
        if field == "n_apm":
            try:
                val = int(raw) if raw else cur_val
                if val < 1:
                    print("    ! n_apm must be >= 1")
                    continue
                sp["n_apm"] = val
                # reset sym if n_apm changed
                if existing and val != existing.get("n_apm"):
                    sp["sym"] = None
                step += 1
            except (ValueError, TypeError):
                print("    ! expected integer")

        elif field == "n_mol":
            try:
                val = int(raw) if raw else cur_val
                if val < 1:
                    print("    ! n_mol must be >= 1")
                    continue
                sp["n_mol"] = val
                step += 1
            except (ValueError, TypeError):
                print("    ! expected integer")

        elif field == "name":
            # blank -> clear; '-' -> explicit clear; any text -> set
            if raw == "-":
                sp["name"] = None
            else:
                sp["name"] = raw if raw else cur_val
            step += 1

    return sp


# ---------------------------------------------------------------------------
# TopologyBuilder
# ---------------------------------------------------------------------------

class TopologyBuilder:
    """Mutable builder for :class:`~mobbrmsd.dataclass.molecular_system`.

    ``TopologyBuilder`` is a pure data-editing layer with no dependency on
    the ``mobbrmsd`` dataclasses.  Conversion to/from ``molecular_system``
    is handled by :mod:`mobbrmsd.run.parse`:

    * :func:`~mobbrmsd.run.parse.to_molecular_system`
    * :func:`~mobbrmsd.run.parse.builder_from_topology`
    * :func:`~mobbrmsd.run.parse.builder_from_mdtraj`

    Examples
    --------
    Programmatic API::

        mols = (
            TopologyBuilder()
            .add(n_apm=3, n_mol=4, name="WAT")
            .add(n_apm=2, n_mol=1, name="HF")
            .sort()
            .build()
        )

    Interactive::

        from mobbrmsd.run.parse import builder_from_topology
        mols = builder_from_topology(top).run()
    """

    def __init__(self) -> None:
        self._species: List[Dict] = []
        self._name: Optional[str] = None
        self._atom_order: Optional[List[int]] = None

    # ------------------------------------------------------------------
    # constructor (duck-typed - no dataclass import)
    # ------------------------------------------------------------------

    @classmethod
    def from_molecular_system(cls, mols) -> "TopologyBuilder":
        """Initialise from any object with ``.mols``, ``.name``, ``.atom_order``."""
        b = cls()
        b._name = mols.name
        b._atom_order = (
            mols.atom_order.tolist() if mols.atom_order is not None else None
        )
        for m in mols.mols:
            b._species.append({
                "n_apm": m.n_apm,
                "n_mol": m.n_mol,
                "sym": copy.deepcopy(m.sym),
                "name": m.name,
            })
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
        """Update fields of species at *index* (n_apm, n_mol, sym, name)."""
        sp = self._species[index]
        for k in ("n_apm", "n_mol", "sym", "name"):
            if k in kwargs:
                sp[k] = kwargs[k]
        return self

    def reorder(self, indices: List[int]) -> "TopologyBuilder":
        """Reorder species by *indices* (a permutation of 0..n-1)."""
        if sorted(indices) != list(range(len(self._species))):
            raise ValueError("indices must be a permutation of 0 .. n_species-1")
        self._species = [self._species[i] for i in indices]
        return self

    def sort(self) -> "TopologyBuilder":
        """Structural sort: ascending n_mol, then descending n_apm."""
        indices = sorted(
            range(len(self._species)),
            key=lambda i: (self._species[i]["n_mol"], -self._species[i]["n_apm"]),
        )
        self._species = [self._species[i] for i in indices]
        return self

    def set_name(self, name: Optional[str]) -> "TopologyBuilder":
        self._name = name
        return self

    def set_atom_order(self, order: Optional[List[int]]) -> "TopologyBuilder":
        self._atom_order = list(order) if order is not None else None
        return self

    # ------------------------------------------------------------------
    # build - delegates to parse.to_molecular_system
    # ------------------------------------------------------------------

    def build(self):
        """Return the final frozen ``molecular_system``."""
        from .parse import to_molecular_system
        return to_molecular_system(self)

    # ------------------------------------------------------------------
    # display
    # ------------------------------------------------------------------

    def _show(self) -> None:
        print(_SEP)
        print(f"  System name : {self._name or '(none)'}")
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
                sym_repr = repr(sp["sym"]) if sp["sym"] else "(none)"
                print(
                    f"  {i:>3}  {sp['n_apm']:>6}  {sp['n_mol']:>6}"
                    f"  {sp['name'] or '(none)':<16}  {sym_repr}"
                )
        print(_SEP)

    # ------------------------------------------------------------------
    # interactive REPL
    # ------------------------------------------------------------------

    _HELP = """\
  Commands:
    list  (l)          - show current species
    add   (a)          - add a new species
    edit  (e) N        - edit species N
    remove(r) N        - remove species N
    move  (mv) N M     - move species N to position M
    sort  (s)          - structural sort (n_mol asc, n_apm desc)
    name  (n) [TEXT]   - set system name (empty to clear)
    order (o) [JSON]   - set atom_order from JSON int list (empty to clear)
    done  (d)          - finish and return molecular_system
    quit  (q)          - abort
    help  (h / ?)      - show this help
"""

    def run(self):
        """Start the interactive REPL.

        Returns the built ``molecular_system`` on ``done``, or ``None`` on
        ``quit``.  Requires ``prompt_toolkit`` for cursor movement and history;
        falls back to plain ``input()`` otherwise.
        """
        if not _PT:
            print("  (prompt_toolkit not found - using plain input; "
                  "install prompt-toolkit for cursor movement)")

        repl_prompt = _PromptFn(use_history=True)

        print(_SEP)
        print("  mobbrmsd topology editor  (type 'help' for commands)")
        self._show()

        while True:
            try:
                raw = repl_prompt("  >> ")
            except _Abort:
                # Ctrl+C in the main loop: clear the line and continue
                print()
                continue
            except EOFError:
                print("\n  Quit.")
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
                field_prompt = _PromptFn(use_history=False)
                try:
                    sp = _wizard_species(field_prompt)
                    self._species.append(sp)
                    self._show()
                except _Abort:
                    print("  add cancelled.")

            elif cmd in ("edit", "e"):
                if len(args) != 1 or not args[0].isdigit():
                    print("  ! usage: edit N")
                    continue
                idx = int(args[0])
                if not (0 <= idx < len(self._species)):
                    print(f"  ! index {idx} out of range")
                    continue
                field_prompt = _PromptFn(use_history=False)
                try:
                    sp = _wizard_species(field_prompt, existing=self._species[idx])
                    self._species[idx] = sp
                    self._show()
                except _Abort:
                    print("  edit cancelled.")

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
                    print(f"  ! indices out of range (0..{n - 1})")
                    continue
                sp = self._species.pop(src)
                self._species.insert(dst, sp)
                self._show()

            elif cmd in ("sort", "s"):
                self.sort()
                self._show()

            elif cmd in ("name", "n"):
                self.set_name(" ".join(args) if args else None)
                print(f"  System name: {self._name!r}")

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
                print("  Quit.")
                return None

            else:
                print(f"  ! unknown command: {cmd!r}. Type 'help'.")


# ---------------------------------------------------------------------------
# standalone entry point
# ---------------------------------------------------------------------------

def _main() -> None:
    import argparse
    import sys

    ap = argparse.ArgumentParser(
        prog="python -m mobbrmsd.run.topology_builder",
        description="Interactive molecular_system editor. Prints JSON on 'done'.",
    )
    ap.add_argument("-t", "--topology", metavar="FILE", default=None,
                    help="MDTraj topology/structure file (PDB, GRO, …).")
    ap.add_argument("-c", "--coords", metavar="FILE", default=None,
                    help="Trajectory for variance-based ordering (needs --topology).")
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
