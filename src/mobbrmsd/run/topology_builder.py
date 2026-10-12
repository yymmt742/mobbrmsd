"""Interactive builder for molecular_system.

Since ``molecular_system`` is a frozen dataclass, this module provides a
mutable ``TopologyBuilder`` that accumulates species edits interactively
and emits a ``molecular_system`` when the user is done.

``TopologyBuilder`` itself has no dependency on :mod:`mobbrmsd.dataclass`.
Conversion to/from ``molecular_system`` is handled by helpers in
:mod:`mobbrmsd.run.parse`.

``atom_names`` (per-atom string labels) are stored in the builder only and
stripped when converting to ``molecular_system``.

Requires ``prompt_toolkit`` for cursor movement; falls back to plain
``input()`` otherwise.
"""

from __future__ import annotations

import copy
import json
from typing import Dict, List, Optional

_SEP = "  " + "-" * 72

# ── prompt_toolkit detection ──────────────────────────────────────────────────

try:
    from prompt_toolkit import PromptSession as _PTSession
    from prompt_toolkit.history import InMemoryHistory as _PTHistory
    from prompt_toolkit.key_binding import KeyBindings as _PTKeyBindings
    _PT = True
except ImportError:
    _PT = False


# ── control-flow exceptions ───────────────────────────────────────────────────

class _Abort(Exception):
    """Cancel the current interactive operation; the REPL continues."""

class _Back(Exception):
    """Go back one step in the current wizard."""


# ── prompt wrappers ───────────────────────────────────────────────────────────

class _ReplPrompt:
    """Main-REPL prompt with history.  Ctrl+C raises ``_Abort``."""

    def __init__(self) -> None:
        if _PT:
            kb = _PTKeyBindings()

            @kb.add('c-c')
            def _(event):
                # store reason so we can detect it after prompt() returns
                event.app._abort_requested = True
                event.app.exit(result="")

            self._session = _PTSession(key_bindings=kb, history=_PTHistory())
        else:
            self._session = None

    def __call__(self, text: str) -> str:
        if _PT:
            try:
                self._session.app  # may not exist yet
                app_ref = [None]
                # capture app reference via pre_run
                def _pre_run():
                    app_ref[0] = self._session.app
                result = self._session.prompt(text, pre_run=_pre_run)
            except (KeyboardInterrupt, EOFError):
                raise _Abort()
            if app_ref[0] is not None and getattr(app_ref[0], '_abort_requested', False):
                app_ref[0]._abort_requested = False
                raise _Abort()
            return result.strip() if isinstance(result, str) else ""
        else:
            try:
                return input(text).strip()
            except (KeyboardInterrupt, EOFError):
                raise _Abort()


class _WizardPrompt:
    """Wizard-step prompt.

    * Up arrow  -> raise ``_Back``
    * Ctrl+C    -> raise ``_Abort``
    * Default shown as ``[val]`` in prompt text; buffer starts empty.
    """

    def __init__(self) -> None:
        if _PT:
            kb = _PTKeyBindings()

            @kb.add('up')
            def _(event):
                event.app._wiz_back = True
                event.app.exit(result="")

            @kb.add('c-c')
            @kb.add('c-d')
            def _(event):
                event.app._wiz_abort = True
                event.app.exit(result="")

            self._session = _PTSession(key_bindings=kb)
        else:
            self._session = None

    def __call__(self, label: str, default: str = "") -> str:
        """Prompt for a wizard field.  Returns the typed string (or *default* if empty)."""
        hint = f" [{default}]" if default else ""
        text = f"    {label}{hint} >> "
        if _PT:
            app_ref = [None]
            def _pre_run():
                a = self._session.app
                a._wiz_back = False
                a._wiz_abort = False
                app_ref[0] = a
            try:
                result = self._session.prompt(text, pre_run=_pre_run)
            except (KeyboardInterrupt, EOFError):
                raise _Abort()
            a = app_ref[0]
            if a and getattr(a, '_wiz_back', False):
                raise _Back()
            if a and getattr(a, '_wiz_abort', False):
                raise _Abort()
            raw = result.strip() if isinstance(result, str) else ""
            return raw if raw else default
        else:
            try:
                raw = input(text).strip()
                if raw.lower() == "!back":
                    raise _Back()
                return raw if raw else default
            except (KeyboardInterrupt, EOFError):
                raise _Abort()


# ── atom name utilities ───────────────────────────────────────────────────────

def _sanitize_atom_name(name: str) -> str:
    """Strip whitespace and truncate to 4 characters."""
    return name.strip()[:4]


def _fill_atom_names(names: List, n_apm: int) -> List[str]:
    """Replace None entries with auto-generated unique names (X001, X002, ...)."""
    result: List = list(names)
    assert len(result) == n_apm
    existing = {n for n in result if n is not None}
    counter = 1
    for i in range(n_apm):
        if result[i] is None:
            while True:
                candidate = f"X{counter:03d}"
                counter += 1
                if candidate not in existing:
                    result[i] = candidate
                    existing.add(candidate)
                    break
    return result


# ── permutation / atom-names display ─────────────────────────────────────────

def _col_width(n_apm: int, atom_names: Optional[List] = None) -> int:
    """Column width that fits both indices and atom names."""
    w = max(2, len(str(n_apm - 1)))
    if atom_names:
        w = max(w, max((len(n) for n in atom_names if n), default=0))
    return max(w, 4 if atom_names else w)


def _render_display(
    n_apm: int,
    mapping: List[int],
    cursor_pos: int,
    atom_names: Optional[List] = None,
) -> str:
    """Render the Pos / Map / (Name) / Cur rows."""
    w = _col_width(n_apm, atom_names)
    rows = []
    if atom_names:
        rows.append("  Name: " + " ".join(
            f"{(n or '????'):<{w}}" for n in atom_names
        ))
    rows.append("  Pos:  " + " ".join(f"{i:>{w}d}" for i in range(n_apm)))
    rows.append("  Map:  " + " ".join(f"{v:>{w}d}" for v in mapping))
    rows.append("  Cur:  " + " " * (cursor_pos * (w + 1)) + "^")
    return "\n".join(rows)


def _render_names_display(
    n_apm: int,
    names: List,
    cursor_pos: int,
) -> str:
    """Render Name / Pos / Cur rows for the atom_names editor."""
    w = _col_width(n_apm, [n for n in names if n])
    rows = [
        "  Name: " + " ".join(f"{(n or '????'):<{w}}" for n in names),
        "  Pos:  " + " ".join(f"{i:>{w}d}" for i in range(n_apm)),
        "  Cur:  " + " " * (cursor_pos * (w + 1)) + "^",
    ]
    return "\n".join(rows)


# ── permutation editor ────────────────────────────────────────────────────────

def _perm_editor_pt(
    n_apm: int,
    existing: Optional[List[int]] = None,
    atom_names: Optional[List] = None,
) -> Optional[List[int]]:
    """Cursor-based permutation editor (prompt_toolkit Application).

    Starts from identity (or *existing*).  Left/Right moves cursor.
    Type a target index + Enter to swap ``mapping[cursor]`` with
    ``mapping[target]``.  Enter on empty buffer confirms.  Ctrl+C cancels.
    """
    from prompt_toolkit import Application
    from prompt_toolkit.buffer import Buffer
    from prompt_toolkit.document import Document
    from prompt_toolkit.key_binding import KeyBindings
    from prompt_toolkit.layout import Layout, HSplit, Window
    from prompt_toolkit.layout.controls import FormattedTextControl, BufferControl

    mapping: List[int] = list(existing) if existing else list(range(n_apm))
    cursor = [0]
    result: List[Optional[List[int]]] = [None]
    msg = [""]
    input_buf = Buffer(name="perm_input", multiline=False)
    kb = KeyBindings()

    def _clear() -> None:
        input_buf.set_document(Document(""))

    @kb.add("left")
    def _(event):
        cursor[0] = max(0, cursor[0] - 1)
        _clear()
        msg[0] = ""

    @kb.add("right")
    def _(event):
        cursor[0] = min(n_apm - 1, cursor[0] + 1)
        _clear()
        msg[0] = ""

    @kb.add("enter")
    def _(event):
        text = input_buf.text.strip()
        if not text:
            result[0] = mapping.copy()
            event.app.exit()
            return
        try:
            target = int(text)
        except ValueError:
            msg[0] = f"! expected integer, got {text!r}"
            _clear()
            return
        if not (0 <= target < n_apm):
            msg[0] = f"! {target} out of range 0..{n_apm - 1}"
            _clear()
            return
        if target == cursor[0]:
            msg[0] = "! same position - nothing to swap"
            _clear()
            return
        mapping[cursor[0]], mapping[target] = mapping[target], mapping[cursor[0]]
        msg[0] = f"swapped pos {cursor[0]} <-> pos {target}"
        _clear()

    @kb.add("c-c")
    @kb.add("c-d")
    def _(event):
        event.app.exit()

    def get_display() -> str:
        lines = [
            "",
            _render_display(n_apm, mapping, cursor[0], atom_names),
        ]
        if msg[0]:
            lines.append(f"  {msg[0]}")
        elif mapping == list(range(n_apm)):
            lines.append("  [identity - swap pairs, or Enter to accept as-is]")
        else:
            lines.append("  [Enter confirm   <- -> move   target+Enter swap   Ctrl+C cancel]")
        lines.append(f"\n  Swap pos {cursor[0]} with: ")
        return "\n".join(lines)

    layout = Layout(
        HSplit([
            Window(content=FormattedTextControl(get_display)),
            Window(content=BufferControl(buffer=input_buf), height=1),
        ])
    )
    Application(layout=layout, key_bindings=kb, full_screen=False, mouse_support=False).run()
    return result[0]


def _perm_editor_simple(
    n_apm: int,
    existing: Optional[List[int]] = None,
    atom_names: Optional[List] = None,
) -> Optional[List[int]]:
    """Sequential permutation editor using plain input() (fallback)."""
    mapping: List[int] = list(existing) if existing else list(range(n_apm))
    cursor = 0

    while True:
        print(_render_display(n_apm, mapping, cursor, atom_names))
        print("  (n/p move   target+Enter swap   blank=confirm   !reset   !cancel)")
        try:
            raw = input(f"  Swap pos {cursor} with: ").strip()
        except (EOFError, KeyboardInterrupt):
            raise _Abort()

        low = raw.lower()
        if low in ("!cancel",):
            raise _Abort()
        if low == "!reset":
            mapping = list(range(n_apm))
            continue
        if low in ("n", "next"):
            cursor = min(n_apm - 1, cursor + 1)
            continue
        if low in ("p", "prev"):
            cursor = max(0, cursor - 1)
            continue
        if not raw:
            return mapping.copy()
        try:
            target = int(raw)
        except ValueError:
            print("  ! expected integer")
            continue
        if not (0 <= target < n_apm):
            print(f"  ! out of range 0..{n_apm - 1}")
            continue
        if target == cursor:
            print("  ! same position")
            continue
        mapping[cursor], mapping[target] = mapping[target], mapping[cursor]
        print(f"  swapped pos {cursor} <-> pos {target}")


def _perm_editor(
    n_apm: int,
    existing: Optional[List[int]] = None,
    atom_names: Optional[List] = None,
) -> Optional[List[int]]:
    if _PT:
        return _perm_editor_pt(n_apm, existing, atom_names)
    return _perm_editor_simple(n_apm, existing, atom_names)


# ── atom_names editor ─────────────────────────────────────────────────────────

def _atom_names_editor_pt(
    n_apm: int, existing: Optional[List] = None
) -> Optional[List[str]]:
    """Cursor-based atom_names editor (prompt_toolkit Application).

    Left/Right moves cursor.  Type a name + Enter to set the current position.
    Empty Enter auto-assigns a name for that position.
    Ctrl+D confirms (auto-fills remaining unset names).
    Up arrow raises ``_Back``.  Ctrl+C raises ``_Abort`` (skip this step).
    """
    from prompt_toolkit import Application
    from prompt_toolkit.buffer import Buffer
    from prompt_toolkit.document import Document
    from prompt_toolkit.key_binding import KeyBindings
    from prompt_toolkit.layout import Layout, HSplit, Window
    from prompt_toolkit.layout.controls import FormattedTextControl, BufferControl

    names: List = list(existing) if existing else [None] * n_apm
    assert len(names) == n_apm

    cursor = [0]
    result: List[Optional[List[str]]] = [None]
    exit_reason: List[str] = [""]
    msg = [""]
    input_buf = Buffer(name="name_input", multiline=False)
    kb = KeyBindings()

    def _clear() -> None:
        input_buf.set_document(Document(""))

    @kb.add("left")
    def _(event):
        cursor[0] = max(0, cursor[0] - 1)
        _clear()
        msg[0] = ""

    @kb.add("right")
    def _(event):
        cursor[0] = min(n_apm - 1, cursor[0] + 1)
        _clear()
        msg[0] = ""

    @kb.add("enter")
    def _(event):
        text = input_buf.text.strip()
        if not text:
            # auto-assign for current position
            tmp = list(names)
            if tmp[cursor[0]] is None:
                tmp[cursor[0]] = "__PLACEHOLDER__"
            filled = _fill_atom_names(tmp, n_apm)
            names[cursor[0]] = filled[cursor[0]]
            msg[0] = f"auto: {names[cursor[0]]}"
        else:
            sanitized = _sanitize_atom_name(text)
            if not sanitized:
                msg[0] = "! empty after sanitize"
                _clear()
                return
            used = [names[i] for i in range(n_apm) if i != cursor[0]]
            if sanitized in used:
                msg[0] = f"! '{sanitized}' already used"
                _clear()
                return
            names[cursor[0]] = sanitized
            msg[0] = f"set: {sanitized}"
        _clear()
        cursor[0] = min(n_apm - 1, cursor[0] + 1)

    @kb.add("c-d")
    def _(event):
        result[0] = _fill_atom_names(list(names), n_apm)
        exit_reason[0] = "confirm"
        event.app.exit()

    @kb.add("up")
    def _(event):
        exit_reason[0] = "back"
        event.app.exit()

    @kb.add("c-c")
    def _(event):
        exit_reason[0] = "abort"
        event.app.exit()

    def get_display() -> str:
        lines = [
            "",
            _render_names_display(n_apm, names, cursor[0]),
        ]
        if msg[0]:
            lines.append(f"  {msg[0]}")
        else:
            unset = sum(1 for n in names if n is None)
            if unset:
                lines.append(f"  [{unset} unset - Enter to auto-assign   Ctrl+D confirm all   Up=back   Ctrl+C skip]")
            else:
                lines.append("  [all set - Ctrl+D to confirm   Up=back   Ctrl+C skip]")
        lines.append(f"\n  Name for pos {cursor[0]}: ")
        return "\n".join(lines)

    layout = Layout(
        HSplit([
            Window(content=FormattedTextControl(get_display)),
            Window(content=BufferControl(buffer=input_buf), height=1),
        ])
    )
    Application(layout=layout, key_bindings=kb, full_screen=False, mouse_support=False).run()

    if exit_reason[0] == "back":
        raise _Back()
    if exit_reason[0] == "abort":
        raise _Abort()
    return result[0]  # None if exited without Ctrl+D


def _atom_names_editor_simple(
    n_apm: int, existing: Optional[List] = None
) -> Optional[List[str]]:
    """Sequential atom_names editor using plain input() (fallback).

    Commands: ``n``/``p`` move; name+Enter set; blank=auto; ``!done`` confirm;
    ``!back`` go back; Ctrl+C skip.
    """
    names: List = list(existing) if existing else [None] * n_apm
    assert len(names) == n_apm
    cursor = 0

    while True:
        print(_render_names_display(n_apm, names, cursor))
        unset = sum(1 for n in names if n is None)
        print(f"  (n/p move   name+Enter set   blank=auto   !done confirm   !back   Ctrl+C skip)")
        try:
            raw = input(f"  Name for pos {cursor}: ").strip()
        except (EOFError, KeyboardInterrupt):
            raise _Abort()

        low = raw.lower()
        if low in ("!back",):
            raise _Back()
        if low in ("!done", "!d"):
            return _fill_atom_names(list(names), n_apm)
        if low in ("n", "next"):
            cursor = min(n_apm - 1, cursor + 1)
            continue
        if low in ("p", "prev"):
            cursor = max(0, cursor - 1)
            continue
        if not raw:
            tmp = list(names)
            if tmp[cursor] is None:
                tmp[cursor] = "__PLACEHOLDER__"
            filled = _fill_atom_names(tmp, n_apm)
            names[cursor] = filled[cursor]
            print(f"  auto: {names[cursor]}")
            cursor = min(n_apm - 1, cursor + 1)
            continue
        sanitized = _sanitize_atom_name(raw)
        if not sanitized:
            print("  ! empty after sanitize")
            continue
        used = [names[i] for i in range(n_apm) if i != cursor]
        if sanitized in used:
            print(f"  ! '{sanitized}' already used")
            continue
        names[cursor] = sanitized
        cursor = min(n_apm - 1, cursor + 1)


def _atom_names_editor(
    n_apm: int, existing: Optional[List] = None
) -> Optional[List[str]]:
    """Open the atom_names editor.  Returns confirmed names or raises _Back/_Abort."""
    if _PT:
        return _atom_names_editor_pt(n_apm, existing)
    return _atom_names_editor_simple(n_apm, existing)


# ── sym session (list of permutations) ───────────────────────────────────────

def _sym_session(
    n_apm: int,
    existing: Optional[List[List[int]]] = None,
    atom_names: Optional[List] = None,
) -> Optional[List[List[int]]]:
    """Interactively manage a list of permutations for sym.

    Returns the list (possibly empty -> None) or raises ``_Abort`` to go back.
    """
    sym: List[List[int]] = list(existing) if existing else []

    def show() -> None:
        print(_SEP)
        print(f"  Sym editor  (n_apm={n_apm})")
        if atom_names:
            print("  Names: " + "  ".join(f"{n}" for n in atom_names))
        print(f"  Identity {list(range(n_apm))} is always implicit.")
        if sym:
            for i, p in enumerate(sym):
                print(f"    {i}: {p}")
        else:
            print("  (no permutations)")
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
            print("  --- edit permutation (identity shown; Enter to confirm) ---")
            perm = _perm_editor(n_apm, atom_names=atom_names)
            if perm is None:
                print("  (cancelled)")
                continue
            if perm == list(range(n_apm)):
                print("  ! identity permutation - not added")
            elif perm in sym:
                print("  ! duplicate - not added")
            else:
                sym.append(perm)
                show()

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
    wiz: _WizardPrompt, existing: Optional[Dict] = None
) -> Dict:
    """Step-through wizard for one species.

    Steps: n_apm -> n_mol -> name -> atom_names -> sym.
    Up arrow goes back one step; Ctrl+C raises ``_Abort``.
    In the atom_names step: Ctrl+C skips (keeps current atom_names), Up goes back.
    """
    sp: Dict = copy.deepcopy(existing) if existing else {
        "n_apm": 1, "n_mol": 1, "name": None, "atom_names": None, "sym": None
    }
    if "atom_names" not in sp:
        sp["atom_names"] = None

    STEPS = ["n_apm", "n_mol", "name", "atom_names", "sym"]
    step = 0

    while step < len(STEPS):
        field = STEPS[step]

        # ── atom_names step ───────────────────────────────────────────────────
        if field == "atom_names":
            n = sp["n_apm"]
            existing_names = sp.get("atom_names") or [None] * n
            if len(existing_names) != n:
                existing_names = [None] * n
            print(f"  --- atom names (n_apm={n}, Ctrl+D confirm, Up=back, Ctrl+C skip) ---")
            try:
                new_names = _atom_names_editor(n, existing_names)
                if new_names is not None:
                    sp["atom_names"] = new_names
                step += 1  # proceed to sym
            except _Back:
                step -= 1  # back to name
            except _Abort:
                step += 1  # Ctrl+C = skip atom_names, proceed to sym
            continue

        # ── sym step ──────────────────────────────────────────────────────────
        if field == "sym":
            print(f"  --- sym (n_apm={sp['n_apm']}) ---")
            try:
                sp["sym"] = _sym_session(
                    sp["n_apm"],
                    existing=sp.get("sym"),
                    atom_names=sp.get("atom_names"),
                )
                step += 1
            except _Abort:
                step -= 1  # back to atom_names
            continue

        # ── text / int steps ─────────────────────────────────────────────────
        cur = sp[field]
        default_str = str(cur) if cur is not None else ""
        try:
            raw = wiz(field, default=default_str)
        except _Back:
            if step > 0:
                step -= 1
            continue
        except _Abort:
            raise

        # ── validate by field ─────────────────────────────────────────────────
        if field == "n_apm":
            try:
                val = int(raw) if raw else (cur or 1)
                if val < 1:
                    print("    ! n_apm must be >= 1")
                    continue
                if existing and val != existing.get("n_apm"):
                    sp["atom_names"] = None  # reset when n_apm changes
                    sp["sym"] = None
                sp["n_apm"] = val
                step += 1
            except (ValueError, TypeError):
                print("    ! expected integer")

        elif field == "n_mol":
            try:
                val = int(raw) if raw else (cur or 1)
                if val < 1:
                    print("    ! n_mol must be >= 1")
                    continue
                sp["n_mol"] = val
                step += 1
            except (ValueError, TypeError):
                print("    ! expected integer")

        elif field == "name":
            # blank keeps current; "-" clears
            if raw == "-":
                sp["name"] = None
            else:
                sp["name"] = raw if raw else cur
            step += 1

    return sp


# ---------------------------------------------------------------------------
# TopologyBuilder
# ---------------------------------------------------------------------------

class TopologyBuilder:
    """Mutable builder for :class:`~mobbrmsd.dataclass.molecular_system`.

    ``TopologyBuilder`` is a pure data-editing layer with no dependency on
    the ``mobbrmsd`` dataclasses.  ``atom_names`` per species are held here
    and stripped when converting to ``molecular_system``.

    Examples
    --------
    Programmatic API::

        mols = (
            TopologyBuilder()
            .add(n_apm=3, n_mol=4, name="WAT",
                 atom_names=["OW", "HW1", "HW2"])
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
    # constructor (duck-typed -- no dataclass import)
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
                "atom_names": None,
            })
        return b

    @classmethod
    def from_dict(cls, d: dict) -> "TopologyBuilder":
        """Reconstruct a builder from a dict (builder JSON format, including atom_names)."""
        b = cls()
        ms = d.get("molecular_system", d)
        if isinstance(ms, str):
            ms = json.loads(ms)
        b._name = ms.get("name")
        b._atom_order = ms.get("atom_order")
        raw_mols = ms.get("molecules", [])
        if isinstance(raw_mols, dict):
            raw_mols = [raw_mols]
        for m in raw_mols:
            b._species.append({
                "n_apm": m["n_apm"],
                "n_mol": m["n_mol"],
                "sym": m.get("sym"),
                "name": m.get("name"),
                "atom_names": m.get("atom_names"),
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
        atom_names: Optional[List[str]] = None,
    ) -> "TopologyBuilder":
        """Append a new species."""
        self._species.append({
            "n_apm": n_apm,
            "n_mol": n_mol,
            "sym": sym,
            "name": name,
            "atom_names": atom_names,
        })
        return self

    def remove(self, index: int) -> "TopologyBuilder":
        self._species.pop(index)
        return self

    def edit(self, index: int, **kwargs) -> "TopologyBuilder":
        sp = self._species[index]
        for k in ("n_apm", "n_mol", "sym", "name", "atom_names"):
            if k in kwargs:
                sp[k] = kwargs[k]
        return self

    def reorder(self, indices: List[int]) -> "TopologyBuilder":
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
    # build -- delegates to parse.to_molecular_system (strips atom_names)
    # ------------------------------------------------------------------

    def build(self):
        """Return the final frozen ``molecular_system`` (atom_names excluded)."""
        from .parse import to_molecular_system
        return to_molecular_system(self)

    # ------------------------------------------------------------------
    # JSON serialization (includes atom_names)
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """Serialize builder state to dict.  Includes ``atom_names`` if present."""
        mol_list = []
        for sp in self._species:
            entry: dict = {"n_apm": sp["n_apm"], "n_mol": sp["n_mol"]}
            if sp.get("sym"):
                entry["sym"] = sp["sym"]
            if sp.get("name"):
                entry["name"] = sp["name"]
            if sp.get("atom_names"):
                entry["atom_names"] = sp["atom_names"]
            mol_list.append(entry)
        inner: dict = {"molecules": mol_list}
        if self._name:
            inner["name"] = self._name
        if self._atom_order is not None:
            inner["atom_order"] = self._atom_order
        return {"molecular_system": inner}

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
            print(f"  {'#':>3}  {'n_apm':>6}  {'n_mol':>6}  {'name':<12}  {'atom_names':<20}  sym")
            for i, sp in enumerate(self._species):
                sym_repr = repr(sp["sym"]) if sp.get("sym") else "(none)"
                an = sp.get("atom_names")
                an_repr = repr(an) if an else "(none)"
                if len(an_repr) > 20:
                    an_repr = an_repr[:17] + "..."
                print(
                    f"  {i:>3}  {sp['n_apm']:>6}  {sp['n_mol']:>6}"
                    f"  {sp.get('name') or '(none)':<12}  {an_repr:<20}  {sym_repr}"
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
        ``quit``.
        """
        if not _PT:
            print("  (prompt_toolkit not found - using plain input)")

        repl_prompt = _ReplPrompt()
        wiz = _WizardPrompt()

        print(_SEP)
        print("  mobbrmsd topology editor  (type 'help' for commands)")
        self._show()

        while True:
            try:
                raw = repl_prompt("  >> ")
            except _Abort:
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
                try:
                    sp = _wizard_species(wiz)
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
                try:
                    sp = _wizard_species(wiz, existing=self._species[idx])
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
                    help="MDTraj topology/structure file (PDB, GRO, ...).")
    ap.add_argument("-c", "--coords", metavar="FILE", default=None,
                    help="Trajectory for variance-based ordering (needs --topology).")
    ap.add_argument("-j", "--json", metavar="FILE", default=None,
                    help="Load builder state from a JSON file (supports atom_names).")
    ap.add_argument("-n", "--name", metavar="NAME", default=None,
                    help="Initial system name.")
    ap.add_argument("-o", "--output", metavar="FILE", default=None,
                    help="Write JSON to FILE instead of stdout.")
    ap.add_argument("--builder-json", action="store_true",
                    help="Output builder JSON (includes atom_names) instead of molecular_system JSON.")
    args = ap.parse_args()

    builder = None

    if args.json is not None:
        try:
            with open(args.json) as f:
                builder = TopologyBuilder.from_dict(json.load(f))
        except Exception as e:
            print(f"ERROR loading JSON: {e}", file=sys.stderr)
            sys.exit(1)

    elif args.topology is not None:
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

    if builder is None:
        builder = TopologyBuilder()
        if args.name:
            builder.set_name(args.name)

    ms = builder.run()
    if ms is None:
        sys.exit(1)

    if args.builder_json:
        output = json.dumps(builder.to_dict(), indent=2)
    else:
        output = json.dumps(ms.to_dict(), indent=2)

    if args.output is not None:
        with open(args.output, "w") as f:
            f.write(output)
            f.write("\n")
        print(f"  Written to {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    _main()
