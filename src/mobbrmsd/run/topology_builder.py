"""Interactive builder for molecular_system."""

from __future__ import annotations

import copy
import json
import shutil
import sys
from typing import Dict, List, Optional

_SEP = "  " + "-" * 72

# ── ANSI helpers ──────────────────────────────────────────────────────────────

def _ansi_up_clear(n: int) -> None:
    if n > 0:
        sys.stdout.write(f"\033[{n}A\033[J")
    else:
        sys.stdout.write("\033[J")
    sys.stdout.flush()


def _cursor_save() -> None:
    sys.stdout.write("\033[s")
    sys.stdout.flush()


def _cursor_restore_clear() -> None:
    sys.stdout.write("\033[u\033[J")
    sys.stdout.flush()


def _trunc(s: str, width: int = 80) -> str:
    """Truncate *s* to *width* chars, appending '...' if shortened."""
    return s if len(s) <= width else s[:width - 3] + "..."


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
    pass

class _Back(Exception):
    pass


# ── sym session REPL prompt ───────────────────────────────────────────────────

class _SymPrompt:
    """sym-session REPL prompt.

    * Up arrow / empty Enter -> raises _Back / returns ""
    * Ctrl+C                 -> raises _Abort
    """

    def __init__(self) -> None:
        if _PT:
            kb = _PTKeyBindings()

            @kb.add('up')
            def _(event):
                event.app._sym_back = True
                event.app.exit(result="")

            @kb.add('c-c')
            def _(event):
                event.app._sym_abort = True
                event.app.exit(result="")

            self._session = _PTSession(key_bindings=kb)

    def __call__(self, text: str) -> str:
        if _PT:
            app_ref = [None]

            def _pre_run():
                a = self._session.app
                a._sym_back = False
                a._sym_abort = False
                app_ref[0] = a

            try:
                result = self._session.prompt(text, pre_run=_pre_run)
            except (KeyboardInterrupt, EOFError):
                raise _Abort()
            a = app_ref[0]
            if a and getattr(a, '_sym_back', False):
                raise _Back()
            if a and getattr(a, '_sym_abort', False):
                raise _Abort()
            return result.strip() if isinstance(result, str) else ""
        else:
            try:
                raw = input(text).strip()
                if raw.lower() in ("!back", "b"):
                    raise _Back()
                return raw
            except (KeyboardInterrupt, EOFError):
                raise _Abort()


# ── prompt wrappers ───────────────────────────────────────────────────────────

class _ReplPrompt:
    def __init__(self) -> None:
        if _PT:
            kb = _PTKeyBindings()

            @kb.add('c-c')
            def _(event):
                event.app._abort_requested = True
                event.app.exit(result="")

            self._session = _PTSession(key_bindings=kb, history=_PTHistory())

    def __call__(self, text: str) -> str:
        if _PT:
            app_ref = [None]

            def _pre_run():
                app_ref[0] = self._session.app

            try:
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
    """Up=back, Ctrl+C=abort, default shown as [val] in prompt, buffer empty."""

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

    def __call__(self, label: str, default: str = "") -> str:
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
    return name.strip()[:4]


def _fill_atom_names(names: List, n_apm: int) -> List[str]:
    """Fill None entries with unique auto-generated names (X001, X002, ...)."""
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


# ── permutation display ───────────────────────────────────────────────────────

def _col_width(n_apm: int, atom_names: Optional[List] = None) -> int:
    w = max(2, len(str(n_apm - 1)))
    if atom_names:
        w = max(w, max((len(n) for n in atom_names if n), default=0))
    return max(w, 4 if atom_names else w)


def _perm_pages(n_apm: int, atom_names: Optional[List] = None) -> tuple:
    """Return (items_per_page, n_pages) for the permutation display."""
    term_w = shutil.get_terminal_size((80, 24)).columns
    w = _col_width(n_apm, atom_names)
    indicator_w = len(f" ({n_apm}/{n_apm})")
    items = max(1, (term_w - 8 - indicator_w) // (w + 1))
    return items, max(1, (n_apm + items - 1) // items)


def _render_display(
    n_apm: int,
    mapping: List[int],
    cursor_pos: int,
    atom_names: Optional[List] = None,
    page: int = 0,
    items_per_page: Optional[int] = None,
    n_pages: int = 1,
) -> str:
    """Render one page of Name / Map / cursor rows.  No Pos row, no 'Cur:' label.

    Appends ``(n/N)`` page indicator when n_pages > 1.
    """
    w = _col_width(n_apm, atom_names)
    if items_per_page is None:
        items_per_page = n_apm
    page_start = page * items_per_page
    page_end = min(page_start + items_per_page, n_apm)
    indicator = f" ({page + 1}/{n_pages})" if n_pages > 1 else ""

    rows = []
    if atom_names:
        rows.append("  Name: " + " ".join(
            f"{(atom_names[j] or '????'):<{w}}" for j in range(page_start, page_end)
        ) + indicator)
    rows.append("  Map:  " + " ".join(
        f"{mapping[j]:>{w}d}" for j in range(page_start, page_end)
    ) + indicator)

    if page_start <= cursor_pos < page_end:
        local = cursor_pos - page_start
        rows.append("        " + " " * (local * (w + 1)) + "^")
    elif n_pages > 1:
        rows.append(f"        (cursor at pos {cursor_pos})")

    return "\n".join(rows)


# ── atom_names display (wrapped, no Pos row, no Cur label) ───────────────────

def _render_names_wrapped(n_apm: int, names: List, cursor_pos: int) -> str:
    """Name row wrapped to terminal width, cursor indicator below the right row."""
    preview = _fill_atom_names(list(names), n_apm)
    term_w = shutil.get_terminal_size((80, 24)).columns
    col_w = max(4, _col_width(n_apm, [n for n in preview if n]))
    prefix = "  Name: "
    cont   = " " * len(prefix)
    items_per_row = max(1, (term_w - len(prefix)) // (col_w + 1))

    lines: List[str] = []
    i = 0
    row_idx = 0
    while i < n_apm:
        chunk_end = min(i + items_per_row, n_apm)
        label = prefix if row_idx == 0 else cont
        lines.append(label + " ".join(f"{preview[j]:<{col_w}}" for j in range(i, chunk_end)))
        if i <= cursor_pos < chunk_end:
            local = cursor_pos - i
            lines.append(cont + " " * (local * (col_w + 1)) + "^")
        i = chunk_end
        row_idx += 1

    return "\n".join(lines)


def _wrap_perm_line(perm: List[int], label: str, suffix: str = "") -> List[str]:
    """Format a permutation with line wrapping at terminal width."""
    term_w = shutil.get_terminal_size((80, 24)).columns
    indent = " " * (len(label) + 2)   # align continuation with '[' after label
    items = [str(x) for x in perm]
    lines: List[str] = []
    current = f"  {label}: ["
    for k, item in enumerate(items):
        sep = "" if k == 0 else ", "
        piece = sep + item
        if k > 0 and len(current) + len(piece) + 1 > term_w:
            lines.append(current + ",")
            current = "  " + indent + item
        else:
            current += piece
    current += "]" + (f"  {suffix}" if suffix else "")
    lines.append(current)
    return lines


# ── permutation editor ────────────────────────────────────────────────────────

def _perm_editor_pt(
    n_apm: int,
    existing: Optional[List[int]] = None,
    atom_names: Optional[List] = None,
) -> Optional[List[int]]:
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

    ipp, n_pages = _perm_pages(n_apm, atom_names)
    display_page = [0]

    def _follow_cursor() -> None:
        display_page[0] = cursor[0] // ipp

    def _clear() -> None:
        input_buf.set_document(Document(""))

    @kb.add("left")
    def _(event):
        cursor[0] = max(0, cursor[0] - 1)
        _follow_cursor(); _clear(); msg[0] = ""

    @kb.add("right")
    def _(event):
        cursor[0] = min(n_apm - 1, cursor[0] + 1)
        _follow_cursor(); _clear(); msg[0] = ""

    @kb.add("up")
    def _(event):
        display_page[0] = max(0, display_page[0] - 1)

    @kb.add("down")
    def _(event):
        display_page[0] = min(n_pages - 1, display_page[0] + 1)

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
            msg[0] = f"! expected integer, got {text!r}"; _clear(); return
        if not (0 <= target < n_apm):
            msg[0] = f"! {target} out of range 0..{n_apm - 1}"; _clear(); return
        if target == cursor[0]:
            msg[0] = "! same position - nothing to swap"; _clear(); return
        mapping[cursor[0]], mapping[target] = mapping[target], mapping[cursor[0]]
        msg[0] = f"swapped pos {cursor[0]} <-> pos {target}"
        _clear()

    @kb.add("c-c")
    @kb.add("c-d")
    def _(event):
        event.app.exit()

    def get_display() -> str:
        dp = min(display_page[0], n_pages - 1)
        lines = ["", _render_display(n_apm, mapping, cursor[0], atom_names, dp, ipp, n_pages)]
        if msg[0]:
            lines.append(f"  {msg[0]}")
        elif mapping == list(range(n_apm)):
            lines.append("  [identity - swap pairs, or Enter to accept as-is]")
        elif n_pages > 1:
            lines.append("  [Enter confirm   <- -> move   Up/Down page   target+Enter swap   Ctrl+C cancel]")
        else:
            lines.append("  [Enter confirm   <- -> move   target+Enter swap   Ctrl+C cancel]")
        lines.append(f"\n  Swap pos {cursor[0]} with: ")
        return "\n".join(lines)

    layout = Layout(HSplit([
        Window(content=FormattedTextControl(get_display)),
        Window(content=BufferControl(buffer=input_buf), height=1),
    ]))
    Application(layout=layout, key_bindings=kb, full_screen=False, mouse_support=False).run()
    return result[0]


def _perm_editor_simple(
    n_apm: int,
    existing: Optional[List[int]] = None,
    atom_names: Optional[List] = None,
) -> Optional[List[int]]:
    mapping: List[int] = list(existing) if existing else list(range(n_apm))
    cursor = 0
    ipp, n_pages = _perm_pages(n_apm, atom_names)
    dp = [0]

    def _follow_cursor():
        dp[0] = cursor // ipp

    while True:
        print(_render_display(n_apm, mapping, cursor, atom_names, dp[0], ipp, n_pages))
        if n_pages > 1:
            print("  (n/p move   +/- page   target+Enter swap   blank=confirm   !reset   !cancel)")
        else:
            print("  (n/p move   target+Enter swap   blank=confirm   !reset   !cancel)")
        try:
            raw = input(f"  Swap pos {cursor} with: ").strip()
        except (EOFError, KeyboardInterrupt):
            raise _Abort()
        low = raw.lower()
        if low == "!cancel":
            raise _Abort()
        if low == "!reset":
            mapping = list(range(n_apm)); continue
        if low in ("n", "next"):
            cursor = min(n_apm - 1, cursor + 1); _follow_cursor(); continue
        if low in ("p", "prev"):
            cursor = max(0, cursor - 1); _follow_cursor(); continue
        if low in ("+", "pg+"):
            dp[0] = min(n_pages - 1, dp[0] + 1); continue
        if low in ("-", "pg-"):
            dp[0] = max(0, dp[0] - 1); continue
        if not raw:
            return mapping.copy()
        try:
            target = int(raw)
        except ValueError:
            print("  ! expected integer"); continue
        if not (0 <= target < n_apm):
            print(f"  ! out of range 0..{n_apm - 1}"); continue
        if target == cursor:
            print("  ! same position"); continue
        mapping[cursor], mapping[target] = mapping[target], mapping[cursor]
        print(f"  swapped pos {cursor} <-> pos {target}")


def _perm_editor(n_apm, existing=None, atom_names=None):
    if _PT:
        return _perm_editor_pt(n_apm, existing, atom_names)
    return _perm_editor_simple(n_apm, existing, atom_names)


# ── atom_names editor ─────────────────────────────────────────────────────────

def _atom_names_editor_pt(
    n_apm: int, existing: Optional[List] = None
) -> Optional[List[str]]:
    """Cursor-based atom_names editor.

    * Left/Right: move cursor
    * name+Enter: set at cursor, advance
    * blank+Enter: auto-assign at cursor, advance
    * Enter at last position / Down / Ctrl+D: confirm all and exit
    * Up: _Back   Ctrl+C: _Abort (skip step)

    On exit: clears Application output, prints 'atomnames >> [...]'.
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
    last_render_lines = [0]
    msg: List[str] = [""]
    input_buf = Buffer(name="name_input", multiline=False)
    kb = KeyBindings()

    def _clear() -> None:
        input_buf.set_document(Document(""))

    def _set_current(text: str) -> bool:
        """Validate and set name at cursor.  Sets msg[0] on failure."""
        msg[0] = ""
        if text:
            stripped = text.strip()
            if not stripped:
                msg[0] = "! empty name"
                return False
            if not stripped.isascii():
                msg[0] = "! non-ASCII characters not allowed"
                return False
            sanitized = stripped[:4]
            if len(stripped) > 4:
                msg[0] = f"! truncated to '{sanitized}'"
                # still accept, just warn — fall through
            used = [names[i] for i in range(n_apm) if i != cursor[0]]
            if sanitized in used:
                msg[0] = f"! '{sanitized}' already used at another position"
                return False
            names[cursor[0]] = sanitized
        else:
            if names[cursor[0]] is None:
                filled = _fill_atom_names(list(names), n_apm)
                names[cursor[0]] = filled[cursor[0]]
        return True

    def _do_confirm(event):
        result[0] = _fill_atom_names(list(names), n_apm)
        exit_reason[0] = "confirm"
        event.app.exit()

    @kb.add("left")
    def _(event):
        cursor[0] = max(0, cursor[0] - 1); msg[0] = ""; _clear()

    @kb.add("right")
    def _(event):
        cursor[0] = min(n_apm - 1, cursor[0] + 1); msg[0] = ""; _clear()

    @kb.add("enter")
    def _(event):
        if not _set_current(input_buf.text.strip()):
            _clear(); return
        _clear()
        if cursor[0] < n_apm - 1:
            cursor[0] += 1
        else:
            _do_confirm(event)

    @kb.add("down")
    @kb.add("c-d")
    def _(event):
        _do_confirm(event)

    @kb.add("up")
    def _(event):
        exit_reason[0] = "back"; event.app.exit()

    @kb.add("c-c")
    def _(event):
        exit_reason[0] = "abort"; event.app.exit()

    def get_display() -> str:
        body = _render_names_wrapped(n_apm, names, cursor[0])
        if msg[0]:
            body += f"\n  {msg[0]}"
        content = "\n" + body
        last_render_lines[0] = content.count('\n') + 1
        return content

    layout = Layout(HSplit([
        Window(content=FormattedTextControl(get_display)),
        Window(content=BufferControl(buffer=input_buf), height=1),
    ]))
    Application(layout=layout, key_bindings=kb, full_screen=False, mouse_support=False).run()

    n_up = last_render_lines[0] + 2
    sys.stdout.write("\n")
    _ansi_up_clear(n_up)

    reason = exit_reason[0]
    if reason == "back":
        raise _Back()
    if reason == "abort":
        raise _Abort()
    if result[0] is not None:
        sys.stdout.write(f"  atomnames >> {result[0]}\n")
        sys.stdout.flush()
    return result[0]


def _atom_names_editor_simple(
    n_apm: int, existing: Optional[List] = None
) -> Optional[List[str]]:
    names: List = list(existing) if existing else [None] * n_apm
    assert len(names) == n_apm
    cursor = 0

    _cursor_save()
    while True:
        print(_render_names_wrapped(n_apm, names, cursor))
        print("  (n/p move   name+Enter set   blank=auto   !done confirm   !back)")
        try:
            raw = input("  >> ").strip()
        except (EOFError, KeyboardInterrupt):
            _cursor_restore_clear(); raise _Abort()

        low = raw.lower()
        if low == "!back":
            _cursor_restore_clear(); raise _Back()
        if low in ("!done", "!d"):
            final = _fill_atom_names(list(names), n_apm)
            _cursor_restore_clear()
            print(f"  atomnames >> {final}")
            return final
        if low in ("n", "next"):
            cursor = min(n_apm - 1, cursor + 1)
            _cursor_restore_clear(); _cursor_save(); continue
        if low in ("p", "prev"):
            cursor = max(0, cursor - 1)
            _cursor_restore_clear(); _cursor_save(); continue

        if not raw:
            if names[cursor] is None:
                filled = _fill_atom_names(list(names), n_apm)
                names[cursor] = filled[cursor]
            if cursor < n_apm - 1:
                cursor += 1
                _cursor_restore_clear(); _cursor_save()
            else:
                final = _fill_atom_names(list(names), n_apm)
                _cursor_restore_clear()
                print(f"  atomnames >> {final}")
                return final
            continue

        sanitized = _sanitize_atom_name(raw)
        if not sanitized:
            _cursor_restore_clear(); _cursor_save(); continue
        used = [names[i] for i in range(n_apm) if i != cursor]
        if sanitized in used:
            _cursor_restore_clear(); _cursor_save(); continue
        names[cursor] = sanitized
        if cursor < n_apm - 1:
            cursor += 1
            _cursor_restore_clear(); _cursor_save()
        else:
            final = _fill_atom_names(list(names), n_apm)
            _cursor_restore_clear()
            print(f"  atomnames >> {final}")
            return final


def _atom_names_editor(n_apm, existing=None):
    if _PT:
        return _atom_names_editor_pt(n_apm, existing)
    return _atom_names_editor_simple(n_apm, existing)


# ── sym session ───────────────────────────────────────────────────────────────

def _sym_session(
    n_apm: int,
    existing: Optional[List[List[int]]] = None,
    atom_names: Optional[List] = None,
) -> Optional[List[List[int]]]:
    """Manage sym permutations interactively.

    Identity is always index 0.  User-added permutations are indexed 1, 2, ...
    Clears output on exit; prints 'sym >> [...]' on done.
    """
    sym: List[List[int]] = list(existing) if existing else []
    identity = list(range(n_apm))
    sym_prompt = _SymPrompt()

    def _show() -> None:
        print(_SEP)
        print(f"  Sym editor  (n_apm={n_apm}  n_sym={1 + len(sym)})")
        if atom_names:
            print(_trunc("  Names: " + "  ".join(str(n) for n in atom_names)))
        print(_trunc(f"    0: {identity}  (identity)"))
        for i, p in enumerate(sym):
            print(_trunc(f"    {i + 1}: {p}"))
        print(_SEP)
        print("  a(dd)  r(emove) N  l(ist)  d(one)/Enter  b(ack)/Up")

    _cursor_save()
    _show()

    def _cleanup(show_result: bool) -> None:
        _cursor_restore_clear()
        if show_result:
            sys.stdout.write(f"  sym >> {sym if sym else []}\n")
            sys.stdout.flush()

    def _do_remove(idx_str: str) -> bool:
        """Remove permutation by 1-based index string. Returns True on success."""
        if not idx_str.isdigit():
            print("  ! usage: remove N  (N >= 1)"); return False
        idx = int(idx_str)
        if idx < 1 or idx > len(sym):
            print(f"  ! index {idx} out of range 1..{len(sym)}"); return False
        sym.pop(idx - 1)
        return True

    while True:
        try:
            raw = sym_prompt("  sym >> ")
        except _Back:
            _cleanup(False); raise _Abort()
        except _Abort:
            _cleanup(False); raise _Abort()

        # Empty Enter = done
        if not raw:
            _cleanup(True)
            return sym or None

        parts = raw.split()
        cmd = parts[0].lower()
        args = parts[1:]

        # Number shorthand: "2" -> remove 2
        if cmd.isdigit():
            args = [cmd]; cmd = "r"

        if cmd in ("done", "d"):
            _cleanup(True)
            return sym or None

        elif cmd in ("!back", "back", "b"):
            _cleanup(False); raise _Abort()

        elif cmd in ("add", "a"):
            perm = _perm_editor(n_apm, atom_names=atom_names)
            if perm is None:
                pass
            elif perm == identity:
                print("  ! identity - not added")
            elif perm in sym:
                print("  ! duplicate - not added")
            else:
                sym.append(perm)
            _cursor_restore_clear(); _cursor_save(); _show()

        elif cmd in ("remove", "r"):
            if len(args) != 1:
                print("  ! usage: remove N  (N >= 1)"); continue
            if _do_remove(args[0]):
                _cursor_restore_clear(); _cursor_save(); _show()

        elif cmd in ("list", "l"):
            _cursor_restore_clear(); _cursor_save(); _show()

        else:
            print(f"  ! unknown: {cmd!r}.  Commands: a / r N / l / d / b")


# ── species wizard ────────────────────────────────────────────────────────────

def _wizard_species(wiz: _WizardPrompt, existing: Optional[Dict] = None) -> Dict:
    """Step-through wizard: n_apm -> n_mol -> name -> atom_names -> sym."""
    sp: Dict = copy.deepcopy(existing) if existing else {
        "n_apm": 1, "n_mol": 1, "name": None, "atom_names": None, "sym": None
    }
    if "atom_names" not in sp:
        sp["atom_names"] = None

    STEPS = ["n_apm", "n_mol", "name", "atom_names", "sym"]
    step = 0

    while step < len(STEPS):
        field = STEPS[step]

        if field == "atom_names":
            n = sp["n_apm"]
            existing_names = sp.get("atom_names") or [None] * n
            if len(existing_names) != n:
                existing_names = [None] * n
            try:
                new_names = _atom_names_editor(n, existing_names)
                if new_names is not None:
                    sp["atom_names"] = new_names
                step += 1
            except _Back:
                step -= 1
            except _Abort:
                step += 1  # Ctrl+C = skip, proceed to sym
            continue

        if field == "sym":
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

        if field == "n_apm":
            try:
                val = int(raw) if raw else (cur or 1)
                if val < 1:
                    print("    ! n_apm must be >= 1"); continue
                if existing and val != existing.get("n_apm"):
                    sp["atom_names"] = None; sp["sym"] = None
                sp["n_apm"] = val; step += 1
            except (ValueError, TypeError):
                print("    ! expected integer")

        elif field == "n_mol":
            try:
                val = int(raw) if raw else (cur or 1)
                if val < 1:
                    print("    ! n_mol must be >= 1"); continue
                sp["n_mol"] = val; step += 1
            except (ValueError, TypeError):
                print("    ! expected integer")

        elif field == "name":
            sp["name"] = None if raw == "-" else (raw if raw else cur)
            step += 1

    return sp


# ---------------------------------------------------------------------------
# TopologyBuilder
# ---------------------------------------------------------------------------

class TopologyBuilder:
    """Mutable builder for molecular_system.

    atom_names per species are stored here and stripped when converting to
    molecular_system.

    n_sym per species = 1 (identity) + number of user-added permutations.
    """

    def __init__(self) -> None:
        self._species: List[Dict] = []
        self._name: Optional[str] = None
        self._atom_order: Optional[List[int]] = None

    @classmethod
    def from_molecular_system(cls, mols) -> "TopologyBuilder":
        b = cls()
        b._name = mols.name
        b._atom_order = (
            mols.atom_order.tolist() if mols.atom_order is not None else None
        )
        for m in mols.mols:
            b._species.append({
                "n_apm": m.n_apm, "n_mol": m.n_mol,
                "sym": copy.deepcopy(m.sym), "name": m.name, "atom_names": None,
            })
        return b

    @classmethod
    def from_dict(cls, d: dict) -> "TopologyBuilder":
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
                "n_apm": m["n_apm"], "n_mol": m["n_mol"],
                "sym": m.get("sym"), "name": m.get("name"), "atom_names": m.get("atom_names"),
            })
        return b

    # ------------------------------------------------------------------
    # fluent API
    # ------------------------------------------------------------------

    def add(self, n_apm, n_mol, sym=None, name=None, atom_names=None) -> "TopologyBuilder":
        self._species.append({
            "n_apm": n_apm, "n_mol": n_mol, "sym": sym, "name": name, "atom_names": atom_names,
        })
        return self

    def remove(self, index: int) -> "TopologyBuilder":
        self._species.pop(index); return self

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
        self._species.sort(key=lambda sp: (sp["n_mol"], -sp["n_apm"]))
        return self

    def set_name(self, name) -> "TopologyBuilder":
        self._name = name; return self

    def set_atom_order(self, order) -> "TopologyBuilder":
        self._atom_order = list(order) if order is not None else None; return self

    def build(self):
        from .parse import to_molecular_system
        return to_molecular_system(self)

    def to_dict(self) -> dict:
        mol_list = []
        for sp in self._species:
            entry: dict = {"n_apm": sp["n_apm"], "n_mol": sp["n_mol"]}
            if sp.get("sym"):  entry["sym"]        = sp["sym"]
            if sp.get("name"): entry["name"]       = sp["name"]
            if sp.get("atom_names"): entry["atom_names"] = sp["atom_names"]
            mol_list.append(entry)
        inner: dict = {"molecules": mol_list}
        if self._name:          inner["name"]       = self._name
        if self._atom_order is not None: inner["atom_order"] = self._atom_order
        return {"molecular_system": inner}

    # ------------------------------------------------------------------
    # display
    # ------------------------------------------------------------------

    def _show(self) -> None:
        print(_SEP)
        print(f"  System name : {self._name or '(none)'}")
        if self._atom_order is not None:
            print(f"  atom_order  : {self._atom_order}")
        print(_SEP)
        if not self._species:
            print("  (no species defined)")
        else:
            print(f"  {'#':>3}  {'n_apm':>6}  {'n_mol':>6}  {'n_sym':>6}  name")
            for i, sp in enumerate(self._species):
                sym = sp.get("sym") or []
                n_sym = 1 + len(sym)
                print(f"  {i:>3}  {sp['n_apm']:>6}  {sp['n_mol']:>6}  {n_sym:>6}  {sp.get('name') or '(none)'}")
        print(_SEP)
        if self._species:
            print("  (use 'list N' for atom_names / sym detail)")

    def _show_species(self, idx: int) -> None:
        """Detailed view of species idx: atom_names and sym with line wrapping."""
        sp = self._species[idx]
        n_apm = sp["n_apm"]
        sym = sp.get("sym") or []
        an = sp.get("atom_names")

        lines: List[str] = []
        lines.append(
            f"  --- Species {idx}  n_apm={n_apm}  n_mol={sp['n_mol']}"
            f"  n_sym={1 + len(sym)}  name={sp.get('name')!r}"
        )

        # atom_names - column-wrapped
        if an:
            preview = _fill_atom_names(list(an), n_apm)
            term_w = shutil.get_terminal_size((80, 24)).columns
            col_w = max(4, max(len(n) for n in preview))
            items_per_row = max(1, (term_w - 4) // (col_w + 1))
            lines.append("  atom_names:")
            for chunk_start in range(0, n_apm, items_per_row):
                chunk = preview[chunk_start:chunk_start + items_per_row]
                lines.append("    " + " ".join(f"{n:<{col_w}}" for n in chunk))
        else:
            lines.append("  atom_names: (none)")

        # sym permutations - perm-wrapped
        lines.append(f"  sym  ({1 + len(sym)} total):")
        lines.extend(_wrap_perm_line(list(range(n_apm)), "0", "(identity)"))
        for j, p in enumerate(sym):
            lines.extend(_wrap_perm_line(p, str(j + 1)))

        if len(lines) >= 10:
            print(f"  ({len(lines)} lines total)")
            try:
                if input("  Show? [Y/n] ").strip().lower() == 'n':
                    return
            except (KeyboardInterrupt, EOFError):
                return

        for line in lines:
            print(line)

    # ------------------------------------------------------------------
    # interactive REPL
    # ------------------------------------------------------------------

    _HELP = """\
  Commands:
    list  (l)          - show species summary
    list  (l) N        - show species N detail (atom_names + sym)
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
        """Interactive REPL.  Returns built molecular_system or None on quit."""
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
                print(); continue
            except EOFError:
                print("\n  Quit."); return None

            if not raw:
                continue

            parts = raw.split()
            cmd = parts[0].lower()
            args = parts[1:]

            if cmd in ("help", "h", "?"):
                print(self._HELP)

            elif cmd in ("list", "l"):
                if args and args[0].isdigit():
                    idx = int(args[0])
                    if 0 <= idx < len(self._species):
                        self._show_species(idx)
                    else:
                        print(f"  ! index {idx} out of range")
                else:
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
                    print("  ! usage: edit N"); continue
                idx = int(args[0])
                if not (0 <= idx < len(self._species)):
                    print(f"  ! index {idx} out of range"); continue
                try:
                    sp = _wizard_species(wiz, existing=self._species[idx])
                    self._species[idx] = sp
                    self._show()
                except _Abort:
                    print("  edit cancelled.")

            elif cmd in ("remove", "r"):
                if len(args) != 1 or not args[0].isdigit():
                    print("  ! usage: remove N"); continue
                idx = int(args[0])
                if not (0 <= idx < len(self._species)):
                    print(f"  ! index {idx} out of range"); continue
                self.remove(idx); self._show()

            elif cmd in ("move", "mv"):
                if len(args) != 2 or not all(a.isdigit() for a in args):
                    print("  ! usage: move N M"); continue
                src, dst = int(args[0]), int(args[1])
                n = len(self._species)
                if not (0 <= src < n and 0 <= dst < n):
                    print(f"  ! indices out of range (0..{n - 1})"); continue
                sp = self._species.pop(src)
                self._species.insert(dst, sp)
                self._show()

            elif cmd in ("sort", "s"):
                self.sort(); self._show()

            elif cmd in ("name", "n"):
                self.set_name(" ".join(args) if args else None)
                print(f"  System name: {self._name!r}")

            elif cmd in ("order", "o"):
                raw_order = " ".join(args)
                if not raw_order:
                    self.set_atom_order(None); print("  atom_order cleared.")
                else:
                    try:
                        order = json.loads(raw_order)
                        if not isinstance(order, list) or not all(isinstance(x, int) for x in order):
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
                print("  Quit."); return None

            else:
                print(f"  ! unknown command: {cmd!r}. Type 'help'.")


# ---------------------------------------------------------------------------
# standalone entry point
# ---------------------------------------------------------------------------

def _main() -> None:
    import argparse

    ap = argparse.ArgumentParser(
        prog="python -m mobbrmsd.run.topology_builder",
        description="Interactive molecular_system editor.",
    )
    ap.add_argument("-t", "--topology", metavar="FILE", default=None)
    ap.add_argument("-c", "--coords",   metavar="FILE", default=None)
    ap.add_argument("-j", "--json",     metavar="FILE", default=None)
    ap.add_argument("-n", "--name",     metavar="NAME", default=None)
    ap.add_argument("-o", "--output",   metavar="FILE", default=None)
    ap.add_argument("--builder-json", action="store_true")
    args = ap.parse_args()

    builder = None

    if args.json is not None:
        try:
            with open(args.json) as f:
                builder = TopologyBuilder.from_dict(json.load(f))
        except Exception as e:
            print(f"ERROR: {e}", file=sys.stderr); sys.exit(1)

    elif args.topology is not None:
        try:
            import mdtraj
        except ImportError:
            print("ERROR: mdtraj required.", file=sys.stderr); sys.exit(1)
        try:
            from .parse import builder_from_topology, builder_from_mdtraj
            if args.coords is not None:
                traj = mdtraj.load(args.coords, top=args.topology)
                builder = builder_from_mdtraj(traj.topology, coords=traj.xyz * 10.0, name=args.name)
            else:
                builder = builder_from_topology(mdtraj.load_topology(args.topology), name=args.name)
        except Exception as e:
            print(f"ERROR: {e}", file=sys.stderr); sys.exit(1)

    if builder is None:
        builder = TopologyBuilder()
        if args.name:
            builder.set_name(args.name)

    ms = builder.run()
    if ms is None:
        sys.exit(1)

    output = json.dumps(builder.to_dict() if args.builder_json else ms.to_dict(), indent=2)
    if args.output is not None:
        with open(args.output, "w") as f:
            f.write(output); f.write("\n")
        print(f"  Written to {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    _main()
