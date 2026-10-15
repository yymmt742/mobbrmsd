"""mobbrmsd visualization utilities (requires matplotlib)."""
from __future__ import annotations

import numpy
import numpy.typing as npt
from pathlib import Path
from typing import Optional, Union

# Known keyword names accepted by plot_superposition (used by demo wrappers
# to filter their own kwargs before forwarding).
_PLOT_KW = frozenset({
    "cmap", "figsize", "lim", "azim_start", "frames", "interval",
    "line_kw", "marker_kw",
})


def _save_figure(output, fig, rot_func=None, rot_frames=None):
    """Save a matplotlib figure or animation to *output*.

    When *rot_func* and *rot_frames* are given the output is animated:
      - .gif / .webp  → single animated file (Pillow writer)
      - .png / .jpg   → numbered frame files  (foo.001.png, foo.002.png, …)
      - other         → static snapshot of the first frame

    When no animation function is given, ``fig.savefig(output)`` is used.
    """
    import sys

    p = Path(output)
    suffix = p.suffix.lower()

    if rot_func is None or rot_frames is None:
        fig.savefig(output)
        return

    rot_frames = list(rot_frames)

    if suffix in (".gif", ".webp"):
        import matplotlib.animation as _anim
        ani = _anim.FuncAnimation(fig, rot_func, frames=rot_frames, interval=50)
        ani.save(output, writer="pillow")
        print(f"  Animation saved to {output}", file=sys.stderr)
    elif suffix in (".png", ".jpg", ".jpeg"):
        width = max(3, len(str(len(rot_frames))))
        for k, frame_val in enumerate(rot_frames):
            rot_func(frame_val)
            frame_path = str(p.parent / f"{p.stem}.{k + 1:0{width}d}{p.suffix}")
            fig.savefig(frame_path)
        print(
            f"  {len(rot_frames)} frames saved as "
            f"{p.stem}.001{p.suffix} … {p.stem}.{len(rot_frames):0{width}d}{p.suffix}",
            file=sys.stderr,
        )
    else:
        rot_func(rot_frames[0])
        fig.savefig(output)
        print(f"  Figure saved to {output}", file=sys.stderr)


def plot_superposition(
    x: npt.NDArray,
    y_before: npt.NDArray,
    y_after: npt.NDArray,
    output: Optional[str] = None,
    *,
    cmap: Union[str, object] = "tab20",
    figsize: tuple = (8, 8 / 1.618),
    lim: float = 2.5,
    azim_start: float = 15.0,
    frames: Optional[npt.NDArray] = None,
    interval: int = 50,
    line_kw: Optional[dict] = None,
    marker_kw: Optional[dict] = None,
    **figure_kw,
):
    """3D animation of molecular superposition: before (left) and after (right).

    Parameters
    ----------
    x         : reference coordinates, shape ``(n_mol, n_apm, 3)``
    y_before  : target before superposition, same shape as *x*
    y_after   : target after superposition, same shape as *x*
    output    : save path (e.g. ``"result.webp"``); ``None`` shows interactively
    cmap      : matplotlib colormap name or object (default ``"tab20"``)
    figsize   : size forwarded to ``plt.figure`` (default ``(8, 8/φ)``)
    lim       : symmetric axis limit; ``lim=2.5`` → ``[-2.5, 2.5]``
    azim_start: starting azimuth angle for the rotating view (default ``15``)
    frames    : iterable of azimuth values for each animation frame
                (default ``numpy.arange(15, 374, 2)``)
    interval  : ms between frames when saving (default ``50``)
    line_kw   : extra kwargs merged into the displacement-line ``ax.plot`` calls;
                base defaults are ``ls=":"`` and ``lw=1.0``
    marker_kw : extra kwargs merged into the molecule-marker ``ax.plot`` calls;
                base defaults are ``marker="o"`` and ``ms=6``
    **figure_kw : forwarded verbatim to ``plt.figure`` (e.g. ``dpi=150``)

    Returns
    -------
    fig : ``matplotlib.figure.Figure``
        The caller is responsible for closing the figure (``plt.close(fig)``).
    """
    import matplotlib.pyplot as plt

    if frames is None:
        frames = numpy.arange(15, 374, 2)
    frames = list(frames)

    _line_kw   = {"ls": ":", "lw": 1.0}
    _marker_kw = {"marker": "o", "ms": 6}
    if line_kw:
        _line_kw.update(line_kw)
    if marker_kw:
        _marker_kw.update(marker_kw)

    cmap_fn = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap
    x        = numpy.asarray(x)
    y_before = numpy.asarray(y_before)
    y_after  = numpy.asarray(y_after)

    fig = plt.figure(figsize=figsize, **figure_kw)
    axes = [
        fig.add_subplot(1, 2, 1, projection="3d"),
        fig.add_subplot(1, 2, 2, projection="3d"),
    ]

    for ax, tgt in zip(axes, [y_before, y_after]):
        # displacement lines between reference and target atoms
        for i, (xi_mol, yi_mol) in enumerate(zip(x, tgt)):
            for xi, yi in zip(xi_mol, yi_mol):
                ax.plot(
                    [xi[0], yi[0]], [xi[1], yi[1]], [xi[2], yi[2]],
                    color=cmap_fn(2 * i),
                    **_line_kw,
                )
        # molecule markers (reference hollow, target filled)
        for coords, fillstyle, c_off in zip([x, tgt], ["none", "full"], [1, 0]):
            for i, mol in enumerate(coords):
                ax.plot(
                    mol[:, 0], mol[:, 1], mol[:, 2],
                    color=cmap_fn(2 * i + c_off),
                    fillstyle=fillstyle,
                    **_marker_kw,
                )
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_xlim([-lim, lim])
        ax.set_ylim([-lim, lim])
        ax.set_zlim([-lim, lim])
        ax.view_init(azim=azim_start)
        ax.set_box_aspect([1, 1, 1])

    plt.tight_layout()

    def _rot(frame):
        for ax in axes:
            ax.view_init(azim=frame)

    if output is None:
        import matplotlib.animation as _anim
        ani = _anim.FuncAnimation(fig, _rot, frames=frames, interval=interval)
        plt.show()
    else:
        _save_figure(output, fig, rot_func=_rot, rot_frames=frames)

    return fig
