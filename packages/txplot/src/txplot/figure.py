"""Figure composition: `Figure`, `GenomicPanel`, `ProteinPanel`."""
from __future__ import annotations

from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np

from .space import Bridge, GenomicSpace, ProteinSpace
from .tracks import Track

__all__ = ["Figure", "GenomicPanel", "ProteinPanel"]


# ---------------------------------------------------------------------------
# Axis formatting
# ---------------------------------------------------------------------------


def _format_bp(pos: int) -> str:
    if abs(pos) >= 1_000_000:
        return f"{pos / 1_000_000:.2f} Mb"
    if abs(pos) >= 1_000:
        return f"{pos / 1_000:.1f} kb"
    return f"{pos} bp"


def _smart_ticks(start: int, end: int, n: int = 5) -> list[int]:
    if end <= start:
        return [start]
    span = end - start
    step_raw = span / (n - 1)
    magnitude = 10 ** int(np.floor(np.log10(step_raw))) if step_raw > 0 else 1
    step = magnitude
    for m in (1, 2, 5, 10):
        if m * magnitude >= step_raw:
            step = m * magnitude
            break
    first = (start // step + 1) * step
    ticks: list[int] = []
    pos = first
    while pos <= end:
        ticks.append(int(pos))
        pos += step
    if not ticks:
        ticks = [start, end]
    return ticks


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------


class _PanelBase:
    def __init__(self, space, tracks: Sequence[Track] | None = None):
        self.space = space
        self.tracks: list[Track] = list(tracks or [])

    def add(self, track: Track) -> "_PanelBase":
        self.tracks.append(track)
        return self

    @property
    def height_ratio(self) -> float:
        return sum(getattr(t, "height", 1.0) for t in self.tracks) or 1.0

    def _layout(self):
        heights = [getattr(t, "height", 1.0) for t in self.tracks]
        total = sum(heights) or 1.0
        slots: list[tuple[float, float]] = []
        cursor = total
        for h in heights:
            slots.append((cursor - h, cursor))
            cursor -= h
        return slots, total


class GenomicPanel(_PanelBase):
    """A horizontal panel sharing a `GenomicSpace`. Tracks stack vertically;
    a genomic-coordinate axis is drawn along the bottom."""

    def __init__(
        self,
        space: GenomicSpace,
        tracks: Sequence[Track] | None = None,
        show_axis: bool = True,
        title: str | None = None,
    ):
        super().__init__(space, tracks)
        self.show_axis = show_axis
        self.title = title

    def render(self, ax) -> None:
        slots, total = self._layout()
        for track, (yb, yt) in zip(self.tracks, slots):
            track.render(ax, self.space, yb, yt)

        ax.set_xlim(0, 1)
        ax.set_ylim(0, total)
        ax.set_yticks([])
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)

        if self.show_axis:
            ticks = _smart_ticks(self.space.start, self.space.end)
            ax.set_xticks([self.space.transform(t) for t in ticks])
            ax.set_xticklabels([_format_bp(t) for t in ticks], fontsize=7)
            ax.tick_params(axis="x", length=3, pad=2)
            ax.spines["bottom"].set_visible(True)
            ax.spines["bottom"].set_linewidth(0.6)
        else:
            ax.set_xticks([])
            ax.spines["bottom"].set_visible(False)

        if self.title:
            ax.set_title(self.title, fontsize=9, loc="left", pad=4)


class ProteinPanel(_PanelBase):
    """A horizontal panel sharing a `ProteinSpace`. An AA-position axis is
    drawn along the bottom."""

    def __init__(
        self,
        space: ProteinSpace,
        tracks: Sequence[Track] | None = None,
        show_axis: bool = True,
        title: str | None = None,
    ):
        super().__init__(space, tracks)
        self.show_axis = show_axis
        self.title = title

    def render(self, ax) -> None:
        slots, total = self._layout()
        for track, (yb, yt) in zip(self.tracks, slots):
            track.render(ax, self.space, yb, yt)

        ax.set_xlim(0, 1)
        ax.set_ylim(0, total)
        ax.set_yticks([])
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)

        if self.show_axis:
            ticks = _smart_ticks(1, self.space.length)
            ax.set_xticks([self.space.transform(t) for t in ticks])
            ax.set_xticklabels([f"{t} aa" for t in ticks], fontsize=7)
            ax.tick_params(axis="x", length=3, pad=2)
            ax.spines["bottom"].set_visible(True)
            ax.spines["bottom"].set_linewidth(0.6)
        else:
            ax.set_xticks([])
            ax.spines["bottom"].set_visible(False)

        if self.title:
            ax.set_title(self.title, fontsize=9, loc="left", pad=4)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------


class Figure:
    """Stacks panels (and bridges) vertically into a matplotlib Figure.

    Example
    -------
    >>> fig = Figure(
    ...     GenomicPanel(gs, [Transcripts(gtf, "BRCA1")]),
    ...     Bridge(gs, ps),
    ...     ProteinPanel(ps, [Domains(df)]),
    ... )
    >>> fig.save("brca1.svg")
    """

    def __init__(
        self,
        *components,
        figsize: tuple[float, float] = (10.0, 5.0),
        bridge_height_ratio: float = 0.8,
    ):
        self.components = list(components)
        self.figsize = figsize
        self.bridge_height_ratio = bridge_height_ratio

    def _ratios(self) -> list[float]:
        ratios: list[float] = []
        for c in self.components:
            if isinstance(c, Bridge):
                ratios.append(self.bridge_height_ratio)
            else:
                ratios.append(getattr(c, "height_ratio", 1.0))
        return ratios

    def render(self):
        n = len(self.components)
        if n == 0:
            raise ValueError("Figure has no components")
        fig, axes = plt.subplots(
            n, 1,
            figsize=self.figsize,
            gridspec_kw={"height_ratios": self._ratios()},
        )
        if n == 1:
            axes = [axes]
        for comp, ax in zip(self.components, axes):
            comp.render(ax)
        fig.subplots_adjust(hspace=0.15, left=0.1, right=0.95, top=0.92, bottom=0.08)
        return fig

    def save(self, path: str, dpi: int = 300, **kwargs) -> None:
        fig = self.render()
        fig.savefig(path, dpi=dpi, bbox_inches="tight", **kwargs)
        plt.close(fig)
