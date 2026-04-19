"""Track classes rendered into a panel Axes.

All tracks accept DataFrame inputs with documented schemas. See the package
README for the expected columns per track.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from . import palette
from ._gtf import parse_gene
from .space import GenomicSpace, ProteinSpace

__all__ = ["Track", "Transcripts", "Lollipops", "Sashimi", "Coverage", "Domains"]


class Track:
    """Tracks declare a relative ``height`` and render into a horizontal slot
    [y_bottom, y_top] of their panel's Axes."""

    height: float = 1.0

    def render(self, ax, space, y_bottom: float, y_top: float) -> None:
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Transcripts
# ---------------------------------------------------------------------------


def _split_exon_by_cds(
    es: int, ee: int, cds_intervals: Sequence[tuple[int, int]]
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Split an exon [es, ee] by overlapping CDS intervals.

    Returns ``(utr_parts, cds_parts)``, two lists of disjoint genomic intervals.
    """
    cds_within = [
        (max(cs, es), min(ce, ee))
        for cs, ce in cds_intervals
        if cs <= ee and ce >= es
    ]
    if not cds_within:
        return [(es, ee)], []
    cds_within.sort()
    utr: list[tuple[int, int]] = []
    p = es
    for cs, ce in cds_within:
        if p < cs:
            utr.append((p, cs - 1))
        p = ce + 1
    if p <= ee:
        utr.append((p, ee))
    return utr, cds_within


class Transcripts(Track):
    """Transcript models with CDS/UTR distinction and directional intron chevrons.

    Parameters
    ----------
    gtf
        Path to a GTF file (or the same one used for the space).
    gene
        Gene name or gene_id to plot.
    cds_utr
        If ``True`` (default), draw CDS as tall rectangles and UTRs as thin
        rectangles. If ``False``, all exons are drawn at the full height.
    chevrons
        If ``True``, stamp ``>`` or ``<`` glyphs along the intron line in the
        direction of transcription.
    label
        If ``True`` (default), label each row with its transcript_id to the left.
    color
        Base color for transcript bodies (default Okabe-Ito blue).
    """

    def __init__(
        self,
        gtf: str,
        gene: str,
        cds_utr: bool = True,
        chevrons: bool = True,
        label: bool = True,
        color: str = "#0072B2",
    ):
        self.gtf = gtf
        self.gene = gene
        self.cds_utr = cds_utr
        self.chevrons = chevrons
        self.label = label
        self.color = color
        self._transcripts: dict | None = None

    def _annots(self):
        if self._transcripts is None:
            self._transcripts = parse_gene(self.gtf, self.gene)
        return self._transcripts

    @property
    def height(self) -> float:
        return max(1, len(self._annots())) * 0.8

    def render(self, ax, space: GenomicSpace, y_bottom: float, y_top: float) -> None:
        transcripts = self._annots()
        if not transcripts:
            return

        tids = sorted(transcripts.keys(), key=lambda x: (transcripts[x].start, x))
        n = len(tids)
        slot_h = (y_top - y_bottom) / n
        cds_h = slot_h * 0.55
        utr_h = slot_h * 0.28

        for i, tid in enumerate(tids):
            t = transcripts[tid]
            center = y_top - (i + 0.5) * slot_h

            left = space.transform(t.start)
            right = space.transform(t.end)
            ax.plot([left, right], [center, center], color="#444", lw=0.8, zorder=1)

            if self.chevrons:
                marker = ">" if t.strand == "+" else "<"
                for (s1, e1), (s2, e2) in zip(
                    sorted(t.exons), sorted(t.exons)[1:]
                ):
                    ds = space.transform(e1)
                    de = space.transform(s2)
                    span = de - ds
                    if span <= 0:
                        continue
                    k = max(1, int(span / 0.025))
                    xs = np.linspace(ds, de, k + 2)[1:-1]
                    ax.plot(
                        xs,
                        [center] * len(xs),
                        linestyle="",
                        marker=marker,
                        color="#444",
                        ms=3.5,
                        zorder=2,
                    )

            cds_set = list(t.cds) if self.cds_utr else []
            for es, ee in t.exons:
                utr_parts, cds_parts = _split_exon_by_cds(es, ee, cds_set)
                if not self.cds_utr:
                    utr_parts, cds_parts = [], [(es, ee)]
                for ps, pe in utr_parts:
                    ds = space.transform(ps)
                    de = space.transform(pe)
                    ax.add_patch(
                        Rectangle(
                            (ds, center - utr_h / 2),
                            max(1e-4, de - ds),
                            utr_h,
                            facecolor=self.color,
                            edgecolor="#222",
                            lw=0.4,
                            zorder=3,
                        )
                    )
                for ps, pe in cds_parts:
                    ds = space.transform(ps)
                    de = space.transform(pe)
                    ax.add_patch(
                        Rectangle(
                            (ds, center - cds_h / 2),
                            max(1e-4, de - ds),
                            cds_h,
                            facecolor=self.color,
                            edgecolor="#222",
                            lw=0.5,
                            zorder=3,
                        )
                    )

            if self.label:
                ax.text(
                    -0.005,
                    center,
                    tid,
                    ha="right",
                    va="center",
                    fontsize=7,
                    color="#222",
                )


# ---------------------------------------------------------------------------
# Lollipops
# ---------------------------------------------------------------------------


class Lollipops(Track):
    """Mutation lollipops.

    Expected DataFrame schema:

    ========  ===================================================
    column    description
    ========  ===================================================
    position  Genomic (1-based) or protein AA position.
    height    Optional numeric (VAF or recurrence count). If
              missing, defaults to 1 per variant.
    category  Optional categorical label for coloring (e.g.
              ``"missense"``, ``"nonsense"``, ...).
    label     Optional text drawn next to the head.
    ========  ===================================================

    Parameters
    ----------
    data
        DataFrame with the columns above.
    coord
        ``"genomic"`` (default) uses the `GenomicSpace`; ``"protein"`` uses
        the `ProteinSpace` of the panel.
    color_by
        Column name to color heads by, or ``None`` for a single color.
    colors
        Mapping from category -> hex color, or ``None`` to use
        :data:`palette.CONSEQUENCE_COLORS` + ``Okabe-Ito`` fallback.
    """

    height = 1.5

    def __init__(
        self,
        data: pd.DataFrame,
        coord: str = "genomic",
        color_by: str | None = "category",
        colors: dict[str, str] | None = None,
    ):
        if coord not in ("genomic", "protein"):
            raise ValueError("coord must be 'genomic' or 'protein'")
        self.data = data
        self.coord = coord
        self.color_by = color_by
        self.colors = colors

    def render(self, ax, space, y_bottom: float, y_top: float) -> None:
        if self.data.empty:
            return

        heights = self.data.get("height", pd.Series([1.0] * len(self.data)))
        max_h = float(max(1.0, heights.max()))
        baseline = y_bottom
        top_y = y_top

        colors = {**palette.CONSEQUENCE_COLORS, **(self.colors or {})}
        fallback = palette.OKABE_ITO

        for i, row in enumerate(self.data.itertuples(index=False)):
            pos = int(getattr(row, "position"))
            h = float(getattr(row, "height", 1.0) or 1.0)
            y = baseline + (h / max_h) * (top_y - baseline) * 0.9
            x = space.transform(pos)

            if self.color_by and hasattr(row, self.color_by):
                cat = getattr(row, self.color_by)
                c = colors.get(str(cat), fallback[i % len(fallback)])
            else:
                c = palette.OKABE_ITO[0]

            ax.plot([x, x], [baseline, y], color="#444", lw=0.6, zorder=2)
            ax.plot(
                [x],
                [y],
                marker="o",
                markersize=5,
                markerfacecolor=c,
                markeredgecolor="#222",
                markeredgewidth=0.4,
                zorder=3,
                linestyle="",
            )
            lbl = getattr(row, "label", None)
            if isinstance(lbl, str) and lbl:
                ax.text(x, y + 0.05 * (top_y - baseline), lbl,
                        ha="center", va="bottom", fontsize=6, color="#222")


# ---------------------------------------------------------------------------
# Sashimi
# ---------------------------------------------------------------------------


class Sashimi(Track):
    """Splice junction arcs with thickness proportional to read count.

    Expected DataFrame schema:

    ========  =============================================
    column    description
    ========  =============================================
    start     Genomic position of the 5' junction end.
    end       Genomic position of the 3' junction end.
    count     Read count supporting the junction.
    ========  =============================================

    Parameters
    ----------
    data
        DataFrame with the columns above.
    scale
        ``"linear"`` or ``"log"`` for mapping counts to arc widths.
    max_lw
        Maximum arc linewidth in points.
    """

    height = 2.0

    def __init__(
        self,
        data: pd.DataFrame,
        scale: str = "linear",
        max_lw: float = 3.0,
        color: str = "#0072B2",
    ):
        if scale not in ("linear", "log"):
            raise ValueError("scale must be 'linear' or 'log'")
        self.data = data
        self.scale = scale
        self.max_lw = max_lw
        self.color = color

    def render(self, ax, space: GenomicSpace, y_bottom: float, y_top: float) -> None:
        if self.data.empty:
            return
        counts = self.data["count"].to_numpy(dtype=float)
        if self.scale == "log":
            transformed = np.log1p(counts)
        else:
            transformed = counts
        max_t = float(max(1.0, transformed.max()))

        baseline = y_bottom + 0.02 * (y_top - y_bottom)
        peak = y_top - 0.02 * (y_top - y_bottom)
        for (_, row), t in zip(self.data.iterrows(), transformed):
            x0 = space.transform(int(row["start"]))
            x1 = space.transform(int(row["end"]))
            if x1 <= x0:
                continue
            lw = 0.5 + (t / max_t) * (self.max_lw - 0.5)
            # Quadratic Bezier: two ends at baseline, peak at midpoint.
            midx = (x0 + x1) / 2
            midy = peak
            ts = np.linspace(0, 1, 40)
            xs = (1 - ts) ** 2 * x0 + 2 * (1 - ts) * ts * midx + ts ** 2 * x1
            ys = (1 - ts) ** 2 * baseline + 2 * (1 - ts) * ts * midy + ts ** 2 * baseline
            ax.plot(xs, ys, color=self.color, lw=lw, zorder=2)
            ax.text(midx, midy, str(int(row["count"])),
                    ha="center", va="bottom", fontsize=6, color="#222", zorder=3)


# ---------------------------------------------------------------------------
# Coverage
# ---------------------------------------------------------------------------


class Coverage(Track):
    """Per-base read depth as a filled step histogram.

    Expected DataFrame schema:

    ========  =============================================
    column    description
    ========  =============================================
    position  Genomic position (1-based).
    depth     Read depth (nonnegative).
    ========  =============================================
    """

    height = 1.2

    def __init__(self, data: pd.DataFrame, color: str = "#56B4E9"):
        self.data = data.sort_values("position").reset_index(drop=True)
        self.color = color

    def render(self, ax, space: GenomicSpace, y_bottom: float, y_top: float) -> None:
        if self.data.empty:
            return
        xs = np.array([space.transform(int(p)) for p in self.data["position"]])
        ys = self.data["depth"].to_numpy(dtype=float)
        if ys.max() <= 0:
            return
        y_norm = y_bottom + (ys / ys.max()) * (y_top - y_bottom) * 0.95
        ax.fill_between(xs, y_bottom, y_norm, color=self.color, alpha=0.7,
                        linewidth=0, zorder=2)
        # Max-depth tick on the right.
        ax.text(1.005, y_top, f"{int(ys.max())}",
                ha="left", va="top", fontsize=6, color="#444",
                transform=ax.get_yaxis_transform() if False else ax.transData)


# ---------------------------------------------------------------------------
# Domains
# ---------------------------------------------------------------------------


class Domains(Track):
    """Protein domain rectangles (Pfam / InterPro style).

    Expected DataFrame schema:

    ========  =============================================
    column    description
    ========  =============================================
    name      Domain name (drawn as label inside the box).
    start     AA start (1-based inclusive).
    end       AA end (inclusive).
    ========  =============================================
    """

    height = 1.0

    def __init__(self, data: pd.DataFrame, colors: dict[str, str] | None = None):
        self.data = data
        self.colors = colors

    def render(self, ax, space: ProteinSpace, y_bottom: float, y_top: float) -> None:
        if self.data.empty:
            return
        h = (y_top - y_bottom) * 0.6
        cy = (y_top + y_bottom) / 2
        # Backbone: thin line representing the whole protein.
        ax.plot([0, 1], [cy, cy], color="#888", lw=0.8, zorder=1)

        names = list(self.data["name"])
        if self.colors:
            color_map = self.colors
        else:
            cs = palette.categorical(len(set(names)))
            color_map = {n: c for n, c in zip(sorted(set(names)), cs)}

        for row in self.data.itertuples(index=False):
            x0 = space.transform(int(row.start))
            x1 = space.transform(int(row.end))
            c = color_map.get(row.name, palette.OKABE_ITO[0])
            ax.add_patch(
                Rectangle(
                    (x0, cy - h / 2),
                    max(1e-4, x1 - x0),
                    h,
                    facecolor=c,
                    edgecolor="#222",
                    lw=0.5,
                    zorder=2,
                )
            )
            ax.text((x0 + x1) / 2, cy, row.name,
                    ha="center", va="center", fontsize=7, color="#111", zorder=3)
