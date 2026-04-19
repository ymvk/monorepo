"""Coordinate systems: `GenomicSpace`, `ProteinSpace`, and `Bridge` connecting them."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ._gtf import TranscriptAnnotation, parse_gene

__all__ = ["GenomicSpace", "ProteinSpace", "Bridge", "Segment"]


@dataclass
class Segment:
    """A piece of the display axis: an exonic genomic region (scaled) or a
    truncated intron region (fixed fractional width)."""
    gstart: int
    gend: int
    dstart: float
    dend: float
    exonic: bool


def _merge(intervals) -> list[tuple[int, int]]:
    ivs = sorted(intervals)
    if not ivs:
        return []
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


class GenomicSpace:
    """Maps genomic positions (bp) onto a [0, 1] figure coordinate.

    Parameters
    ----------
    seqname, start, end
        The locus being displayed (1-based inclusive).
    exon_regions
        Optional list of (start, end) intervals representing the union of
        exons for the gene. Required when ``compress_introns=True``.
    compress_introns
        If ``True``, introns between exonic regions are collapsed to a fixed
        fraction of total width (`intron_width_frac` each), and exonic
        regions share the remainder proportionally.
    intron_width_frac
        Fraction of [0, 1] reserved per truncated intron (default 0.04).
    """

    def __init__(
        self,
        seqname: str,
        start: int,
        end: int,
        exon_regions: Sequence[tuple[int, int]] | None = None,
        compress_introns: bool = True,
        intron_width_frac: float = 0.04,
    ):
        self.seqname = seqname
        self.start = start
        self.end = end
        self.compress_introns = compress_introns and bool(exon_regions)
        self.transcripts: dict[str, TranscriptAnnotation] = {}

        if not self.compress_introns:
            self._segments = [Segment(start, end, 0.0, 1.0, True)]
            return

        merged = _merge(exon_regions or [])
        n_introns = max(0, len(merged) - 1)
        intron_total = n_introns * intron_width_frac
        exonic_total = max(0.1, 1.0 - intron_total)
        exonic_bp = sum(e - s + 1 for s, e in merged)

        segments: list[Segment] = []
        cursor = 0.0
        prev_end: int | None = None
        for s, e in merged:
            if prev_end is not None:
                segments.append(Segment(prev_end + 1, s - 1, cursor, cursor + intron_width_frac, False))
                cursor += intron_width_frac
            w = (e - s + 1) / exonic_bp * exonic_total
            segments.append(Segment(s, e, cursor, cursor + w, True))
            cursor += w
            prev_end = e
        # Renormalize to exactly [0, 1].
        if cursor > 0 and cursor != 1.0:
            for seg in segments:
                seg.dstart /= cursor
                seg.dend /= cursor
        self._segments = segments

    @classmethod
    def from_gtf(
        cls,
        gtf: str,
        gene: str,
        compress_introns: bool = True,
        intron_width_frac: float = 0.04,
        pad: int = 0,
    ) -> "GenomicSpace":
        transcripts = parse_gene(gtf, gene)
        if not transcripts:
            raise ValueError(f"No transcripts found for gene {gene!r} in {gtf}")
        first = next(iter(transcripts.values()))
        starts = [t.start for t in transcripts.values()]
        ends = [t.end for t in transcripts.values()]
        all_exons: list[tuple[int, int]] = []
        for t in transcripts.values():
            all_exons.extend(t.exons)
        instance = cls(
            seqname=first.seqname,
            start=min(starts) - pad,
            end=max(ends) + pad,
            exon_regions=all_exons,
            compress_introns=compress_introns,
            intron_width_frac=intron_width_frac,
        )
        instance.transcripts = transcripts
        return instance

    @property
    def segments(self) -> list[Segment]:
        return list(self._segments)

    def transform(self, pos: int) -> float:
        """Map a genomic position to display coordinate in [0, 1]."""
        if pos <= self._segments[0].gstart:
            return self._segments[0].dstart
        if pos >= self._segments[-1].gend:
            return self._segments[-1].dend
        for seg in self._segments:
            if seg.gstart <= pos <= seg.gend:
                gspan = max(1, seg.gend - seg.gstart + 1)
                dspan = seg.dend - seg.dstart
                return seg.dstart + (pos - seg.gstart) / gspan * dspan
        # In an uncovered region (e.g. between merged exons when compressing):
        # find the nearest segment.
        for prev, nxt in zip(self._segments, self._segments[1:]):
            if prev.gend < pos < nxt.gstart:
                return (prev.dend + nxt.dstart) / 2
        return 0.0


class ProteinSpace:
    """Maps 1-based amino acid positions onto [0, 1]."""

    def __init__(
        self,
        length: int,
        cds_intervals: Sequence[tuple[int, int]] | None = None,
        strand: str = "+",
    ):
        if length < 1:
            raise ValueError("ProteinSpace length must be >= 1")
        self.length = length
        self.cds_intervals = list(cds_intervals or [])
        self.strand = strand

    @classmethod
    def from_transcript(cls, gs: GenomicSpace, transcript_id: str) -> "ProteinSpace":
        t = gs.transcripts.get(transcript_id)
        if t is None:
            raise ValueError(f"Transcript {transcript_id!r} not found on GenomicSpace")
        if not t.cds:
            raise ValueError(f"Transcript {transcript_id!r} has no CDS features")
        total = sum(e - s + 1 for s, e in t.cds)
        return cls(length=total // 3, cds_intervals=t.cds, strand=t.strand)

    def transform(self, aa: int) -> float:
        if self.length == 1:
            return 0.5
        clamped = max(1, min(self.length, aa))
        return (clamped - 1) / (self.length - 1)

    def cds_to_aa(self) -> list[tuple[int, int, int, int]]:
        """Yield (genomic_start, genomic_end, aa_start, aa_end) per CDS chunk,
        respecting strand (for - strand, largest genomic coord is aa=1)."""
        ordered = sorted(self.cds_intervals, reverse=(self.strand == "-"))
        out: list[tuple[int, int, int, int]] = []
        cumulative = 0
        for s, e in ordered:
            length_bp = e - s + 1
            aa_start = cumulative // 3 + 1
            aa_end = (cumulative + length_bp - 1) // 3 + 1
            out.append((s, e, aa_start, aa_end))
            cumulative += length_bp
        return out


class Bridge:
    """Connects exon boundaries in a `GenomicSpace` to corresponding amino
    acid ranges in a `ProteinSpace` via filled trapezoids. Renders into its
    own Axes between the two panels."""

    def __init__(self, genomic: GenomicSpace, protein: ProteinSpace, alpha: float = 0.25, color: str = "#56B4E9"):
        self.genomic = genomic
        self.protein = protein
        self.alpha = alpha
        self.color = color

    def render(self, ax) -> None:
        from matplotlib.patches import Polygon

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")

        for gs, ge, aas, aae in self.protein.cds_to_aa():
            gx0 = self.genomic.transform(gs)
            gx1 = self.genomic.transform(ge)
            px0 = self.protein.transform(aas)
            px1 = self.protein.transform(aae)
            # Trapezoid: top edge = genomic range, bottom edge = protein range
            verts = [(gx0, 1.0), (gx1, 1.0), (px1, 0.0), (px0, 0.0)]
            ax.add_patch(
                Polygon(verts, closed=True, facecolor=self.color,
                        edgecolor="none", alpha=self.alpha)
            )
