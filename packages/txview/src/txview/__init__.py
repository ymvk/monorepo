"""Plot transcript models for a gene from a GTF annotation, in the terminal."""
from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

__all__ = ["Transcript", "Segment", "parse_gtf", "layout", "render"]

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')

EXON_CHAR = "█"
INTRON_CHAR = "─"
TRUNC_CHAR = "╌"
TSS_PLUS_CHAR = "▶"
TSS_MINUS_CHAR = "◀"
TTS_CHAR = "│"


@dataclass
class Transcript:
    id: str
    seqname: str
    strand: str
    gene_name: str
    exons: list[tuple[int, int]] = field(default_factory=list)

    @property
    def start(self) -> int:
        return min(s for s, _ in self.exons)

    @property
    def end(self) -> int:
        return max(e for _, e in self.exons)


@dataclass
class Segment:
    """A chunk of the display: either an exonic genomic region (scaled) or a
    truncated intron region (fixed width)."""

    gstart: int
    gend: int
    dstart: int
    dend: int
    exonic: bool


def _open(path: str | Path):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt")
    return open(p, "rt")


def _parse_attrs(s: str) -> dict[str, str]:
    return dict(_ATTR_RE.findall(s))


def parse_gtf(path: str | Path, gene: str) -> dict[str, Transcript]:
    """Return all transcripts whose gene_id or gene_name matches `gene`."""
    transcripts: dict[str, Transcript] = {}
    with _open(path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            attrs = _parse_attrs(parts[8])
            if attrs.get("gene_id") != gene and attrs.get("gene_name") != gene:
                continue
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            t = transcripts.setdefault(
                tid,
                Transcript(
                    id=tid,
                    seqname=parts[0],
                    strand=parts[6],
                    gene_name=attrs.get("gene_name", attrs.get("gene_id", "")),
                ),
            )
            t.exons.append((int(parts[3]), int(parts[4])))
    for t in transcripts.values():
        t.exons.sort()
    return transcripts


def _merge(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    ivs = sorted(intervals)
    if not ivs:
        return []
    out = [ivs[0]]
    for s, e in ivs[1:]:
        ps, pe = out[-1]
        if s <= pe + 1:
            out[-1] = (ps, max(pe, e))
        else:
            out.append((s, e))
    return out


def layout(
    transcripts: dict[str, Transcript],
    total_width: int = 120,
    intron_width: int = 5,
    min_exon_width: int = 1,
) -> tuple[list[Segment], int]:
    """Map genomic coordinates to display columns, compressing introns to a
    fixed `intron_width` while keeping exonic regions roughly proportional."""
    all_exons: list[tuple[int, int]] = []
    for t in transcripts.values():
        all_exons.extend(t.exons)
    merged = _merge(all_exons)
    if not merged:
        return [], 0

    n_introns = max(0, len(merged) - 1)
    exonic_bp = sum(e - s + 1 for s, e in merged)
    intron_cols = n_introns * intron_width
    exon_cols_budget = max(len(merged) * min_exon_width, total_width - intron_cols)
    scale = exon_cols_budget / exonic_bp

    segments: list[Segment] = []
    col = 0
    prev_end: int | None = None
    for s, e in merged:
        if prev_end is not None:
            segments.append(Segment(prev_end + 1, s - 1, col, col + intron_width, False))
            col += intron_width
        width = max(min_exon_width, round((e - s + 1) * scale))
        segments.append(Segment(s, e, col, col + width, True))
        col += width
        prev_end = e
    return segments, col


def _project(pos: int, seg: Segment) -> int:
    """Map a genomic position to a display column within `seg`."""
    gspan = max(1, seg.gend - seg.gstart + 1)
    dspan = seg.dend - seg.dstart
    frac = (pos - seg.gstart) / gspan
    return seg.dstart + min(dspan - 1, max(0, int(frac * dspan)))


def _find_exonic_seg(pos: int, segments: Sequence[Segment]) -> Segment | None:
    for seg in segments:
        if seg.exonic and seg.gstart <= pos <= seg.gend:
            return seg
    return None


def render(
    transcripts: dict[str, Transcript],
    segments: Sequence[Segment],
    total_width: int,
) -> str:
    if not transcripts:
        return "(no transcripts)"

    tids = sorted(transcripts.keys(), key=lambda x: (transcripts[x].start, x))
    name_width = max(len(t) for t in tids)

    first = next(iter(transcripts.values()))
    gstart = min(t.start for t in transcripts.values())
    gend = max(t.end for t in transcripts.values())
    header = f"{first.gene_name} [{first.seqname}:{gstart}-{gend}] ({first.strand})"

    lines = [header]
    for tid in tids:
        t = transcripts[tid]
        row = [" "] * total_width

        seg_l = _find_exonic_seg(t.start, segments)
        seg_r = _find_exonic_seg(t.end, segments)
        tleft = _project(t.start, seg_l) if seg_l else 0
        tright = (_project(t.end, seg_r) + 1) if seg_r else total_width

        for c in range(tleft, tright):
            row[c] = INTRON_CHAR

        for seg in segments:
            if seg.exonic or t.start > seg.gend or t.end < seg.gstart:
                continue
            for c in range(seg.dstart, seg.dend):
                row[c] = TRUNC_CHAR

        for es, ee in t.exons:
            for seg in segments:
                if not seg.exonic or es > seg.gend or ee < seg.gstart:
                    continue
                ov_s = max(es, seg.gstart)
                ov_e = min(ee, seg.gend)
                ds = _project(ov_s, seg)
                de = _project(ov_e, seg) + 1
                for c in range(ds, de):
                    row[c] = EXON_CHAR

        # TSS/TTS markers: TSS follows strand direction (5' end), TTS is 3' end.
        if t.strand == "-":
            tss_pos, tts_pos = t.end, t.start
            tss_mark = TSS_MINUS_CHAR
        else:
            tss_pos, tts_pos = t.start, t.end
            tss_mark = TSS_PLUS_CHAR
        tss_seg = _find_exonic_seg(tss_pos, segments)
        tts_seg = _find_exonic_seg(tts_pos, segments)
        if tss_seg:
            row[_project(tss_pos, tss_seg)] = tss_mark
        if tts_seg:
            row[_project(tts_pos, tts_seg)] = TTS_CHAR

        lines.append(f"{tid:<{name_width}}  {''.join(row).rstrip()}")

    return "\n".join(lines)
