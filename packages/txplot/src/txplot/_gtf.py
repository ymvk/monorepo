"""Minimal GTF parser. Extracts exon and CDS features for a named gene."""
from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path

__all__ = ["TranscriptAnnotation", "parse_gene"]

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


@dataclass
class TranscriptAnnotation:
    id: str
    seqname: str
    strand: str
    gene_name: str
    exons: list[tuple[int, int]] = field(default_factory=list)
    cds: list[tuple[int, int]] = field(default_factory=list)

    @property
    def start(self) -> int:
        return min(s for s, _ in self.exons)

    @property
    def end(self) -> int:
        return max(e for _, e in self.exons)


def _open(path: str | Path):
    p = Path(path)
    return gzip.open(p, "rt") if p.suffix == ".gz" else open(p, "rt")


def _parse_attrs(s: str) -> dict[str, str]:
    return dict(_ATTR_RE.findall(s))


def parse_gene(path: str | Path, gene: str) -> dict[str, TranscriptAnnotation]:
    """Return all transcripts whose gene_id or gene_name matches `gene`.
    Parses exon and CDS feature rows; ignores everything else."""
    transcripts: dict[str, TranscriptAnnotation] = {}
    with _open(path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature not in ("exon", "CDS"):
                continue
            attrs = _parse_attrs(parts[8])
            if attrs.get("gene_id") != gene and attrs.get("gene_name") != gene:
                continue
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            t = transcripts.setdefault(
                tid,
                TranscriptAnnotation(
                    id=tid,
                    seqname=parts[0],
                    strand=parts[6],
                    gene_name=attrs.get("gene_name", attrs.get("gene_id", "")),
                ),
            )
            interval = (int(parts[3]), int(parts[4]))
            if feature == "exon":
                t.exons.append(interval)
            else:
                t.cds.append(interval)
    for t in transcripts.values():
        t.exons.sort()
        t.cds.sort()
    return transcripts
