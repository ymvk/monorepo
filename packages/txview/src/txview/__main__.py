"""CLI: python -m txview, or `txview` after install."""
from __future__ import annotations

import argparse
import shutil
import sys

from . import layout, parse_gtf, render


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="txview",
        description="Plot transcript models for a gene in the terminal.",
    )
    p.add_argument("gtf", help="Path to GTF file (plain or .gz).")
    p.add_argument("gene", help="Gene name or gene_id to plot.")
    p.add_argument(
        "--width",
        type=int,
        default=None,
        help="Display width in columns (default: terminal width minus label margin).",
    )
    p.add_argument(
        "--intron-width",
        type=int,
        default=5,
        help="Columns reserved for each truncated intron region (default: 5).",
    )
    args = p.parse_args(argv)

    width = args.width or max(60, shutil.get_terminal_size((120, 24)).columns - 30)

    transcripts = parse_gtf(args.gtf, args.gene)
    if not transcripts:
        print(
            f"No transcripts found for gene {args.gene!r} in {args.gtf}",
            file=sys.stderr,
        )
        return 1

    segments, used = layout(transcripts, total_width=width, intron_width=args.intron_width)
    print(render(transcripts, segments, used))
    return 0


if __name__ == "__main__":
    sys.exit(main())
