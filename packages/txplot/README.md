# txplot

Transcript model figures built on matplotlib. Distinct from
existing libraries in that **coordinate systems are first-class objects**:
genomic and protein spaces compose, and the joint genomic-protein "exon bridge"
view is a renderable object rather than a side-effect.

## Install

```bash
pip install "txplot @ git+https://github.com/ymvk/monorepo@main#subdirectory=packages/txplot"
```

## Quick start

```python
import pandas as pd
import txplot as tx

gs = tx.GenomicSpace.from_gtf("annotation.gtf", gene="BRCA1")
ps = tx.ProteinSpace.from_transcript(gs, "ENST00000357654")

fig = tx.Figure(
    tx.GenomicPanel(gs, [
        tx.Coverage(coverage_df),
        tx.Sashimi(junction_df, scale="log"),
        tx.Transcripts("annotation.gtf", "BRCA1", cds_utr=True, chevrons=True),
        tx.Lollipops(mutations_df, coord="genomic", color_by="category"),
    ], title="BRCA1"),
    tx.Bridge(gs, ps),
    tx.ProteinPanel(ps, [
        tx.Domains(pfam_df),
        tx.Lollipops(mutations_df, coord="protein", color_by="category"),
    ], title="BRCA1 protein"),
    figsize=(12, 7),
)
fig.save("brca1.svg")
```

## Design

Three coordinate objects:

| Object           | Maps                                    | Used by                         |
| ---------------- | --------------------------------------- | ------------------------------- |
| `GenomicSpace`   | genomic bp to [0, 1] figure coord       | `GenomicPanel` + genomic tracks |
| `ProteinSpace`   | 1-based AA position to [0, 1]           | `ProteinPanel` + protein tracks |
| `Bridge`         | genomic CDS exons to AA ranges          | Own renderable panel            |

`GenomicSpace` optionally compresses introns to a fixed fractional width per
intron; exonic regions share the remaining width proportionally. Disable with
`compress_introns=False` for a linear axis.

Tracks render into a horizontal slot of their panel's Axes. They have a
`height` attribute that determines their relative vertical allocation.

## Tracks

Each track accepts a `pandas.DataFrame` with a documented schema. No parsing
of BAM/VCF is done in v0.1; users bring their data as DataFrames and can
pipeline via pysam, bcbio, polars, etc.

### `Transcripts(gtf, gene, cds_utr=True, chevrons=True)`

Reads exon and CDS records from the GTF. Draws tall rectangles for CDS
regions and thin rectangles for UTR regions; overlays directional `>`/`<`
chevrons on intron lines.

### `Lollipops(data, coord="genomic" | "protein", color_by="category")`

DataFrame columns:

| Column     | Type   | Required | Description                                     |
| ---------- | ------ | -------- | ----------------------------------------------- |
| `position` | int    | yes      | Genomic (if `coord="genomic"`) or AA position   |
| `height`   | float  | no       | VAF or recurrence count (default 1)             |
| `category` | str    | no       | Consequence label for color                     |
| `label`    | str    | no       | Text drawn next to the head                     |

Consequence colors default to a colorblind-safe palette (missense=vermillion,
nonsense=black, synonymous=sky, splice=reddish purple, ...).

### `Sashimi(data, scale="linear" | "log")`

| Column   | Type | Description                             |
| -------- | ---- | --------------------------------------- |
| `start`  | int  | Junction 5' genomic position            |
| `end`    | int  | Junction 3' genomic position            |
| `count`  | int  | Read count supporting the junction      |

Renders quadratic-Bezier arcs between `start` and `end`; arc linewidth
proportional to `count` (with optional `log1p` scaling for wide dynamic range).

### `Coverage(data)`

| Column     | Type | Description         |
| ---------- | ---- | ------------------- |
| `position` | int  | Genomic position    |
| `depth`    | int  | Read depth          |

Filled histogram. The max depth is annotated at the right edge.

### `Domains(data)` (ProteinSpace only)

| Column  | Type | Description                              |
| ------- | ---- | ---------------------------------------- |
| `name`  | str  | Domain name, drawn as label in the box   |
| `start` | int  | AA start (1-based inclusive)             |
| `end`   | int  | AA end (inclusive)                       |

## Palettes

`txplot.palette.OKABE_ITO` is the default for small categorical sets (up to 8).
Larger sets delegate to [glasbey](https://github.com/lmcinnes/glasbey) via
`txplot.palette.categorical(n, colorblind_safe=True)`.

## What's not in v0.1

Deferred to v0.2+:

- Conservation track (phyloP / phastCons; needs bigWig support).
- Direct BAM/VCF parsing helpers (`txplot.io`).
- Protein hydrophobicity / secondary-structure overlays.
- Multi-sample / multi-condition faceting.
