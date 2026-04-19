# txview

Terminal plotter for transcript models of a gene from a GTF annotation.

## Install

Latest from `main`:

```bash
pip install "txview @ git+https://github.com/ymvk/monorepo@main#subdirectory=packages/txview"
```

Pin to a version tag:

```bash
pip install "txview @ git+https://github.com/ymvk/monorepo@txview-v0.1.0#subdirectory=packages/txview"
```

## CLI

```bash
txview path/to/annotation.gtf BRCA1
```

Accepts `.gtf` or `.gtf.gz`. The gene may be specified by `gene_name` or
`gene_id`.

Flags:

| Flag              | Default                     | Description                                 |
| ----------------- | --------------------------- | ------------------------------------------- |
| `--width`         | terminal width - 30         | Display width in columns                    |
| `--intron-width`  | 5                           | Columns reserved per truncated intron       |

## Example

```
BRCA1 [chr17:43044295-43125483] (-)
ENST00000357654  │████████████████████████████╌╌╌╌╌█╌╌╌╌╌█╌╌╌╌╌█╌╌╌╌╌██╌╌╌╌╌██╌╌╌╌╌██╌╌╌╌╌███████████████████████████◀
ENST00000471181  │████████████████████████████╌╌╌╌╌─╌╌╌╌╌─╌╌╌╌╌─╌╌╌╌╌██╌╌╌╌╌──╌╌╌╌╌──╌╌╌╌╌███████████████████████◀
ENST00000586385                                    │╌╌╌╌╌─╌╌╌╌╌█╌╌╌╌╌──╌╌╌╌╌██╌╌╌╌╌──╌╌╌╌╌█◀
```

Legend:

| Char | Meaning                                                                   |
| ---- | ------------------------------------------------------------------------- |
| `█`  | Exon                                                                      |
| `─`  | Intron within a transcript's span (e.g. alternative splicing)             |
| `╌`  | Truncated intron region (intergenic-like compression for display)         |
| `▶`  | TSS (5' end) on `+` strand; `◀` on `-` strand                             |
| `│`  | TTS (3' end)                                                              |

## Python API

```python
from txview import parse_gtf, layout, render

transcripts = parse_gtf("annotation.gtf", "BRCA1")
segments, width = layout(transcripts, total_width=120, intron_width=5)
print(render(transcripts, segments, width))
```

`parse_gtf` returns `dict[str, Transcript]` keyed by `transcript_id`.
