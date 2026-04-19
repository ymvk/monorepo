# monorepo

Personal monorepo for packages, scripts, and learning material.

## Packages

### [txview](https://github.com/ymvk/monorepo/tree/main/packages/txview) (2026-04-19)

Terminal plotter for transcript models of a gene from a GTF annotation. Parses the GTF for the requested gene (by name or ID), compresses introns to a fixed visual width so long genes stay on one line, and renders each transcript with Unicode blocks for exons plus strand-aware TSS/TTS markers.

### [txplot](https://github.com/ymvk/monorepo/tree/main/packages/txplot) (2026-04-19)

Matplotlib-based library for rendering transcript model figures. Treats genomic and protein coordinate systems as first-class composable objects, so joint genomic-protein views with "exon bridge" polygons fall out of the API rather than requiring bolt-on logic. Ships tracks for transcripts (CDS/UTR with directional chevrons), sashimi-style splice arcs, coverage, mutation lollipops (plottable in either genomic or protein space), and Pfam-style domains, with colorblind-safe Okabe-Ito and glasbey palettes by default.
