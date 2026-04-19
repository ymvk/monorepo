import textwrap

import matplotlib

matplotlib.use("Agg")

import pandas as pd
import pytest

from txplot import (
    Bridge,
    Coverage,
    Domains,
    Figure,
    GenomicPanel,
    GenomicSpace,
    Lollipops,
    ProteinPanel,
    ProteinSpace,
    Sashimi,
    Transcripts,
)
from txplot.tracks import _split_exon_by_cds

TINY_GTF = textwrap.dedent('''\
chr1\thavana\texon\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T1";
chr1\thavana\tCDS\t150\t200\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T1";
chr1\thavana\texon\t300\t400\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T1";
chr1\thavana\tCDS\t300\t400\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T1";
chr1\thavana\texon\t500\t600\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T1";
chr1\thavana\tCDS\t500\t550\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T1";
chr1\thavana\texon\t150\t200\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T2";
chr1\thavana\tCDS\t160\t200\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T2";
chr1\thavana\texon\t500\t600\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T2";
chr1\thavana\tCDS\t500\t540\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T2";
''')


@pytest.fixture
def gtf_path(tmp_path):
    p = tmp_path / "tiny.gtf"
    p.write_text(TINY_GTF)
    return str(p)


# ---------------------------------------------------------------------------
# GenomicSpace
# ---------------------------------------------------------------------------


def test_genomic_space_linear():
    gs = GenomicSpace("chr1", 100, 600, compress_introns=False)
    assert gs.transform(100) == 0.0
    assert gs.transform(600) == 1.0
    # Monotonic.
    xs = [gs.transform(p) for p in range(100, 601, 50)]
    assert all(a <= b for a, b in zip(xs, xs[1:]))


def test_genomic_space_compressed_ends_are_bounds():
    gs = GenomicSpace(
        "chr1", 100, 600,
        exon_regions=[(100, 200), (300, 400), (500, 600)],
        compress_introns=True,
        intron_width_frac=0.05,
    )
    assert gs.transform(100) == pytest.approx(0.0)
    assert gs.transform(600) == pytest.approx(1.0)


def test_genomic_space_from_gtf(gtf_path):
    gs = GenomicSpace.from_gtf(gtf_path, gene="FOO")
    assert gs.seqname == "chr1"
    assert gs.start == 100
    assert gs.end == 600
    assert "T1" in gs.transcripts and "T2" in gs.transcripts


# ---------------------------------------------------------------------------
# ProteinSpace
# ---------------------------------------------------------------------------


def test_protein_space_from_transcript(gtf_path):
    gs = GenomicSpace.from_gtf(gtf_path, gene="FOO")
    ps = ProteinSpace.from_transcript(gs, "T1")
    # Total CDS bp for T1 = (200-150+1) + (400-300+1) + (550-500+1) = 51 + 101 + 51 = 203
    # AA length = 203 // 3 = 67
    assert ps.length == 67
    assert ps.transform(1) == pytest.approx(0.0)
    assert ps.transform(ps.length) == pytest.approx(1.0)


def test_protein_space_missing_cds_errors(gtf_path, tmp_path):
    # Transcript with no CDS should raise.
    only_exons = 'chr1\thavana\texon\t100\t200\t.\t+\t.\tgene_id "G2"; gene_name "BAR"; transcript_id "T3";\n'
    p = tmp_path / "no_cds.gtf"
    p.write_text(only_exons)
    gs = GenomicSpace.from_gtf(str(p), gene="BAR")
    with pytest.raises(ValueError):
        ProteinSpace.from_transcript(gs, "T3")


# ---------------------------------------------------------------------------
# CDS/UTR split
# ---------------------------------------------------------------------------


def test_split_exon_all_utr():
    utr, cds = _split_exon_by_cds(100, 200, [])
    assert utr == [(100, 200)]
    assert cds == []


def test_split_exon_all_cds():
    utr, cds = _split_exon_by_cds(100, 200, [(100, 200)])
    assert utr == []
    assert cds == [(100, 200)]


def test_split_exon_mixed():
    utr, cds = _split_exon_by_cds(100, 200, [(150, 200)])
    assert utr == [(100, 149)]
    assert cds == [(150, 200)]


def test_split_exon_cds_in_middle():
    utr, cds = _split_exon_by_cds(100, 200, [(120, 180)])
    assert utr == [(100, 119), (181, 200)]
    assert cds == [(120, 180)]


# ---------------------------------------------------------------------------
# Integration: render a full figure without errors
# ---------------------------------------------------------------------------


def test_figure_renders(gtf_path, tmp_path):
    gs = GenomicSpace.from_gtf(gtf_path, gene="FOO")
    ps = ProteinSpace.from_transcript(gs, "T1")

    mutations = pd.DataFrame({
        "position": [180, 350, 520],
        "height": [3, 7, 2],
        "category": ["missense", "nonsense", "synonymous"],
        "label": ["p.R60Q", "p.W100*", ""],
    })
    junctions = pd.DataFrame({
        "start": [200, 400],
        "end": [300, 500],
        "count": [120, 45],
    })
    coverage = pd.DataFrame({
        "position": list(range(100, 601, 5)),
        "depth": [i % 30 + 5 for i in range(100, 601, 5)],
    })
    aa_muts = pd.DataFrame({
        "position": [10, 30, 55],
        "height": [2, 5, 1],
        "category": ["missense", "nonsense", "missense"],
    })
    domains = pd.DataFrame({
        "name": ["RING", "BRCT"],
        "start": [5, 40],
        "end": [25, 65],
    })

    fig = Figure(
        GenomicPanel(gs, [
            Coverage(coverage),
            Sashimi(junctions, scale="log"),
            Transcripts(gtf_path, "FOO"),
            Lollipops(mutations, coord="genomic", color_by="category"),
        ], title="FOO"),
        Bridge(gs, ps),
        ProteinPanel(ps, [
            Domains(domains),
            Lollipops(aa_muts, coord="protein", color_by="category"),
        ], title="T1 protein"),
        figsize=(12, 7),
    )
    out = tmp_path / "fig.svg"
    fig.save(str(out))
    assert out.exists() and out.stat().st_size > 0


def test_transcripts_track_no_cds_still_renders(gtf_path, tmp_path):
    gs = GenomicSpace.from_gtf(gtf_path, gene="FOO")
    fig = Figure(
        GenomicPanel(gs, [Transcripts(gtf_path, "FOO", cds_utr=False)]),
        figsize=(8, 2),
    )
    out = tmp_path / "tx.png"
    fig.save(str(out))
    assert out.exists()
