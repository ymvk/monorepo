import textwrap

from txview import (
    EXON_CHAR,
    TRUNC_CHAR,
    TSS_MINUS_CHAR,
    TSS_PLUS_CHAR,
    TTS_CHAR,
    layout,
    parse_gtf,
    render,
)

TINY_GTF = textwrap.dedent('''\
chr1\thavana\texon\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T1";
chr1\thavana\texon\t300\t400\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T1";
chr1\thavana\texon\t500\t600\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T1";
chr1\thavana\texon\t150\t200\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T2";
chr1\thavana\texon\t500\t600\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T2";
chr1\thavana\texon\t1000\t2000\t.\t+\t.\tgene_id "G2"; gene_name "BAR"; transcript_id "T3";
''')


def _write(tmp_path, body=TINY_GTF):
    gtf = tmp_path / "x.gtf"
    gtf.write_text(body)
    return gtf


def test_parse_filters_by_gene_name(tmp_path):
    ts = parse_gtf(_write(tmp_path), "FOO")
    assert set(ts.keys()) == {"T1", "T2"}
    assert ts["T1"].exons == [(100, 200), (300, 400), (500, 600)]
    assert ts["T2"].exons == [(150, 200), (500, 600)]
    assert ts["T1"].strand == "+"
    assert ts["T1"].seqname == "chr1"


def test_parse_filters_by_gene_id(tmp_path):
    ts = parse_gtf(_write(tmp_path), "G2")
    assert set(ts.keys()) == {"T3"}


def test_parse_unknown_gene_returns_empty(tmp_path):
    ts = parse_gtf(_write(tmp_path), "NOPE")
    assert ts == {}


def test_parse_ignores_comments_and_non_exon(tmp_path):
    body = (
        "# header comment\n"
        'chr1\thavana\tgene\t100\t600\t.\t+\t.\tgene_id "G1"; gene_name "FOO";\n'
        'chr1\thavana\ttranscript\t100\t600\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T1";\n'
        'chr1\thavana\texon\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T1";\n'
    )
    ts = parse_gtf(_write(tmp_path, body), "FOO")
    assert ts["T1"].exons == [(100, 200)]


def test_layout_produces_exonic_and_truncated_segments(tmp_path):
    ts = parse_gtf(_write(tmp_path), "FOO")
    segs, total = layout(ts, total_width=40, intron_width=3)
    exonic = [s for s in segs if s.exonic]
    introns = [s for s in segs if not s.exonic]
    # Merged exons for FOO: (100,200), (300,400), (500,600) → 3 exonic, 2 introns
    assert len(exonic) == 3
    assert len(introns) == 2
    assert all(s.dend - s.dstart == 3 for s in introns)
    # Segments are contiguous.
    for a, b in zip(segs, segs[1:]):
        assert a.dend == b.dstart
    assert total == segs[-1].dend


def test_render_contains_transcript_ids_and_blocks(tmp_path):
    ts = parse_gtf(_write(tmp_path), "FOO")
    segs, total = layout(ts, total_width=60, intron_width=3)
    out = render(ts, segs, total)
    assert "T1" in out
    assert "T2" in out
    assert "FOO" in out
    assert EXON_CHAR in out
    assert TRUNC_CHAR in out


def test_render_empty():
    assert render({}, [], 0) == "(no transcripts)"


def test_render_marks_tss_and_tts_plus_strand(tmp_path):
    ts = parse_gtf(_write(tmp_path), "FOO")
    segs, total = layout(ts, total_width=60, intron_width=3)
    lines = render(ts, segs, total).splitlines()
    t1_line = next(l for l in lines if l.startswith("T1"))
    # + strand: arrow at 5' (left), bar at 3' (right)
    assert TSS_PLUS_CHAR in t1_line
    assert TTS_CHAR in t1_line
    assert t1_line.index(TSS_PLUS_CHAR) < t1_line.index(TTS_CHAR)


def test_render_marks_tss_and_tts_minus_strand(tmp_path):
    body = (
        'chr1\thavana\texon\t100\t200\t.\t-\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T1";\n'
        'chr1\thavana\texon\t500\t600\t.\t-\t.\tgene_id "G1"; gene_name "FOO"; transcript_id "T1";\n'
    )
    ts = parse_gtf(_write(tmp_path, body), "FOO")
    segs, total = layout(ts, total_width=60, intron_width=3)
    lines = render(ts, segs, total).splitlines()
    t1_line = next(l for l in lines if l.startswith("T1"))
    # - strand: bar at 3' (left, lower coord), arrow at 5' (right, higher coord)
    assert TSS_MINUS_CHAR in t1_line
    assert TTS_CHAR in t1_line
    assert t1_line.index(TTS_CHAR) < t1_line.index(TSS_MINUS_CHAR)
