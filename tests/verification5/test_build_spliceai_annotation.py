"""
Unit tests for scripts/build_spliceai_annotation.py

The annotation is what makes a locally computed SpliceAI score describe the transcript we
report rather than GENCODE v24's canonical. A silent omission here would reproduce exactly the
defect the work exists to remove, so the tests are mostly about what happens when a transcript
cannot be found or maps to more than one structure.
"""
import gzip
import importlib.util
import io
import os
import sys

import pytest

_spec = importlib.util.spec_from_file_location(
    "build_spliceai_annotation",
    os.path.join(os.path.dirname(__file__), "..", "..", "scripts",
                 "build_spliceai_annotation.py"))
bsa = importlib.util.module_from_spec(_spec)
sys.modules["build_spliceai_annotation"] = bsa
_spec.loader.exec_module(bsa)

HEADER = "#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n"


def row(gene, enst, canonical, biotype, refseqs, chrom="chr1", strand="+",
        start="100", end="900", ex_s="100,500,", ex_e="200,900,"):
    name = f"{gene}---ENSG0001---{enst}---{'yes' if canonical else 'no'}---{biotype}---{refseqs}"
    return "\t".join([name, chrom, strand, start, end, ex_s, ex_e]) + "\n"


def annotation(tmp_path, rows, gz=True):
    p = tmp_path / ("ann.txt.gz" if gz else "ann.txt")
    data = HEADER + "".join(rows)
    if gz:
        with gzip.open(p, "wt") as fh:
            fh.write(data)
    else:
        p.write_text(data)
    return str(p)


def build(tmp_path, rows, requested, gz=True):
    out = io.StringIO()
    rep = io.StringIO()
    written, missing, ambiguous = bsa.build(
        annotation(tmp_path, rows, gz), bsa.read_requested(requested), out, report=rep)
    return out.getvalue(), written, missing, ambiguous, rep.getvalue()


def wanted(tmp_path, *accessions):
    p = tmp_path / "want.txt"
    p.write_text("\n".join(accessions) + "\n")
    return str(p)


class TestParseName:
    def test_unpacks_the_six_fields(self):
        m = bsa.parse_name("BRAF---ENSG0001---ENST0002---yes---protein_coding---NM_004333.6")
        assert m["gene"] == "BRAF" and m["enst"] == "ENST0002"
        assert m["canonical"] is True and m["biotype"] == "protein_coding"
        assert m["refseqs"] == ["NM_004333.6"]

    def test_several_refseqs(self):
        m = bsa.parse_name("X---G---T---no---protein_coding---NM_1.1,NM_2.2")
        assert m["refseqs"] == ["NM_1.1", "NM_2.2"]

    def test_a_dot_is_not_an_accession(self):
        assert bsa.parse_name("X---G---T---no---protein_coding---.")["refseqs"] == []

    def test_a_name_without_the_six_fields_is_rejected(self):
        assert bsa.parse_name("BRAF") is None


class TestReadRequested:
    def test_keyed_without_the_version_but_remembers_it(self, tmp_path):
        got = bsa.read_requested(wanted(tmp_path, "NM_004333.6", "NM_000546.6"))
        assert got == {"NM_004333": "NM_004333.6", "NM_000546": "NM_000546.6"}

    def test_blank_lines_and_comments_are_ignored(self, tmp_path):
        p = tmp_path / "w.txt"
        p.write_text("# a comment\n\nNM_1.1\n   \n")
        assert bsa.read_requested(str(p)) == {"NM_1": "NM_1.1"}


class TestBuild:
    def test_writes_the_requested_transcript(self, tmp_path):
        out, written, missing, amb, _ = build(
            tmp_path, [row("BRAF", "ENST1", True, "protein_coding", "NM_004333.6")],
            wanted(tmp_path, "NM_004333.6"))
        assert written == 1 and missing == [] and amb == []
        assert out.startswith("#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n")
        assert "BRAF--NM_004333.6\tchr1\t+\t100\t900\t100,500,\t200,900," in out

    def test_the_name_carries_the_transcript(self, tmp_path):
        # SpliceAI reports NAME as its SYMBOL, which is how a score becomes traceable
        out, *_ = build(tmp_path, [row("MET", "E1", False, "protein_coding", "NM_001127500.3")],
                        wanted(tmp_path, "NM_001127500.3"))
        assert "MET--NM_001127500.3" in out

    def test_matching_ignores_the_version(self, tmp_path):
        out, written, missing, *_ = build(
            tmp_path, [row("BRAF", "E1", True, "protein_coding", "NM_004333.9")],
            wanted(tmp_path, "NM_004333.6"))
        assert written == 1 and missing == []
        assert "NM_004333.9" in out        # the annotation's spelling, not ours

    def test_unrequested_transcripts_are_left_out(self, tmp_path):
        out, written, *_ = build(
            tmp_path, [row("BRAF", "E1", True, "protein_coding", "NM_004333.6"),
                       row("TP53", "E2", True, "protein_coding", "NM_000546.6")],
            wanted(tmp_path, "NM_004333.6"))
        assert written == 1 and "TP53" not in out

    def test_a_transcript_not_in_the_annotation_is_reported(self, tmp_path):
        out, written, missing, amb, report = build(
            tmp_path, [row("BRAF", "E1", True, "protein_coding", "NM_004333.6")],
            wanted(tmp_path, "NM_004333.6", "NM_999999.1"))
        assert written == 1
        assert missing == ["NM_999999"]
        assert "MISSING NM_999999.1" in report   # named, never dropped quietly

    def test_two_structures_resolve_to_the_canonical(self, tmp_path):
        out, written, missing, amb, report = build(
            tmp_path, [row("G", "E_noncanon", False, "protein_coding", "NM_1.1", end="500"),
                       row("G", "E_canon", True, "protein_coding", "NM_1.1", end="900")],
            wanted(tmp_path, "NM_1.1"))
        assert written == 1
        assert amb and amb[0][2] == "E_canon"
        assert "\t900\t" in out
        assert "canonical" in report

    def test_two_structures_with_no_canonical_say_so(self, tmp_path):
        out, written, missing, amb, report = build(
            tmp_path, [row("G", "E_a", False, "protein_coding", "NM_1.1"),
                       row("G", "E_b", False, "protein_coding", "NM_1.1")],
            wanted(tmp_path, "NM_1.1"))
        assert written == 1 and amb
        assert amb[0][3] is False
        assert "none canonical" in report

    def test_one_row_per_transcript(self, tmp_path):
        # SpliceAI scores every annotation row overlapping the variant; one row per
        # transcript is what keeps the output unambiguous
        out, written, *_ = build(
            tmp_path, [row("CDKN2A", "E1", True, "protein_coding", "NM_000077.5"),
                       row("CDKN2A", "E2", False, "protein_coding", "NM_058195.4")],
            wanted(tmp_path, "NM_000077.5", "NM_058195.4"))
        assert written == 2
        assert "CDKN2A--NM_000077.5" in out and "CDKN2A--NM_058195.4" in out

    def test_a_plain_text_annotation_is_read_too(self, tmp_path):
        out, written, *_ = build(
            tmp_path, [row("BRAF", "E1", True, "protein_coding", "NM_004333.6")],
            wanted(tmp_path, "NM_004333.6"), gz=False)
        assert written == 1

    def test_exon_coordinates_are_copied_unchanged(self, tmp_path):
        # the point of reusing the published annotation: the structures are the ones the
        # published tool uses, so a score computed against them is comparable with it
        out, *_ = build(
            tmp_path, [row("G", "E1", True, "protein_coding", "NM_1.1",
                           ex_s="11,22,33,", ex_e="111,222,333,")],
            wanted(tmp_path, "NM_1.1"))
        assert "\t11,22,33,\t111,222,333," in out


class TestTheRealPanelAnnotation:
    """Against the committed panel list, if a SpliceAI-lookup deployment is present."""

    ANN = ("/Users/monkiky/Documents/external/refs/SpliceAI-lookup/references_files/"
           "annotations/gencode.v43.annotation.txt.gz")
    PANEL = os.path.join(os.path.dirname(__file__), "..", "..", "reference_lists",
                         "TSO500_transcripts_list.txt")

    @pytest.mark.skipif(not os.path.exists(ANN), reason="no SpliceAI-lookup annotation present")
    def test_both_cdkn2a_transcripts_are_addressable(self, tmp_path):
        out, written, missing, *_ = build_real(self.ANN, self.PANEL)
        assert "CDKN2A--NM_000077.5" in out    # p16INK4a
        assert "CDKN2A--NM_058195.4" in out    # p14ARF, the one we report
        assert written > 500


def build_real(ann, panel):
    out, rep = io.StringIO(), io.StringIO()
    written, missing, amb = bsa.build(ann, bsa.read_requested(panel), out, report=rep)
    return out.getvalue(), written, missing, amb, rep.getvalue()
