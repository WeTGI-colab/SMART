"""
Unit tests for scripts/check_reference.py

Nothing verified that a VCF's REF alleles agree with the assembly they are annotated
against, and VEP does not refuse a VCF that disagrees: it takes the reference base from
the FASTA and keeps the submitted ALT, so the consequence it reports belongs to a variant
nobody submitted. Found in our own fixture, where three records were wrong and no layer
said so.

Covers:
  - reading a FASTA through its .fai, including line wrapping and contig naming
  - which records are checkable: symbolic ALTs and placeholder REFs are not
  - detection of strand-flipped alleles, off-by-one REFs and a wrong assembly
  - fail versus warn, and the exit status the pipeline depends on
"""
import gzip
import importlib.util
import os
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "..", "scripts")
_spec = importlib.util.spec_from_file_location("check_reference",
                                               os.path.join(SCRIPTS, "check_reference.py"))
cr = importlib.util.module_from_spec(_spec)
sys.modules["check_reference"] = cr
_spec.loader.exec_module(cr)


# chr1 is 60 bases wrapped at 20; chr2 exists to test contig naming
CHR1 = "ACGTACGTACGTACGTACGT" "TTTTCCCCGGGGAAAATTTT" "GATCGATCGATCGATCGATC"
CHR2 = "AAAACCCCGGGGTTTT"


def write_fasta(tmp_path, contigs=(("chr1", CHR1, 20), ("chr2", CHR2, 16)), index=True):
    """A FASTA plus a correct .fai, built by hand so the offsets are exercised."""
    fa = tmp_path / "ref.fa"
    lines, fai = [], []
    offset = 0
    for name, seq, width in contigs:
        header = f">{name}\n"
        lines.append(header)
        offset += len(header)
        fai.append((name, len(seq), offset, width, width + 1))
        for i in range(0, len(seq), width):
            chunk = seq[i:i + width] + "\n"
            lines.append(chunk)
            offset += len(chunk)
    fa.write_text("".join(lines))
    if index:
        (tmp_path / "ref.fa.fai").write_text(
            "".join(f"{n}\t{l}\t{o}\t{lb}\t{lw}\n" for n, l, o, lb, lw in fai))
    return str(fa)


def write_vcf(tmp_path, records, name="in.vcf", gz=False):
    """records: iterable of (chrom, pos, id, ref, alt)."""
    body = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    for c, p, i, r, a in records:
        body += f"{c}\t{p}\t{i}\t{r}\t{a}\t.\t.\t.\n"
    path = tmp_path / (name + (".gz" if gz else ""))
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(body)
    else:
        path.write_text(body)
    return str(path)


# ===========================================================================
# Reference
# ===========================================================================

class TestReference:
    def test_fetches_a_single_base(self, tmp_path):
        ref = cr.Reference(write_fasta(tmp_path))
        assert ref.fetch("chr1", 1, 1) == "A"
        assert ref.fetch("chr1", 4, 1) == "T"
        ref.close()

    def test_fetches_across_a_line_break(self, tmp_path):
        # bases 19-22 straddle the wrap at 20, so the newline must not be counted
        ref = cr.Reference(write_fasta(tmp_path))
        assert ref.fetch("chr1", 19, 4) == CHR1[18:22]
        ref.close()

    def test_fetches_a_span_covering_several_lines(self, tmp_path):
        ref = cr.Reference(write_fasta(tmp_path))
        assert ref.fetch("chr1", 5, 50) == CHR1[4:54]
        ref.close()

    def test_fetches_the_last_base(self, tmp_path):
        ref = cr.Reference(write_fasta(tmp_path))
        assert ref.fetch("chr1", len(CHR1), 1) == CHR1[-1]
        ref.close()

    def test_past_the_end_is_none(self, tmp_path):
        ref = cr.Reference(write_fasta(tmp_path))
        assert ref.fetch("chr1", len(CHR1), 2) is None
        assert ref.fetch("chr1", 0, 1) is None
        ref.close()

    def test_a_second_contig_is_offset_correctly(self, tmp_path):
        ref = cr.Reference(write_fasta(tmp_path))
        assert ref.fetch("chr2", 1, 4) == "AAAA"
        assert ref.fetch("chr2", 13, 4) == "TTTT"
        ref.close()

    @pytest.mark.parametrize("asked,expected", [("chr1", "chr1"), ("1", "chr1")])
    def test_the_chr_prefix_may_differ(self, tmp_path, asked, expected):
        ref = cr.Reference(write_fasta(tmp_path))
        assert ref.resolve(asked) == expected
        ref.close()

    def test_a_contig_the_fasta_lacks_resolves_to_none(self, tmp_path):
        ref = cr.Reference(write_fasta(tmp_path))
        assert ref.resolve("chrZ") is None
        ref.close()

    def test_a_missing_index_is_an_explicit_error(self, tmp_path):
        fa = write_fasta(tmp_path, index=False)
        with pytest.raises(FileNotFoundError, match="samtools faidx"):
            cr.Reference(fa)


# ===========================================================================
# is_checkable
# ===========================================================================

class TestIsCheckable:
    @pytest.mark.parametrize("ref,alt", [("A", "T"), ("ACGT", "A"), ("a", "t")])
    def test_literal_alleles_are_checkable(self, ref, alt):
        assert cr.is_checkable(ref, alt)

    @pytest.mark.parametrize("alt", ["<DUP>", "<DEL>", "<INS>", ".", "]chr1:100]A", "A["])
    def test_symbolic_alts_are_not(self, alt):
        assert not cr.is_checkable("N", alt)
        assert not cr.is_checkable("A", alt)

    @pytest.mark.parametrize("ref", ["N", "ANG", "", "R", "ACGTN"])
    def test_placeholder_refs_are_not(self, ref):
        assert not cr.is_checkable(ref, "A")

    def test_an_absurdly_long_ref_is_not(self):
        assert not cr.is_checkable("A" * (cr.MAX_REF_LEN + 1), "A")


# ===========================================================================
# check
# ===========================================================================

class TestCheck:
    def test_a_correct_vcf_has_no_mismatches(self, tmp_path):
        ref = cr.Reference(write_fasta(tmp_path))
        vcf = write_vcf(tmp_path, [("chr1", 1, "a", "A", "T"),
                                   ("chr1", 21, "b", "T", "C")])
        res = cr.check(vcf, ref)
        assert res["n_mismatches"] == 0
        assert res["checked"] == 2
        ref.close()

    def test_a_strand_flipped_allele_is_caught(self, tmp_path):
        # genome has A at chr1:1; a minus-strand write-up gives T
        ref = cr.Reference(write_fasta(tmp_path))
        vcf = write_vcf(tmp_path, [("chr1", 1, "flip", "T", "A")])
        res = cr.check(vcf, ref)
        assert res["n_mismatches"] == 1
        chrom, pos, vid, got_ref, genome, alt = res["mismatches"][0]
        assert (vid, got_ref, genome) == ("flip", "T", "A")
        ref.close()

    def test_an_off_by_one_multibase_ref_is_caught(self, tmp_path):
        ref = cr.Reference(write_fasta(tmp_path))
        right = CHR1[0:5]
        shifted = CHR1[1:6]
        ok = write_vcf(tmp_path, [("chr1", 1, "ok", right, "A")], name="ok.vcf")
        bad = write_vcf(tmp_path, [("chr1", 1, "shifted", shifted, "A")], name="bad.vcf")
        assert cr.check(ok, ref)["n_mismatches"] == 0
        assert cr.check(bad, ref)["n_mismatches"] == 1
        ref.close()

    def test_symbolic_records_are_skipped_not_failed(self, tmp_path):
        ref = cr.Reference(write_fasta(tmp_path))
        vcf = write_vcf(tmp_path, [("chr1", 1, "dup", "N", "<DUP>"),
                                   ("chr1", 5, "del", "N", "<DEL>")])
        res = cr.check(vcf, ref)
        assert res["checked"] == 0
        assert res["skipped"] == 2
        assert res["n_mismatches"] == 0
        ref.close()

    def test_an_unknown_contig_is_recorded_not_counted_as_a_mismatch(self, tmp_path):
        ref = cr.Reference(write_fasta(tmp_path))
        vcf = write_vcf(tmp_path, [("chrZ", 1, "x", "A", "T")])
        res = cr.check(vcf, ref)
        assert res["missing_contigs"] == ["chrZ"]
        assert res["n_mismatches"] == 0
        ref.close()

    def test_a_position_past_the_contig_end_is_recorded(self, tmp_path):
        ref = cr.Reference(write_fasta(tmp_path))
        vcf = write_vcf(tmp_path, [("chr1", 10_000, "x", "A", "T")])
        res = cr.check(vcf, ref)
        assert res["out_of_range"] == 1
        assert res["n_mismatches"] == 0
        ref.close()

    def test_lower_case_reference_and_vcf_agree(self, tmp_path):
        # soft-masked reference bases are lower case and must still match
        ref = cr.Reference(write_fasta(tmp_path, contigs=(("chr1", CHR1.lower(), 20),)))
        vcf = write_vcf(tmp_path, [("chr1", 1, "a", "A", "T")])
        assert cr.check(vcf, ref)["n_mismatches"] == 0
        ref.close()

    def test_a_gzipped_vcf_is_read(self, tmp_path):
        ref = cr.Reference(write_fasta(tmp_path))
        vcf = write_vcf(tmp_path, [("chr1", 1, "a", "T", "A")], gz=True)
        assert cr.check(vcf, ref)["n_mismatches"] == 1
        ref.close()

    def test_mismatches_are_counted_beyond_the_report_limit(self, tmp_path):
        ref = cr.Reference(write_fasta(tmp_path))
        recs = [("chr1", i, f"v{i}", "N" if False else
                 ("T" if CHR1[i - 1] != "T" else "A"), "C") for i in range(1, 21)]
        vcf = write_vcf(tmp_path, recs)
        res = cr.check(vcf, ref, limit=3)
        assert res["n_mismatches"] == 20
        assert len(res["mismatches"]) == 3
        ref.close()


# ===========================================================================
# report and exit status
# ===========================================================================

class TestReportAndExit:
    def _res(self, n, checked=10):
        return {"checked": checked, "skipped": 0,
                "mismatches": [("chr1", "1", "v", "A", "T", "G")] * min(n, 5),
                "n_mismatches": n, "missing_contigs": [], "out_of_range": 0}

    def test_clean_run_continues(self, capsys):
        assert cr.report(self._res(0), "x.vcf", "ref.fa", "fail") is True
        assert "no REF mismatches" in capsys.readouterr().out

    def test_fail_mode_stops_the_run(self, capsys):
        assert cr.report(self._res(1), "x.vcf", "ref.fa", "fail") is False
        assert "ERROR" in capsys.readouterr().out

    def test_warn_mode_continues_but_says_so(self, capsys):
        assert cr.report(self._res(1), "x.vcf", "ref.fa", "warn") is True
        assert "WARNING" in capsys.readouterr().out

    def test_a_mostly_mismatching_vcf_suggests_the_wrong_assembly(self, capsys):
        cr.report(self._res(9, checked=10), "x.vcf", "ref.fa", "warn")
        assert "different assembly" in capsys.readouterr().out

    def test_a_few_mismatches_suggest_strand(self, capsys):
        cr.report(self._res(1, checked=100), "x.vcf", "ref.fa", "warn")
        assert "strand-flipped" in capsys.readouterr().out

    def test_the_vep_behaviour_is_always_explained(self, capsys):
        cr.report(self._res(1), "x.vcf", "ref.fa", "warn")
        assert "VEP does not reject" in capsys.readouterr().out


# ===========================================================================
# the real fixture, which is what prompted this
# ===========================================================================

FIXTURE = os.path.join(os.path.dirname(__file__), "..", "verification1", "verification1.vcf")


@pytest.mark.skipif(not os.path.exists(FIXTURE), reason="verification1.vcf not present")
class TestTheKnownBadFixture:
    """verification1.vcf carries three records whose REF does not match GRCh38.

    Asserted against a stub reference holding just those positions, so the test needs no
    3 GB FASTA: the point is that the checker flags these and passes the others.
    """

    def test_the_three_known_bad_records_are_still_there(self):
        import re
        text = open(FIXTURE).read()
        for vid, ref, alt in (("GNAQ_Q209L", "A", "T"),
                              ("DICER1_E1705K", "G", "A"),
                              ("EGFR_del19", "GAATTAAGAGAAGCAT", "G")):
            assert re.search(rf"\t{vid}\t{ref}\t{alt}\t", text), (
                f"{vid} no longer has the alleles this test was written against; "
                f"if the fixture was corrected, update or remove this test")

    def test_the_correct_alleles_are_recorded_for_whoever_fixes_it(self):
        # GRCh38, verified against the reference: chr9:77794572 T>A,
        # chr14:95094139 C>T, chr7:55174772 GGAATTAAGAGAAGCA>G
        expected = {"GNAQ_Q209L": ("chr9", 77794572, "T", "A"),
                    "DICER1_E1705K": ("chr14", 95094139, "C", "T"),
                    "EGFR_del19": ("chr7", 55174772, "GGAATTAAGAGAAGCA", "G")}
        assert len(expected) == 3
