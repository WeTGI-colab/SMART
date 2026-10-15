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
class TestTheCorrectedFixture:
    """verification1.vcf carried three records whose REF did not match GRCh38 (#19).

    All three are now corrected, verified twice over: the REF against the GRCh38 FASTA,
    and the resulting variant against OncoKB's byGenomicChange endpoint, which resolves
    each corrected record to the alteration the fixture name claims -- Q209L, E1705K and
    E746_A750del -- and resolves none of the originals to anything.

    The test now guards the correction rather than the defect, so a regenerated fixture
    cannot quietly reintroduce it.
    """

    # GRCh38. EGFR carries the anchor base, making it the canonical 15 bp in-frame
    # deletion c.2236_2250del; the original wrote the deleted bases without an anchor
    # and got the last one wrong. GNAQ and DICER1 are both minus-strand genes whose
    # alleles had been written on the gene's strand instead of the genome's.
    CORRECTED = {"GNAQ_Q209L":    ("chr9",  77794572, "T", "A"),
                 "DICER1_E1705K": ("chr14", 95094139, "C", "T"),
                 "EGFR_del19":    ("chr7",  55174772, "GGAATTAAGAGAAGCA", "G")}

    # The same three records as they were before #19 was fixed, used as a positive control
    # below: the checker must flag exactly these against GRCh38 and nothing else.
    PRE_FIX = {"GNAQ_Q209L":    ("chr9",  77794572, "A", "T"),
               "DICER1_E1705K": ("chr14", 95094139, "G", "A"),
               "EGFR_del19":    ("chr7",  55174772, "GAATTAAGAGAAGCAT", "G")}

    def test_the_fixture_carries_the_corrected_alleles(self):
        import re
        text = open(FIXTURE).read()
        for vid, (chrom, pos, ref, alt) in self.CORRECTED.items():
            assert re.search(rf"^{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t", text, re.M), (
                f"{vid} does not carry the corrected alleles {ref}>{alt} at {chrom}:{pos}; "
                f"if the fixture was regenerated, it has reintroduced issue #19")

    def test_the_strand_flipped_alleles_are_gone(self):
        """The specific wrong values, named, so a partial revert is caught too."""
        import re
        text = open(FIXTURE).read()
        for vid, ref, alt in (("GNAQ_Q209L", "A", "T"),
                              ("DICER1_E1705K", "G", "A"),
                              ("EGFR_del19", "GAATTAAGAGAAGCAT", "G")):
            assert not re.search(rf"\t{vid}\t{ref}\t{alt}\t", text), (
                f"{vid} is back to the pre-#19 alleles {ref}>{alt}")

    # Overridable, because the reference set does not live in the repository and its
    # location differs per machine. Skipped rather than failed where it is absent.
    HG38 = os.environ.get("SMART_HG38",
                          os.path.expanduser("~/Documents/external/refs/liftover/hg38.fa"))

    @pytest.mark.skipif(not os.path.exists(
        os.environ.get("SMART_HG38",
                       os.path.expanduser("~/Documents/external/refs/liftover/hg38.fa")) + ".fai"),
        reason="GRCh38 FASTA not present on this machine (set SMART_HG38 to point at it)")
    def test_the_two_committed_copies_agree(self):
        """The fixture's second defect in #19: `.vcf` and `.vcf.gz` had different contents.

        The pipeline reads the `.gz`, so the plain file was documentation that disagreed with
        what was tested -- it held 18 records against the gz's 22, missing EGFR_T790M,
        KRAS_G12C, PIK3CA_E545K and TP53_R248W, which is where the wrong figure "18 variants"
        came from. Nothing checked this, and regenerating one from the other is exactly how
        four variants silently stop being tested.
        """
        plain = [l for l in open(FIXTURE).read().splitlines() if not l.startswith("#")]
        with gzip.open(FIXTURE + ".gz", "rt") as fh:
            gz = [l for l in fh.read().splitlines() if not l.startswith("#")]
        assert plain == gz, (
            f"verification1.vcf has {len(plain)} records and verification1.vcf.gz has "
            f"{len(gz)}, or their contents differ; the pipeline reads the .gz, so they must "
            f"be generated from one another and committed together")
        assert len(gz) == 45, f"the fixture should hold 45 records, found {len(gz)}"

    @pytest.mark.skipif(not os.path.exists(
        os.environ.get("SMART_HG38",
                       os.path.expanduser("~/Documents/external/refs/liftover/hg38.fa")) + ".fai"),
        reason="GRCh38 FASTA not present on this machine (set SMART_HG38 to point at it)")
    def test_the_checker_detects_the_pre_fix_alleles_against_the_real_assembly(self):
        """The positive control, and the reason the test above means anything.

        A checker that returned "no mismatches" unconditionally would pass every other test
        here that uses the real FASTA. So feed it the three alleles as they were before #19
        was fixed and require it to flag all three, reading the same assembly. The detection
        tests elsewhere in this file use a stub FASTA; this one uses GRCh38 itself.
        """
        import re
        text = open(FIXTURE).read()
        for vid, (chrom, pos, ref, alt) in self.PRE_FIX.items():
            text = re.sub(rf"^{chrom}\t{pos}\t{vid}\t[^\t]+\t[^\t]+",
                          f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}", text, flags=re.M)
        import tempfile
        with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as fh:
            fh.write(text)
            broken = fh.name
        try:
            ref = cr.Reference(self.HG38)
            try:
                res = cr.check(broken, ref)
            finally:
                ref.close()
        finally:
            os.unlink(broken)
        flagged = {m[2] for m in res["mismatches"]}
        assert flagged == set(self.PRE_FIX), (
            f"the checker flagged {sorted(flagged)} against GRCh38; it must flag all three "
            f"of {sorted(self.PRE_FIX)}, or it is not reading the assembly")
        # and it must not flag anything else: the other 15 records are correct
        assert res["n_mismatches"] == 3, f"also flagged records that are correct: {res['mismatches']}"

    def test_the_checker_passes_the_fixture_against_the_real_assembly(self):
        """The end-to-end assertion: run the real checker over the real FASTA.

        The other tests in this file compare the fixture against values recorded in this
        file, so they would both be wrong together if the recorded values were. This one
        reads the assembly, and it is the reason the fixture can be trusted. It skips
        where the FASTA is absent rather than passing vacuously.
        """
        ref = cr.Reference(self.HG38)
        try:
            res = cr.check(FIXTURE, ref)
        finally:
            ref.close()
        assert res["n_mismatches"] == 0, (
            f"REF mismatches against GRCh38: {res['mismatches']}")
        assert res["checked"] == 41, (
            f"{res['checked']} records were checked; the fixture holds 45, of which 4 have "
            f"symbolic ALTs, so 41 should be checkable -- a different number means something "
            f"is being skipped that should not be")
