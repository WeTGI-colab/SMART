"""
Unit tests for tests/verify_local_sources.py

The verifier compares SMART's annotation against the files it was annotated from, closing
the gap left by verification1 (which can only reach sources with a REST API). Its own
comparison logic is what these tests cover -- the normalisation in particular, because a
naive comparison produced 350 phantom mismatches before it was right.
"""
import gzip
import importlib.util
import os
import sys

import pytest

HERE = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    "verify_local_sources", os.path.join(HERE, "..", "verify_local_sources.py"))
vls = importlib.util.module_from_spec(_spec)
sys.modules["verify_local_sources"] = vls
_spec.loader.exec_module(vls)


class TestContig:
    @pytest.mark.parametrize("given,style,expected", [
        ("chr1", "chr", "chr1"), ("1", "chr", "chr1"),
        ("chr1", "plain", "1"), ("1", "plain", "1"),
        ("chrX", "plain", "X"), ("X", "chr", "chrX"),
    ])
    def test_both_naming_conventions(self, given, style, expected):
        # ClinVar, SpliceAI, REVEL and dbNSFP use '1'; Cancer Hotspots uses 'chr1'
        assert vls.contig(given, style) == expected


class TestNorm:
    @pytest.mark.parametrize("value", ["", ".", "-", "nan", "NA", "None", "[]"])
    def test_the_empty_forms_are_one_thing(self, value):
        assert vls.norm(value) == ""

    def test_numbers_compare_as_numbers(self):
        assert vls.norm("0.10") == vls.norm("0.1") == "0.1"

    @pytest.mark.parametrize("value", ["True", "true", "yes", "present"])
    def test_the_unambiguous_spellings_of_presence_fold_together(self, value):
        assert vls.norm(value) == "True"

    @pytest.mark.parametrize("value", ["1", "1.0"])
    def test_one_stays_a_number_here(self, value):
        # 1 is both how the table writes a flag's presence and a real score phastCons
        # reaches. norm() cannot know which, so it keeps the number and verdict() decides
        # the flag case where both sides are in view. See TestFlagsAreNotNumbers.
        assert vls.norm(value) == "1"

    def test_list_separators_are_equivalent(self):
        # VEP rewrites ',' and '|' to '&' inside a custom-annotation value
        assert vls.norm("a,b,c") == vls.norm("a|b|c") == vls.norm("a&b&c")

    def test_list_order_does_not_matter(self):
        assert vls.norm("b&a") == vls.norm("a&b")


class TestSameNumber:
    def test_three_and_four_significant_figures_agree(self):
        # the VEP cache writes 1.16e-05 where the REST API writes 1.159e-05
        assert vls._same_number("1.16e-05", "1.159e-05")

    def test_a_real_difference_does_not(self):
        assert not vls._same_number("0.08", "0.21")

    def test_zero_against_zero(self):
        assert vls._same_number("0", "0.0")

    def test_text_is_not_a_number(self):
        assert not vls._same_number("Pathogenic", "0.1")


class TestGnomadKey:
    @pytest.mark.parametrize("col,key", [
        ("gnomADe_AF", "gnomade"),
        ("gnomADg_AF", "gnomadg"),
        ("gnomADe_NFE_AF", "gnomade_nfe"),
        ("gnomADg_REMAINING_AF", "gnomadg_remaining"),
        ("gnomADe_AFR_AF", "gnomade_afr"),
        # the 1000 Genomes populations carry no prefix
        ("AFR_AF", "afr"), ("EUR_AF", "eur"), ("SAS_AF", "sas"),
    ])
    def test_maps_a_column_to_the_api_key(self, col, key):
        assert vls.gnomad_key(col) == key


class TestVerdict:
    def test_equal_values_pass(self):
        assert vls.verdict("Pathogenic", "Pathogenic")[0] == "PASS"

    def test_both_empty_passes(self):
        assert vls.verdict("", set())[0] == "PASS"
        assert vls.verdict("nan", "")[0] == "PASS"

    def test_a_value_the_source_lacks_mismatches(self):
        assert vls.verdict("0.931", set())[0] == "MISMATCH"

    def test_any_of_several_source_rows_passes(self):
        res, why = vls.verdict("0.5", {"0.5", "0.7"})
        assert res == "PASS" and "2 source rows" in why

    def test_rounding_passes_and_says_so(self):
        res, why = vls.verdict("1.16e-05", "1.159e-05")
        assert res == "PASS" and "rounding" in why

    def test_a_different_value_mismatches(self):
        assert vls.verdict("0.08", "0.21")[0] == "MISMATCH"


class TestFetchFrequencies:
    def test_indels_are_skipped_rather_than_queried(self, tmp_path):
        # the region endpoint takes substitutions; an indel needs a different call
        got, note = vls.fetch_frequencies("chr7", 55174772, "GGAATTAAGAGAAGCA", "G",
                                          str(tmp_path))
        assert got is None and "substitutions only" in note

    def test_symbolic_alts_are_skipped(self, tmp_path):
        got, note = vls.fetch_frequencies("chr7", 1, "N", "<DUP>", str(tmp_path))
        assert got is None and note.startswith("SKIP")

    def test_a_cached_response_is_reused_without_a_request(self, tmp_path):
        import json
        d = tmp_path / "cache"
        d.mkdir()
        (d / "1_100_A_T.json").write_text(json.dumps(
            [{"colocated_variants": [{"id": "rs1", "frequencies": {"T": {"gnomade": 0.5}}}]}]))
        got, note = vls.fetch_frequencies("chr1", 100, "A", "T", str(d))
        assert note is None and got == {"gnomade": 0.5}

    def test_a_variant_with_no_frequencies_is_answered_not_errored(self, tmp_path):
        import json
        d = tmp_path / "cache"
        d.mkdir()
        (d / "1_100_A_T.json").write_text(json.dumps(
            [{"colocated_variants": [{"id": "COSV1", "somatic": 1}]}]))
        got, note = vls.fetch_frequencies("chr1", 100, "A", "T", str(d))
        assert note is None and got == {}


class TestSources:
    def test_every_declared_source_has_a_contig_style(self):
        for name, spec in vls.SOURCES.items():
            assert "path" in spec
            assert "contig" in spec, f"{name} does not declare its contig naming"

    def test_clinvar_and_hotspots_disagree_on_the_chr_prefix(self):
        # the thing most likely to make a lookup silently find nothing
        assert vls.SOURCES["clinvar"]["contig"] == "plain"
        assert vls.SOURCES["hotspots"]["contig"] == "chr"


class TestSeparatorsIncludeSemicolon:
    def test_dbnsfp_packs_values_with_a_semicolon(self):
        # dbNSFP writes one value per transcript joined by ';' and VEP rewrites it to '&';
        # without this the whole dbNSFP module reported every multi-transcript row as wrong
        assert vls.norm("T;T;T") == vls.norm("T&T&T")
        assert vls.norm(".;LB;B") == vls.norm(".&LB&B")

    def test_all_four_separators_are_one_thing(self):
        assert vls.norm("a;b") == vls.norm("a|b") == vls.norm("a,b") == vls.norm("a&b")


class TestFlagsAreNotNumbers:
    def test_a_flag_matches_the_one_the_table_writes(self):
        # CancerHotspots_HOTSPOT is a VCF flag; the table renders its presence as 1.0
        assert vls.verdict("1.0", "True")[0] == "PASS"
        assert vls.verdict("1", "True")[0] == "PASS"

    def test_a_score_of_one_stays_a_number(self):
        # phastCons reaches 1.0 often, and it is a measurement, not a presence
        assert vls.norm("1.0") == "1"
        assert vls.verdict("1.0", "1.000000")[0] == "PASS"

    def test_a_flag_does_not_match_an_unrelated_number(self):
        assert vls.verdict("0.5", "True")[0] == "MISMATCH"


class TestGenie:
    def test_haematological_lineages_are_the_classifier_s(self):
        # must stay the same set genie_annotation.py uses, or the recomputation is not a check
        assert vls.GENIE_HAEM == {"Myeloid", "Lymphoid"}

    def _counts(self, tmp_path, rows):
        import gzip as gz
        d = tmp_path / "GENIE"
        d.mkdir()
        hdr = ("gene\tvariant_class\tchrom\tstart\tend\tref\talt\tref_aa\talt_aa\t"
               "oncotree_code\ttissue\tlineage\tn_patients\n")
        with gz.open(d / "genie_o4_counts.tsv.gz", "wt") as fh:
            fh.write(hdr)
            for r in rows:
                fh.write("\t".join(r) + "\n")
        return vls.load_genie_counts(str(tmp_path))

    def test_counts_are_summed_across_rows_of_one_variant(self, tmp_path):
        counts = self._counts(tmp_path, [
            ("BRAF", "missense", "chr7", "140753336", "140753336", "A", "T",
             "V", "E", "COAD", "Bowel", "Solid", "1401"),
            ("BRAF", "missense", "chr7", "140753336", "140753336", "A", "T",
             "V", "E", "LUAD", "Lung", "Solid", "464"),
        ])
        row = dict(Chromosome="chr7", Start_Position=140753336,
                   Reference_Allele="A", Tumor_Seq_Allele2="T")
        got = dict((c, v) for c, v, _ in vls.check_genie(
            row, ["GENIE_SAME_NUCLEOTIDE_ALL", "GENIE_SAME_NUCLEOTIDE_SOLID",
                  "GENIE_SAME_NUCLEOTIDE_HAEM"], counts))
        assert got["GENIE_SAME_NUCLEOTIDE_ALL"] == "1865"
        assert got["GENIE_SAME_NUCLEOTIDE_SOLID"] == "1865"
        assert got["GENIE_SAME_NUCLEOTIDE_HAEM"] == "0"

    def test_a_patient_with_two_lineages_counts_in_both(self, tmp_path):
        counts = self._counts(tmp_path, [
            ("X", "missense", "chr1", "100", "100", "A", "T", "V", "E",
             "A|B", "t1|t2", "Lymphoid|Solid", "3"),
        ])
        row = dict(Chromosome="chr1", Start_Position=100,
                   Reference_Allele="A", Tumor_Seq_Allele2="T")
        got = dict((c, v) for c, v, _ in vls.check_genie(
            row, ["GENIE_SAME_NUCLEOTIDE_ALL", "GENIE_SAME_NUCLEOTIDE_HAEM",
                  "GENIE_SAME_NUCLEOTIDE_SOLID"], counts))
        assert got == {"GENIE_SAME_NUCLEOTIDE_ALL": "3",
                       "GENIE_SAME_NUCLEOTIDE_HAEM": "3",
                       "GENIE_SAME_NUCLEOTIDE_SOLID": "3"}

    def test_a_variant_absent_from_the_table_counts_zero(self, tmp_path):
        counts = self._counts(tmp_path, [])
        row = dict(Chromosome="chr1", Start_Position=999,
                   Reference_Allele="A", Tumor_Seq_Allele2="T")
        got = dict((c, v) for c, v, _ in vls.check_genie(
            row, ["GENIE_SAME_NUCLEOTIDE_ALL"], counts))
        assert got["GENIE_SAME_NUCLEOTIDE_ALL"] == "0"

    def test_indels_are_skipped_rather_than_recomputed(self, tmp_path):
        # matching an indel needs the left-alignment the pipeline applies; an independent
        # attempt would differ for reasons of its own rather than finding a defect
        counts = self._counts(tmp_path, [])
        row = dict(Chromosome="chr13", Start_Position=48342703,
                   Reference_Allele="CTA", Tumor_Seq_Allele2="C")
        notes = [n for _, _, n in vls.check_genie(
            row, ["GENIE_SAME_NUCLEOTIDE_ALL"], counts)]
        assert notes and "left-alignment" in notes[0]

    def test_the_derived_o4_columns_are_not_claimed_to_be_recomputed(self, tmp_path):
        counts = self._counts(tmp_path, [])
        row = dict(Chromosome="chr1", Start_Position=1,
                   Reference_Allele="A", Tumor_Seq_Allele2="T")
        notes = {c: n for c, _, n in vls.check_genie(
            row, ["GENIE_O4_MATCH_TYPE", "GENIE_TOP_TUMOUR_TYPES"], counts)}
        for c, n in notes.items():
            assert n and n.startswith("SKIP"), f"{c} must not claim an independent check"

    def test_a_subset_count_above_the_whole_is_caught(self, tmp_path):
        counts = self._counts(tmp_path, [])
        row = dict(Chromosome="chr1", Start_Position=1, Reference_Allele="A",
                   Tumor_Seq_Allele2="T", GENIE_O4_MATCH_ALL="10",
                   GENIE_O4_MATCH_HAEM="12")
        got = dict((c, v) for c, v, _ in vls.check_genie(
            row, ["GENIE_O4_MATCH_HAEM"], counts))
        assert got["GENIE_O4_MATCH_HAEM"] == "<= 10.0"

    def test_an_unannotated_row_is_not_held_to_the_invariant(self, tmp_path):
        counts = self._counts(tmp_path, [])
        row = dict(Chromosome="chr1", Start_Position=1, Reference_Allele="A",
                   Tumor_Seq_Allele2="T", GENIE_O4_MATCH_ALL="nan",
                   GENIE_O4_MATCH_HAEM="nan")
        notes = [n for _, _, n in vls.check_genie(row, ["GENIE_O4_MATCH_HAEM"], counts)]
        assert notes and notes[0].startswith("SKIP")

    def test_a_missing_counts_table_is_skipped_not_failed(self, tmp_path):
        row = dict(Chromosome="chr1", Start_Position=1,
                   Reference_Allele="A", Tumor_Seq_Allele2="T")
        notes = [n for _, _, n in vls.check_genie(row, ["GENIE_SAME_NUCLEOTIDE_ALL"], None)]
        assert notes and "absent" in notes[0]


class TestPatcherCreatesAbsentColumns:
    """In compute mode VEP runs without the SpliceAI plugin, so the columns do not exist.

    The patcher was first written for the A/B measurement, where a lookup-mode run had already
    created them, and it refused a MAF without them. That made the whole compute path fail at
    the merge step — loudly, but completely. These tests pin the created-column behaviour so it
    cannot regress back to replace-only.
    """

    @staticmethod
    def _patcher():
        import importlib.util
        import os
        import sys
        spec = importlib.util.spec_from_file_location(
            "spliceai_patch_maf",
            os.path.join(os.path.dirname(__file__), "..", "..", "scripts",
                         "spliceai_patch_maf.py"))
        m = importlib.util.module_from_spec(spec)
        sys.modules["spliceai_patch_maf"] = m
        spec.loader.exec_module(m)
        return m

    def _files(self, tmp_path, maf_header, maf_row, info):
        maf = tmp_path / "in.maf"
        maf.write_text("#SMART_VERSION 2.1.0\n" + "\t".join(maf_header) + "\n"
                       + "\t".join(maf_row) + "\n")
        vcf = tmp_path / "in.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n"
                       "##contig=<ID=chr7,length=159345973>\n"
                       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                       f"chr7\t140753336\t.\tA\tT\t.\t.\t{info}\n")
        return str(maf), str(vcf)

    HDR = ["Hugo_Symbol", "Chromosome", "Start_Position",
           "Reference_Allele", "Tumor_Seq_Allele2"]
    ROW = ["BRAF", "chr7", "140753336", "A", "T"]
    INFO = "SpliceAI=T|BRAF--NM_004333.6|0.01|0.02|0.03|0.04|1|2|3|4"

    def test_a_maf_without_the_columns_gains_them(self, tmp_path):
        m = self._patcher()
        maf, vcf = self._files(tmp_path, self.HDR, self.ROW, self.INFO)
        local, _ = m.read_local(vcf)
        out = tmp_path / "out.maf"
        import io as _io
        m.patch(maf, local, str(out), report=_io.StringIO())
        lines = out.read_text().splitlines()
        hdr = lines[1].split("\t")
        for c in ("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_DL",
                  "SpliceAI_pred_DP_AG", "SpliceAI_pred_SYMBOL"):
            assert c in hdr, f"{c} was not created"
        row = dict(zip(hdr, lines[2].split("\t")))
        assert row["SpliceAI_pred_DS_AG"] == "0.01"
        assert row["SpliceAI_pred_DS_DL"] == "0.04"
        assert row["SpliceAI_pred_SYMBOL"] == "BRAF--NM_004333.6"

    def test_the_leading_metadata_line_survives(self, tmp_path):
        m = self._patcher()
        maf, vcf = self._files(tmp_path, self.HDR, self.ROW, self.INFO)
        local, _ = m.read_local(vcf)
        out = tmp_path / "out.maf"
        import io as _io
        m.patch(maf, local, str(out), report=_io.StringIO())
        assert out.read_text().splitlines()[0] == "#SMART_VERSION 2.1.0"

    def test_existing_columns_are_replaced_not_duplicated(self, tmp_path):
        m = self._patcher()
        hdr = self.HDR + ["SpliceAI_pred_DS_AG"]
        maf, vcf = self._files(tmp_path, hdr, self.ROW + ["0.99"], self.INFO)
        local, _ = m.read_local(vcf)
        out = tmp_path / "out.maf"
        import io as _io
        m.patch(maf, local, str(out), report=_io.StringIO())
        lines = out.read_text().splitlines()
        cols = lines[1].split("\t")
        assert cols.count("SpliceAI_pred_DS_AG") == 1
        assert dict(zip(cols, lines[2].split("\t")))["SpliceAI_pred_DS_AG"] == "0.01"


class TestComputeModeSpliceAI:
    """Compute mode has no precomputed source, so the checks are internal consistency.

    The module used to skip these rows entirely, which took it from 137 checks to none and
    stopped looking at the one thing that can go wrong: `spliceai_patch_maf.py` merges the
    model's VCF into the table on CHROM+POS+REF+ALT, and a key matching the wrong row would
    put one variant's deltas on another.
    """

    @staticmethod
    def _row(**kw):
        base = {"Hugo_Symbol": "NRAS", "NM_Transcript": "NM_002524",
                "SpliceAI_pred_SYMBOL": "NRAS--NM_002524.5",
                "SpliceAI_pred_DS_AG": "0.12", "SpliceAI_pred_DP_AG": "-34"}
        base.update(kw)
        return base

    def _verdicts(self, row, cols):
        out = {}
        for col, src, note in vls.check_spliceai_computed(row, cols):
            if note:
                out[col] = note.split(":")[0]
            else:
                out[col] = vls.verdict(row.get(col, ""), src)[0]
        return out

    def test_compute_mode_is_recognised_from_the_symbol(self):
        """`GENE--REFSEQ` can only come from a local run: the tables hold a bare symbol."""
        assert vls.computed_locally(self._row())
        assert not vls.computed_locally(self._row(SpliceAI_pred_SYMBOL="NRAS"))
        assert not vls.computed_locally(self._row(SpliceAI_pred_SYMBOL=""))

    def test_a_matching_transcript_passes(self):
        v = self._verdicts(self._row(), ["SpliceAI_pred_SYMBOL"])
        assert v["SpliceAI_pred_SYMBOL"] == "PASS"

    def test_the_version_suffix_is_not_what_is_being_tested(self):
        row = self._row(NM_Transcript="NM_002524.4")
        v = self._verdicts(row, ["SpliceAI_pred_SYMBOL"])
        assert v["SpliceAI_pred_SYMBOL"] == "PASS"

    def test_another_genes_scores_on_this_row_are_caught(self):
        """The misaligned-merge failure, which is why this check exists."""
        row = self._row(SpliceAI_pred_SYMBOL="KRAS--NM_004985.5")
        v = self._verdicts(row, ["SpliceAI_pred_SYMBOL"])
        assert v["SpliceAI_pred_SYMBOL"] == "MISMATCH"

    def test_the_right_gene_but_the_wrong_transcript_is_caught(self):
        row = self._row(SpliceAI_pred_SYMBOL="NRAS--NM_999999.1")
        v = self._verdicts(row, ["SpliceAI_pred_SYMBOL"])
        assert v["SpliceAI_pred_SYMBOL"] == "MISMATCH"

    def test_a_row_with_no_reported_transcript_cannot_pass_vacuously(self):
        row = self._row(NM_Transcript="")
        v = self._verdicts(row, ["SpliceAI_pred_SYMBOL"])
        assert v["SpliceAI_pred_SYMBOL"] == "MISMATCH"

    @pytest.mark.parametrize("ds", ["0", "0.5", "1", "1.0"])
    def test_a_delta_in_range_passes(self, ds):
        v = self._verdicts(self._row(SpliceAI_pred_DS_AG=ds), ["SpliceAI_pred_DS_AG"])
        assert v["SpliceAI_pred_DS_AG"] == "PASS"

    @pytest.mark.parametrize("ds", ["1.3", "-0.2", "12"])
    def test_a_delta_outside_zero_to_one_is_caught(self, ds):
        """A delta is a probability; anything else is a corrupted or mis-parsed field."""
        v = self._verdicts(self._row(SpliceAI_pred_DS_AG=ds), ["SpliceAI_pred_DS_AG"])
        assert v["SpliceAI_pred_DS_AG"] == "MISMATCH"

    def test_a_non_numeric_delta_is_caught(self):
        v = self._verdicts(self._row(SpliceAI_pred_DS_AG="D"), ["SpliceAI_pred_DS_AG"])
        assert v["SpliceAI_pred_DS_AG"] == "MISMATCH"

    @pytest.mark.parametrize("dp", ["0", "-500", "500", "106"])
    def test_a_position_within_the_window_passes(self, dp):
        v = self._verdicts(self._row(SpliceAI_pred_DP_AG=dp), ["SpliceAI_pred_DP_AG"])
        assert v["SpliceAI_pred_DP_AG"] == "PASS"

    def test_a_position_beyond_the_window_is_caught(self):
        """The pipeline runs with -D 500, so 900 cannot have come from it."""
        v = self._verdicts(self._row(SpliceAI_pred_DP_AG="900"), ["SpliceAI_pred_DP_AG"])
        assert v["SpliceAI_pred_DP_AG"] == "MISMATCH"

    @pytest.mark.parametrize("zero", [0.0, 0, "0.0", "0", "0.00"])
    def test_a_delta_of_exactly_zero_is_a_score_not_a_blank(self, zero):
        """The most common value SpliceAI reports, and the one a falsy test loses.

        Written after `str(row.get(col, "") or "")` turned every 0.0 into "" -- 0.0 is falsy
        in Python -- so 95 real checks were reported as "the model returned no score". The
        same shape of mistake as treating NaN as a value: the test has to distinguish a
        value that happens to be zero from no value at all.
        """
        v = self._verdicts(self._row(SpliceAI_pred_DS_AG=zero), ["SpliceAI_pred_DS_AG"])
        assert v["SpliceAI_pred_DS_AG"] == "PASS", "a delta of zero was not checked"

    @pytest.mark.parametrize("zero", [0, "0", 0.0])
    def test_a_delta_position_of_zero_is_a_position_not_a_blank(self, zero):
        v = self._verdicts(self._row(SpliceAI_pred_DP_AG=zero), ["SpliceAI_pred_DP_AG"])
        assert v["SpliceAI_pred_DP_AG"] == "PASS"

    def test_a_none_cell_is_still_a_skip(self):
        """None means the column was absent, which is different from a zero."""
        v = self._verdicts(self._row(SpliceAI_pred_DS_AG=None), ["SpliceAI_pred_DS_AG"])
        assert v["SpliceAI_pred_DS_AG"] == "SKIP"

    def test_an_unscored_row_is_skipped_not_failed(self):
        """The model scores no MNV and no indel; the patcher leaves those cells empty."""
        row = self._row(SpliceAI_pred_DS_AG="", SpliceAI_pred_DP_AG="")
        v = self._verdicts(row, ["SpliceAI_pred_DS_AG", "SpliceAI_pred_DP_AG"])
        assert v["SpliceAI_pred_DS_AG"] == "SKIP"
        assert v["SpliceAI_pred_DP_AG"] == "SKIP"

    def test_the_cutoff_is_still_not_a_source_value(self):
        v = self._verdicts(self._row(), ["SpliceAI_cutoff"])
        assert v["SpliceAI_cutoff"] == "SKIP"


class TestLoeufIsTranscriptSpecific:
    """LOEUF is described everywhere as gene-level. In the dataset it is per transcript.

    Measured on loeuf_dataset_grch38.tsv.gz, counting only the ENST rows: of the 14,247 genes
    with more than one transcript, 13,236 (93%) carry a different LOEUF on different
    transcripts. KAT6B is 0.057 on ENST00000650232 (CDS 2,974 bp, 14 coding exons) and 1.891
    on ENST00000649119 (CDS 23 bp, 1 exon) -- a short transcript has few possible
    loss-of-function sites, so the expected count is small and the upper bound of the
    observed/expected interval, which is what LOEUF is, comes out large.

    So looking LOEUF up by gene returns an arbitrary transcript's value and is wrong for 93%
    of multi-transcript genes, sometimes by the distance between "among the most constrained
    genes in the genome" and "no constraint signal at all". These tests exist so that
    `match_by=transcript`, and the by-transcript lookup in check_loeuf, cannot be simplified
    away as redundant.
    """

    REFS = os.environ.get("SMART_REFS",
                          os.path.expanduser("~/Documents/external/refs"))
    LOEUF = os.path.join(REFS, "gnomAD_constraints", "loeuf_dataset_grch38.tsv.gz")

    def test_the_lookup_is_keyed_on_transcript_not_gene(self):
        """A static check, so it holds with no reference data present."""
        src = open(os.path.join(os.path.dirname(__file__), "..",
                                "verify_local_sources.py")).read()
        start = src.index("def load_loeuf(")
        body = src[start:src.index("\ndef ", start + 10)]
        assert '"transcript"' in body, (
            "load_loeuf no longer keys on the transcript column; a gene-keyed table would "
            "return an arbitrary transcript's LOEUF")
        start = src.index("def check_loeuf(")
        body = src[start:src.index("\ndef ", start + 10)]
        assert 'row.get("Feature"' in body, (
            "check_loeuf no longer looks up by the row's reported transcript")

    @pytest.mark.skipif(not os.path.exists(LOEUF),
                        reason="LOEUF dataset not present (set SMART_REFS)")
    def test_the_dataset_really_does_vary_by_transcript(self):
        """The fact itself, asserted against the file rather than taken from a comment."""
        import collections
        per = collections.defaultdict(set)
        with gzip.open(self.LOEUF, "rt") as fh:
            hdr = fh.readline().rstrip("\n").lstrip("#").split("\t")
            gi, ti, li = (hdr.index("gene"), hdr.index("transcript"),
                          hdr.index("lof.oe_ci.upper"))
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= li or not f[ti].startswith("ENST"):
                    continue
                try:
                    per[f[gi]].add(round(float(f[li]), 4))
                except ValueError:
                    continue
        multi = {g: v for g, v in per.items() if len(v) >= 1}
        varying = [g for g, v in per.items() if len(v) > 1]
        assert len(multi) > 10000, f"only {len(multi)} genes read; is the file truncated?"
        assert len(varying) > 5000, (
            f"only {len(varying)} genes have a transcript-dependent LOEUF. If this has "
            f"genuinely changed, the by-transcript lookup may no longer be necessary -- but "
            f"check the file before relaxing it")

    @pytest.mark.skipif(not os.path.exists(LOEUF),
                        reason="LOEUF dataset not present (set SMART_REFS)")
    def test_a_named_gene_spans_the_whole_constraint_range(self):
        """KAT6B, so the test names a case a reader can look up and check."""
        vals = {}
        with gzip.open(self.LOEUF, "rt") as fh:
            hdr = fh.readline().rstrip("\n").lstrip("#").split("\t")
            gi, ti, li = (hdr.index("gene"), hdr.index("transcript"),
                          hdr.index("lof.oe_ci.upper"))
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > li and f[gi] == "KAT6B" and f[ti].startswith("ENST"):
                    try:
                        vals[f[ti]] = float(f[li])
                    except ValueError:
                        pass
        assert len(vals) > 1, f"KAT6B should have several transcripts, found {sorted(vals)}"
        assert max(vals.values()) - min(vals.values()) > 1.0, (
            f"KAT6B's LOEUF no longer spans more than 1.0 across transcripts: {vals}")


class TestRetryBackoff:
    """An ERROR is a service that did not answer, not a disagreement, so it is worth retrying.

    On the 37-variant cohort the previous four attempts still lost 52 checks to Ensembl
    returning 500 under load, and those counts varied run to run — 130, 26, 52 — which made
    the same cohort look differently verified on different days.
    """

    def test_the_wait_doubles_and_is_capped(self):
        waits = [vls._backoff(a) for a in range(8)]
        assert waits[0] < waits[1] < waits[2], "the wait does not grow"
        assert all(w <= 30.5 for w in waits), f"an uncapped wait would stall a run: {waits}"

    def test_the_wait_carries_jitter(self):
        """Without it, several clients backing off together resynchronise on every round."""
        seen = {round(vls._backoff(3), 4) for _ in range(40)}
        assert len(seen) > 1, "the backoff is deterministic, so clients will synchronise"

    def test_six_attempts_are_configured(self):
        assert vls.VEP_ATTEMPTS >= 6

    def test_the_error_message_states_how_many_tries(self, tmp_path, monkeypatch):
        """So a reader can tell a brief outage from a service that is simply down."""
        import urllib.request
        def boom(*a, **kw):
            raise OSError("down")
        monkeypatch.setattr(urllib.request, "urlopen", boom)
        monkeypatch.setattr(vls.time, "sleep", lambda *_: None)
        got, note = vls.fetch_frequencies("chr1", 100, "A", "G", str(tmp_path), sleep=0)
        assert got is None
        assert f"after {vls.VEP_ATTEMPTS} tries" in note, note

    def test_a_permanent_rejection_is_not_retried(self, tmp_path, monkeypatch):
        """A 400 is an answer: the service understood and has nothing. Retrying wastes a run."""
        import urllib.error, urllib.request
        calls = []
        def reject(*a, **kw):
            calls.append(1)
            raise urllib.error.HTTPError("u", 400, "Bad Request", {}, None)
        monkeypatch.setattr(urllib.request, "urlopen", reject)
        monkeypatch.setattr(vls.time, "sleep", lambda *_: None)
        got, note = vls.fetch_frequencies("chr1", 100, "A", "G", str(tmp_path), sleep=0)
        assert got is None and "rejected" in note
        assert len(calls) == 1, f"a 400 was retried {len(calls)} times"

    def test_a_cached_success_is_never_refetched(self, tmp_path, monkeypatch):
        """The cache is what makes patience affordable: a success is paid for once."""
        import json as _json, urllib.request
        cache = tmp_path
        (cache / "1_100_A_G.json").write_text(_json.dumps([{"colocated_variants": []}]))
        def boom(*a, **kw):
            raise AssertionError("the network was used despite a cached response")
        monkeypatch.setattr(urllib.request, "urlopen", boom)
        vls.fetch_frequencies("chr1", 100, "A", "G", str(cache), sleep=0)


class TestNonAsciiIsNoLongerFolded:
    """While issue 27 was open, the comparison folded both sides to ASCII, because the
    annotator deleted those characters and a difference was unavoidable. The image now
    removes that line, so folding would compare two mutilated values and find them equal --
    the same shape of mistake as ONCOKB_VUS, where verifier and pipeline shared a typo and
    the check passed for years on two empty strings.
    """

    @pytest.mark.parametrize("value", ["PI3Kβ", "α-ketoglutarate", "≥6 months",
                                       "25%–56%", "Müller"])
    def test_characters_are_preserved(self, value):
        assert vls.ascii_fold(value) == value, "a non-ASCII character was dropped"

    def test_a_difference_in_those_characters_is_now_caught(self):
        """The point of removing the fold: this had to pass before and must fail now."""
        assert vls.verdict("DOR of ≥6 months", {"DOR of 6 months"})[0] == "MISMATCH"
        assert vls.verdict("PI3Kβ", {"PI3K"})[0] == "MISMATCH"

    def test_identical_values_still_agree(self):
        assert vls.verdict("PI3Kβ", {"PI3Kβ"})[0] == "PASS"


class TestEveryCheckNamesItsVariant:
    """`local_source_results.tsv` had `Variant_ID = "."` on all 5,986 rows.

    The expression was `row.get("Variant_ID") or row.get("dbSNP_RS") or "."` and neither
    column exists in the MAF, so the fallback was the only branch that ever ran. 4,996
    passing checks were no comfort: the file exists to be read when something disagrees,
    and there was nothing in it to say which variant had disagreed.
    """

    def test_the_maf_id_is_used(self):
        assert vls.variant_label({"ID": "TERT_C228T"}) == "TERT_C228T"

    def test_neither_absent_column_is_required(self):
        """Both columns the old expression looked for are absent from the MAF."""
        row = {"ID": "NRAS_Q61R", "Chromosome": "chr1", "Start_Position": "114713908"}
        assert "Variant_ID" not in row and "dbSNP_RS" not in row
        assert vls.variant_label(row) == "NRAS_Q61R"

    def test_it_falls_back_to_coordinates_rather_than_a_placeholder(self):
        row = {"Chromosome": "chr17", "Start_Position": "7675994",
               "Reference_Allele": "G", "Tumor_Seq_Allele2": "A"}
        assert vls.variant_label(row) == "chr17:7675994_G>A"

    def test_coordinates_without_alleles_still_identify_the_site(self):
        assert vls.variant_label({"Chromosome": "chrX", "Start_Position": "101356177"}) \
            == "chrX:101356177"

    @pytest.mark.parametrize("placeholder", [".", "-", "nan", "", "   "])
    def test_a_placeholder_id_does_not_win_over_coordinates(self, placeholder):
        """A VCF writing "." in ID is the normal DRAGEN case, and was the reason this
        produced "." for every row."""
        row = {"ID": placeholder, "Chromosome": "chr7", "Start_Position": "140753336",
               "Reference_Allele": "A", "Tumor_Seq_Allele2": "T"}
        assert vls.variant_label(row) == "chr7:140753336_A>T"

    def test_dbsnp_is_still_honoured_when_it_is_what_the_input_has(self):
        assert vls.variant_label({"dbSNP_RS": "rs113488022"}) == "rs113488022"

    def test_a_row_with_nothing_at_all_is_the_only_case_that_yields_a_dot(self):
        assert vls.variant_label({}) == "."

    def test_the_committed_results_file_names_its_variants(self):
        """The regression this guards: one distinct Variant_ID across thousands of rows."""
        import csv as _csv
        path = os.path.join(os.path.dirname(__file__), "..", "verification1",
                            "local_source_results.tsv")
        if not os.path.exists(path):
            pytest.skip("localsources has not been run")
        _csv.field_size_limit(10 ** 7)
        with open(path, encoding="utf-8") as fh:
            ids = {r["Variant_ID"] for r in _csv.DictReader(fh, delimiter="\t")}
        assert ids - {"."}, "every row still says '.': no check can be traced to a variant"
        assert len(ids) > 1, f"only {len(ids)} distinct variant label(s) across the file"
