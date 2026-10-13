"""
Unit tests for tests/verify_local_sources.py

The verifier compares SMART's annotation against the files it was annotated from, closing
the gap left by verification1 (which can only reach sources with a REST API). Its own
comparison logic is what these tests cover -- the normalisation in particular, because a
naive comparison produced 350 phantom mismatches before it was right.
"""
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
        got, note = vls.fetch_frequencies("chr7", 55174772, "GAATTAAGAGAAGCAT", "G",
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
