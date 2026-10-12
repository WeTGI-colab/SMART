"""
Unit tests for the annotation flags in scripts/post_analysis.py.

Issues: github.com/WeTGI-colab/SMART/issues/14 (SpliceAI evidence on LOW/MODIFIER rows)
        github.com/WeTGI-colab/SMART/issues/15 (the TERT promoter hotspot positions)

A consequence- or impact-based filter removes Silent/Intron/5'Flank rows before anyone
reaches the SpliceAI columns, so the evidence is present and the row is gone. Two columns
fix that without touching Variant_Classification, which is a defined MAF vocabulary:
SPLICEAI_MAX_DS carries the evidence as one sortable number, ANNOTATION_FLAGS carries the
reason to keep the row.

Covers:
  - _max_score()              — one SpliceAI cell to a number, and the empty forms
  - load_canonical_genomic()  — genomic entries only, read from the classifier's own list
  - add_annotation_flags()    — the 0.2 threshold, the LOW/MODIFIER scope, TERT, both at once
  - what reaches the MAF      — float dtype, empty cells, Variant_Classification untouched
"""
import io
import os

import pandas as pd
import pytest

from conftest import post

TERT_C228T = dict(Chromosome="chr5", Start_Position=1295113,
                  Reference_Allele="G", Tumor_Seq_Allele2="A")
TERT_C250T = dict(Chromosome="chr5", Start_Position=1295135,
                  Reference_Allele="G", Tumor_Seq_Allele2="A")

DS = ("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
      "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")


def row(impact="MODIFIER", scores=(".", ".", ".", "."), **kw):
    """One MAF-like row; scores are the four SpliceAI deltas in DS order."""
    d = dict(Chromosome="chr7", Start_Position=100, Reference_Allele="C",
             Tumor_Seq_Allele2="T", IMPACT=impact)
    d.update(dict(zip(DS, scores)))
    d.update(kw)
    return d


def flags(*rows, canonical=None):
    df = post.add_annotation_flags(pd.DataFrame(list(rows)), canonical=canonical)
    return df


# ===========================================================================
# _max_score
# ===========================================================================

class TestMaxScore:
    @pytest.mark.parametrize("value,expected", [
        ("0.95", 0.95),
        ("0", 0.0),
        ("1", 1.0),
        ("0.0123", 0.0123),
        ("1e-3", 0.001),
    ])
    def test_reads_a_single_score(self, value, expected):
        assert post._max_score(value) == pytest.approx(expected)

    @pytest.mark.parametrize("value", ["", ".", "-", "nan", "NaN", "None", "NA", "   "])
    def test_the_empty_forms_are_none(self, value):
        assert post._max_score(value) is None

    def test_a_cell_with_several_predictions_keeps_the_highest(self):
        # SpliceAI writes one value per field, but a cell that carried more than one
        # must not lose its evidence by failing to parse as a single float
        assert post._max_score("0.11&0.95&0.40") == pytest.approx(0.95)

    def test_text_without_a_number_is_none(self):
        assert post._max_score("no_prediction") is None


# ===========================================================================
# load_canonical_genomic
# ===========================================================================

class TestLoadCanonicalGenomic:
    def test_reads_the_two_tert_positions_from_the_repository_list(self):
        got = post.load_canonical_genomic()
        assert ("5", "1295113", "G", "A") in got
        assert ("5", "1295135", "G", "A") in got

    def test_protein_matched_entries_are_excluded(self):
        # the list is mostly protein-matched; a flag can only test a position
        assert len(post.load_canonical_genomic()) < 10

    def test_the_chr_prefix_is_stripped(self, tmp_path):
        p = tmp_path / "canon.tsv"
        p.write_text("gene\tmatch_rule\nTERT\tgenomic:chr5:1295113:G:A\n")
        assert post.load_canonical_genomic(str(p)) == {("5", "1295113", "G", "A")}

    def test_a_missing_file_is_a_warning_not_a_crash(self, tmp_path, capsys):
        got = post.load_canonical_genomic(str(tmp_path / "absent.tsv"))
        assert got == set()
        assert "not found" in capsys.readouterr().out

    def test_a_malformed_rule_is_skipped_and_reported(self, tmp_path, capsys):
        p = tmp_path / "canon.tsv"
        p.write_text("gene\tmatch_rule\nX\tgenomic:chr5:1295113\nTERT\tgenomic:chr5:1:G:A\n")
        got = post.load_canonical_genomic(str(p))
        assert got == {("5", "1", "G", "A")}
        assert "malformed" in capsys.readouterr().out

    def test_entries_extend_automatically(self, tmp_path):
        # the point of reading the list rather than hard-coding two positions
        p = tmp_path / "canon.tsv"
        p.write_text("gene\tmatch_rule\nTERT\tgenomic:chr5:1295113:G:A\n"
                     "FOO\tgenomic:chr9:12345:C:T\nBAR\tprotein\n")
        assert post.load_canonical_genomic(str(p)) == {("5", "1295113", "G", "A"),
                                                      ("9", "12345", "C", "T")}


# ===========================================================================
# SPLICEAI_MAX_DS
# ===========================================================================

class TestSpliceaiMaxDs:
    def test_takes_the_highest_of_the_four_deltas(self):
        df = flags(row(scores=("0.11", "0.95", "0.40", "0.02")))
        assert df["SPLICEAI_MAX_DS"].iloc[0] == pytest.approx(0.95)

    def test_ignores_the_empty_forms(self):
        df = flags(row(scores=(".", "0.33", "NA", "")))
        assert df["SPLICEAI_MAX_DS"].iloc[0] == pytest.approx(0.33)

    def test_is_empty_when_no_delta_has_a_value(self):
        df = flags(row(scores=(".", ".", ".", ".")))
        assert pd.isna(df["SPLICEAI_MAX_DS"].iloc[0])

    def test_is_computed_for_high_impact_too(self):
        # the flag is scoped to LOW/MODIFIER; the score is for everyone
        df = flags(row(impact="HIGH", scores=("0.99", ".", ".", ".")))
        assert df["SPLICEAI_MAX_DS"].iloc[0] == pytest.approx(0.99)

    def test_is_float_so_it_can_be_sorted_and_compared(self):
        df = flags(row(scores=("0.95", ".", ".", ".")), row(scores=(".", ".", ".", ".")))
        assert df["SPLICEAI_MAX_DS"].dtype == "float64"

    def test_missing_spliceai_columns_warn_rather_than_crash(self, capsys):
        df = post.add_annotation_flags(
            pd.DataFrame([dict(Chromosome="chr7", Start_Position=100,
                               Reference_Allele="C", Tumor_Seq_Allele2="T", IMPACT="LOW")]),
            canonical=set())
        assert pd.isna(df["SPLICEAI_MAX_DS"].iloc[0])
        assert "no SpliceAI delta columns" in capsys.readouterr().out


# ===========================================================================
# the splice_predicted flag
# ===========================================================================

class TestSplicePredictedFlag:
    @pytest.mark.parametrize("impact", ["LOW", "MODIFIER"])
    def test_set_for_low_and_modifier_above_the_threshold(self, impact):
        df = flags(row(impact=impact, scores=("0.95", ".", ".", ".")), canonical=set())
        assert df["ANNOTATION_FLAGS"].iloc[0] == "splice_predicted(0.95)"

    def test_set_exactly_at_the_threshold(self):
        # 0.2 is where the classifier starts awarding O6, so 0.2 must reach it
        df = flags(row(scores=("0.2", ".", ".", ".")), canonical=set())
        assert df["ANNOTATION_FLAGS"].iloc[0] == "splice_predicted(0.2)"

    def test_not_set_just_below_the_threshold(self):
        df = flags(row(scores=("0.19", ".", ".", ".")), canonical=set())
        assert df["ANNOTATION_FLAGS"].iloc[0] == ""

    @pytest.mark.parametrize("impact", ["HIGH", "MODERATE"])
    def test_not_set_for_high_and_moderate(self, impact):
        # an impact filter keeps these anyway; flagging them would flag most of the sheet
        df = flags(row(impact=impact, scores=("0.99", ".", ".", ".")), canonical=set())
        assert df["ANNOTATION_FLAGS"].iloc[0] == ""

    def test_not_set_without_a_score(self):
        df = flags(row(impact="LOW"), canonical=set())
        assert df["ANNOTATION_FLAGS"].iloc[0] == ""

    def test_impact_case_does_not_matter(self):
        df = flags(row(impact="low", scores=("0.95", ".", ".", ".")), canonical=set())
        assert "splice_predicted" in df["ANNOTATION_FLAGS"].iloc[0]

    def test_the_score_is_carried_in_the_flag(self):
        df = flags(row(scores=(".", "0.77", ".", ".")), canonical=set())
        assert df["ANNOTATION_FLAGS"].iloc[0] == "splice_predicted(0.77)"


# ===========================================================================
# the svig_canonical flag
# ===========================================================================

class TestSvigCanonicalFlag:
    @pytest.mark.parametrize("variant", [TERT_C228T, TERT_C250T])
    def test_set_for_both_tert_promoter_hotspots(self, variant):
        df = flags(row(**variant))
        assert df["ANNOTATION_FLAGS"].iloc[0] == "svig_canonical"

    def test_matching_is_genomic_not_by_hgvs(self):
        # c.-124C>T depends on which TERT transcript is annotated; the position does not
        df = flags(row(HGVSc="c.-124C>T", **TERT_C228T))
        assert df["ANNOTATION_FLAGS"].iloc[0] == "svig_canonical"

    def test_not_set_for_a_different_alt_at_the_same_position(self):
        v = dict(TERT_C228T, Tumor_Seq_Allele2="C")
        assert flags(row(**v))["ANNOTATION_FLAGS"].iloc[0] == ""

    def test_not_set_for_a_neighbouring_position(self):
        v = dict(TERT_C228T, Start_Position=1295114)
        assert flags(row(**v))["ANNOTATION_FLAGS"].iloc[0] == ""

    def test_set_whether_or_not_the_chromosome_carries_the_chr_prefix(self):
        v = dict(TERT_C228T, Chromosome="5")
        assert flags(row(**v))["ANNOTATION_FLAGS"].iloc[0] == "svig_canonical"

    def test_set_even_with_no_spliceai_evidence(self):
        df = flags(row(scores=(".", ".", ".", "."), **TERT_C228T))
        assert df["ANNOTATION_FLAGS"].iloc[0] == "svig_canonical"


# ===========================================================================
# both reasons, and the column contract
# ===========================================================================

class TestAnnotationFlagsColumn:
    def test_two_reasons_are_semicolon_separated(self):
        df = flags(row(impact="MODIFIER", scores=("0.95", ".", ".", "."), **TERT_C228T))
        assert df["ANNOTATION_FLAGS"].iloc[0] == "splice_predicted(0.95);svig_canonical"

    def test_the_column_always_exists_even_when_nothing_is_flagged(self):
        df = flags(row(impact="MODERATE"))
        assert "ANNOTATION_FLAGS" in df.columns
        assert df["ANNOTATION_FLAGS"].iloc[0] == ""

    def test_variant_classification_is_not_rewritten(self):
        df = flags(row(impact="LOW", scores=("0.99", ".", ".", "."),
                       Variant_Classification="Silent"), canonical=set())
        assert df["Variant_Classification"].iloc[0] == "Silent"

    def test_unflagged_rows_are_empty_not_the_word_none(self):
        # an object column holding None writes the literal "None" in some paths
        df = flags(row(impact="MODERATE"), canonical=set())
        buf = io.StringIO()
        df[["SPLICEAI_MAX_DS", "ANNOTATION_FLAGS"]].to_csv(buf, sep="\t", index=False)
        body = buf.getvalue().splitlines()[1]
        assert body == "\t"
        assert "None" not in body

    def test_the_count_of_flagged_rows_is_reported(self, capsys):
        flags(row(impact="LOW", scores=("0.95", ".", ".", ".")),
              row(**TERT_C228T),
              row(impact="MODERATE"), canonical=None)
        assert "2 row(s)" in capsys.readouterr().out

    def test_rows_are_judged_independently(self):
        df = flags(row(impact="LOW", scores=("0.95", ".", ".", ".")),
                   row(impact="LOW", scores=("0.01", ".", ".", ".")),
                   row(**TERT_C228T))
        assert list(df["ANNOTATION_FLAGS"]) == ["splice_predicted(0.95)", "",
                                                "svig_canonical"]

    def test_an_empty_frame_gains_both_columns(self):
        df = post.add_annotation_flags(pd.DataFrame(columns=["IMPACT", *DS]), canonical=set())
        assert {"SPLICEAI_MAX_DS", "ANNOTATION_FLAGS"} <= set(df.columns)


# ===========================================================================
# wiring
# ===========================================================================

class TestWiring:
    def test_both_columns_are_described_in_the_shipped_config(self):
        import yaml
        root = os.path.join(os.path.dirname(__file__), "..", "..")
        cfg = yaml.safe_load(open(os.path.join(root, "Config.yaml")))["fields"]
        for col in ("SPLICEAI_MAX_DS", "ANNOTATION_FLAGS"):
            assert col in cfg, f"{col} undescribed: it would be dropped before Tier 3"
            assert cfg[col]["tier"] == 3, f"{col} must reach the sheet scientists read"

    def test_the_threshold_matches_the_o6_threshold(self):
        assert post.SPLICEAI_FLAG_THRESHOLD == 0.2

    def test_the_flag_is_scoped_to_the_impacts_a_filter_removes(self):
        assert set(post.SPLICEAI_FLAG_IMPACTS) == {"LOW", "MODIFIER"}


# ===========================================================================
# the dtypes real data arrives in
# ===========================================================================

class TestRealDtypes:
    """A real MAF gives the SpliceAI columns as float64 with NaN, not as '.' strings.

    pandas types a numeric column on read, so the empty cells arrive as NaN and the
    scores as floats. The string forms above are what a hand-built frame gives; these
    are what the pipeline actually passes in.
    """

    def test_float_scores_are_read(self):
        df = post.add_annotation_flags(pd.DataFrame([
            dict(Chromosome="chr7", Start_Position=100, Reference_Allele="C",
                 Tumor_Seq_Allele2="T", IMPACT="LOW",
                 **dict(zip(DS, [0.95, float("nan"), 0.11, float("nan")])))]),
            canonical=set())
        assert df["SPLICEAI_MAX_DS"].iloc[0] == pytest.approx(0.95)
        assert df["ANNOTATION_FLAGS"].iloc[0] == "splice_predicted(0.95)"

    def test_all_nan_is_no_score_and_no_flag(self):
        df = post.add_annotation_flags(pd.DataFrame([
            dict(Chromosome="chr7", Start_Position=100, Reference_Allele="C",
                 Tumor_Seq_Allele2="T", IMPACT="LOW",
                 **dict(zip(DS, [float("nan")] * 4)))]), canonical=set())
        assert pd.isna(df["SPLICEAI_MAX_DS"].iloc[0])
        assert df["ANNOTATION_FLAGS"].iloc[0] == ""

    def test_a_float_zero_is_a_score_not_an_absence(self):
        # 0.00 is SpliceAI saying "no effect", which is not the same as no prediction
        df = post.add_annotation_flags(pd.DataFrame([
            dict(Chromosome="chr7", Start_Position=100, Reference_Allele="C",
                 Tumor_Seq_Allele2="T", IMPACT="LOW",
                 **dict(zip(DS, [0.0, float("nan"), float("nan"), float("nan")])))]),
            canonical=set())
        assert df["SPLICEAI_MAX_DS"].iloc[0] == 0.0
        assert not pd.isna(df["SPLICEAI_MAX_DS"].iloc[0])

    def test_an_integer_position_matches_the_canonical_list(self):
        # Start_Position is int64 in a real MAF, and the list holds it as text
        df = post.add_annotation_flags(pd.DataFrame([
            dict(Chromosome="chr5", Start_Position=1295113, Reference_Allele="G",
                 Tumor_Seq_Allele2="A", IMPACT="MODIFIER",
                 **dict(zip(DS, [float("nan")] * 4)))]))
        assert df["ANNOTATION_FLAGS"].iloc[0] == "svig_canonical"


# ===========================================================================
# against a real MAF committed to the repository
# ===========================================================================

REAL_MAF = os.path.join(os.path.dirname(__file__), "..", "verification1", "output",
                        "FINAL_Table", "verification1.oncoKB.maf")


@pytest.mark.skipif(not os.path.exists(REAL_MAF), reason="verification1 output not present")
class TestAgainstARealMaf:
    """The unit tests above build a frame of nine columns. A real MAF has over a
    thousand, with the dtypes pandas infers and the values VEP and OncoKB wrote, and it
    passes through every other transformation before reaching the flags."""

    @staticmethod
    def _real():
        return pd.read_csv(REAL_MAF, sep="\t", comment="#", low_memory=False)

    def test_every_column_the_flags_depend_on_is_present(self):
        cols = self._real().columns
        for c in ("IMPACT", "Chromosome", "Start_Position", "Reference_Allele",
                  "Tumor_Seq_Allele2", *DS):
            assert c in cols, f"{c} absent from a real MAF: the flags would be silent"

    def test_the_columns_are_added_to_a_real_frame(self):
        df = post.add_annotation_flags(self._real())
        assert df["SPLICEAI_MAX_DS"].dtype == "float64"
        assert (df["ANNOTATION_FLAGS"] == "").all(), \
            "no row of this sample qualifies; a flag here would be a false positive"

    def test_scores_are_recovered_from_the_real_cells(self):
        df = post.add_annotation_flags(self._real())
        assert df["SPLICEAI_MAX_DS"].notna().sum() > 0, \
            "real SpliceAI cells parsed to nothing"

    def test_a_real_row_made_eligible_is_flagged(self):
        # the sample contains a Silent/LOW row for which SpliceAI returned nothing;
        # give it a score and the row must be kept, with everything else untouched
        df = self._real()
        low = df.index[df["IMPACT"].astype(str).str.upper() == "LOW"]
        assert len(low), "no LOW row in this sample"
        i = low[0]
        before = df.at[i, "Variant_Classification"]
        df.at[i, "SpliceAI_pred_DS_AG"] = 0.95
        out = post.add_annotation_flags(df)
        assert out.at[i, "ANNOTATION_FLAGS"] == "splice_predicted(0.95)"
        assert out.at[i, "Variant_Classification"] == before
        assert (out["ANNOTATION_FLAGS"] != "").sum() == 1

    def test_a_real_row_moved_to_a_canonical_position_is_flagged(self):
        df = self._real()
        i = df.index[0]
        df.at[i, "Chromosome"] = "chr5"
        df.at[i, "Start_Position"] = 1295113
        df.at[i, "Reference_Allele"] = "G"
        df.at[i, "Tumor_Seq_Allele2"] = "A"
        out = post.add_annotation_flags(df)
        assert "svig_canonical" in out.at[i, "ANNOTATION_FLAGS"]

    def test_the_flags_survive_the_whole_per_sample_transform(self):
        # add_annotation_flags runs inside _transform_sample, before drop_columns
        root = os.path.join(os.path.dirname(__file__), "..", "..")
        field_config, field_patterns = post.load_field_config(
            os.path.join(root, "Config.yaml"))
        out, _, _ = post._transform_sample(self._real(), field_config)
        assert "SPLICEAI_MAX_DS" in out.columns, "dropped by a later step"
        assert "ANNOTATION_FLAGS" in out.columns, "dropped by a later step"
        t2, t3 = post.split_by_tier(out, field_config, field_patterns)
        for c in ("SPLICEAI_MAX_DS", "ANNOTATION_FLAGS"):
            assert c in t3.columns, f"{c} does not reach the Tier 3 sheet"
            assert c in t2.columns
