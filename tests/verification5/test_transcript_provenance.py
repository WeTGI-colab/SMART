"""
Unit tests for scripts/transcript_provenance.py.

A user supplies `--transcripts-file` and expects the output to describe that transcript. The
report answers, per column, whether it does -- and the answer is not uniformly yes:

  resolved     the value carries the transcript or the step matched on it
  positional   the source is keyed on position and has no transcript to offer
  unresolved   the source is transcript-specific, and the output does not resolve it to ours

The third category is what this exists for. dbNSFP returns one value per transcript and the
plugin runs without `transcript_match=1`, so the cell keeps all of them and marks none. On
verification1 that is 29 columns and 233 instances where the reported transcript's value
differs from the cell's first element -- in both directions: for TP53 the first element is a
different transcript's score, and for MLH1 the first element is `.` while the reported
transcript has a real value at index 6.

The detection is structural, not a hardcoded column list: a column parallel to
`Ensembl_transcriptid` on every row that has such a list is one value per transcript. A new
dbNSFP column is therefore covered the day it appears.
"""
import importlib.util
import os
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "..", "scripts")
_spec = importlib.util.spec_from_file_location(
    "transcript_provenance", os.path.join(SCRIPTS, "transcript_provenance.py"))
tp = importlib.util.module_from_spec(_spec)
sys.modules["transcript_provenance"] = tp
_spec.loader.exec_module(tp)

REAL_MAF = os.path.join(os.path.dirname(__file__), "..", "verification1",
                        "output", "FINAL_Table", "verification1.oncoKB.maf")


def write_maf(tmp_path, header, rows):
    p = tmp_path / "t.maf"
    with open(p, "w") as fh:
        fh.write("#SMART_VERSION 2.1.0\n")
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return str(p)


class TestFindingThePerTranscriptColumns:
    HDR = ["Hugo_Symbol", "Feature", "Ensembl_transcriptid", "AlphaMissense_score",
           "DOMAINS", "CADD_PHRED"]

    def test_a_parallel_column_is_found(self, tmp_path):
        p = write_maf(tmp_path, self.HDR, [
            ["TP53", "ENST2", "ENST1&ENST2&ENST3", "0.1&0.2&0.3", "a&b", "22.5"]])
        h, r = tp.read_table(p)
        idx = {c: i for i, c in enumerate(h)}
        cols, n = tp.per_transcript_columns(h, r, idx)
        assert "AlphaMissense_score" in cols
        assert n == 1

    def test_a_list_of_a_different_length_is_not_parallel(self, tmp_path):
        """DOMAINS and PUBMED are lists by nature, not one entry per transcript."""
        p = write_maf(tmp_path, self.HDR, [
            ["TP53", "ENST2", "ENST1&ENST2&ENST3", "0.1&0.2&0.3", "a&b", "22.5"]])
        h, r = tp.read_table(p)
        idx = {c: i for i, c in enumerate(h)}
        cols, _ = tp.per_transcript_columns(h, r, idx)
        assert "DOMAINS" not in cols

    def test_a_coincidence_on_one_row_is_not_enough(self, tmp_path):
        """Several unrelated columns match the element count by chance on a single row.

        Requiring it on every row with a transcript list is what separates a real
        per-transcript column from an accident.
        """
        p = write_maf(tmp_path, self.HDR, [
            ["TP53", "ENST2", "ENST1&ENST2", "0.1&0.2", "x&y", "22.5"],
            ["KRAS", "ENSTB", "ENSTA&ENSTB&ENSTC", "0.3&0.4&0.5", "z", "19.0"]])
        h, r = tp.read_table(p)
        idx = {c: i for i, c in enumerate(h)}
        cols, n = tp.per_transcript_columns(h, r, idx)
        assert n == 2
        assert "AlphaMissense_score" in cols
        assert "DOMAINS" not in cols, "a one-row coincidence was taken as per-transcript"

    def test_a_single_valued_column_is_not_parallel(self, tmp_path):
        p = write_maf(tmp_path, self.HDR, [
            ["TP53", "ENST2", "ENST1&ENST2", "0.1&0.2", "d", "22.5"]])
        h, r = tp.read_table(p)
        idx = {c: i for i, c in enumerate(h)}
        cols, _ = tp.per_transcript_columns(h, r, idx)
        assert "CADD_PHRED" not in cols


class TestClassification:
    def test_values_carrying_the_transcript_are_resolved(self):
        for col in ("HGVSc", "HGVSp", "SpliceAI_pred_SYMBOL"):
            assert tp.classify(col, []) == "resolved"

    def test_loeuf_is_resolved_because_the_step_matches_on_transcript(self):
        """entrypoint.sh runs the plugin with match_by=transcript."""
        assert tp.classify("LOEUF", []) == "resolved"

    def test_revel_is_unresolved_because_it_runs_with_no_match(self):
        assert tp.classify("REVEL", []) == "unresolved"
        assert "no_match=1" in tp.UNRESOLVED_REASONS["REVEL"]

    @pytest.mark.parametrize("col", ["ClinVar_CLNSIG", "CIViC_CSQ_SYMBOL", "COSMIC_CNT",
                                     "CancerHotspots_HOTSPOT", "CADD_PHRED", "gnomAD_AF"])
    def test_position_keyed_sources_are_positional(self, col):
        """Not a gap: these sources have no transcript to honour."""
        assert tp.classify(col, []) == "positional"

    def test_a_per_transcript_column_is_unresolved(self):
        assert tp.classify("AlphaMissense_score", ["AlphaMissense_score"]) == "unresolved"

    def test_the_same_column_would_be_other_if_it_were_not_parallel(self):
        """The classification follows the data, not the column's name."""
        assert tp.classify("AlphaMissense_score", []) == "other"


class TestDivergences:
    HDR = ["Hugo_Symbol", "Feature", "Ensembl_transcriptid", "AlphaMissense_score"]

    def _div(self, tmp_path, rows):
        p = write_maf(tmp_path, self.HDR, rows)
        h, r = tp.read_table(p)
        idx = {c: i for i, c in enumerate(h)}
        per, _ = tp.per_transcript_columns(h, r, idx)
        return tp.divergences(h, r, idx, per)

    def test_our_value_differing_from_the_first_is_reported(self, tmp_path):
        """The TP53 case: the first element belongs to another transcript."""
        d = self._div(tmp_path, [["TP53", "ENST3", "ENST1&ENST2&ENST3", "0.1&0.2&0.9"]])
        assert len(d) == 1
        assert d[0]["index"] == 2 and d[0]["ours"] == "0.9" and d[0]["first"] == "0.1"

    def test_a_first_element_with_no_value_is_reported(self, tmp_path):
        """The MLH1 case, and the worse direction: the reader sees '.' and concludes
        there is no prediction, while the reported transcript has one."""
        d = self._div(tmp_path, [["MLH1", "ENST3", "ENST1&ENST2&ENST3", ".&.&0.9976"]])
        assert len(d) == 1
        assert d[0]["ours"] == "0.9976" and d[0]["first"] == "."

    def test_agreement_is_not_reported(self, tmp_path):
        d = self._div(tmp_path, [["TP53", "ENST1", "ENST1&ENST2", "0.5&0.5"]])
        assert d == []

    def test_our_transcript_being_absent_is_reported(self, tmp_path):
        """Worth knowing on its own: the source has no value for the transcript we report."""
        d = self._div(tmp_path, [["TP53", "ENST9", "ENST1&ENST2", "0.1&0.2"]])
        assert len(d) == 1 and d[0]["note"] == "reported transcript absent"

    def test_the_id_columns_are_not_reported_as_divergences(self, tmp_path):
        """Ensembl_transcriptid differing from its own first element is not a finding."""
        hdr = self.HDR + ["Ensembl_proteinid", "genename"]
        p = write_maf(tmp_path, hdr,
                      [["TP53", "ENST2", "ENST1&ENST2", "0.1&0.2", "ENSP1&ENSP2", "TP53&TP53"]])
        h, r = tp.read_table(p)
        idx = {c: i for i, c in enumerate(h)}
        per, _ = tp.per_transcript_columns(h, r, idx)
        cols = {d["column"] for d in tp.divergences(h, r, idx, per)}
        assert "Ensembl_transcriptid" not in cols
        assert "Ensembl_proteinid" not in cols


class TestTheRealOutput:
    @pytest.mark.skipif(not os.path.exists(REAL_MAF),
                        reason="verification1 has not been run")
    def test_the_predictors_are_no_longer_per_transcript_lists(self):
        """The fix for #29: the pipeline now runs dbNSFP with transcript_match=1.

        This test previously asserted the opposite -- that AlphaMissense and the rest *were*
        parallel to Ensembl_transcriptid -- because that was the state it was written to
        describe. Turning the flag on broke it, which is the right way round: a test that
        encodes a defect should fail when the defect is fixed, not quietly keep passing.
        """
        h, r = tp.read_table(REAL_MAF)
        idx = {c: i for i, c in enumerate(h)}
        per, _ = tp.per_transcript_columns(h, r, idx)
        for col in ("AlphaMissense_score", "MetaRNN_score", "PROVEAN_score", "VEST4_score",
                    "ESM1b_score", "DEOGEN2_score"):
            assert col not in per, (
                f"{col} is a per-transcript list again; is the pipeline running with "
                f"--dbnsfp-transcript-match 0, or did the flag stop being passed?")

    @pytest.mark.skipif(not os.path.exists(REAL_MAF),
                        reason="verification1 has not been run")
    def test_mutation_taster_is_still_a_list_and_that_is_documented(self):
        """The one limit the flag cannot fix, asserted so it is not forgotten.

        dbNSFP keys MutationTaster on a different transcript set, so the plugin's own
        documentation says transcript_match has no effect on those fields. They stay
        unresolved, and the report should keep saying so.
        """
        h, r = tp.read_table(REAL_MAF)
        idx = {c: i for i, c in enumerate(h)}
        # Not per_transcript_columns: with transcript_match=1 the anchor it needs
        # (Ensembl_transcriptid as a list) is itself resolved to one transcript, so that
        # detector correctly finds nothing. A column still holding a list in that state is
        # what unresolved_remnants is for.
        rem = tp.unresolved_remnants(h, r, idx)
        assert any(c.startswith("MutationTaster") for c in rem), (
            f"MutationTaster is no longer a list -- if dbNSFP changed how it keys those "
            f"fields, the caveat in entrypoint.sh and issue #29 should be removed. "
            f"Remnants found: {sorted(rem)}")
        assert all(c.startswith("MutationTaster") for c in rem), (
            f"something other than MutationTaster is still an unresolved list: "
            f"{sorted(c for c in rem if not c.startswith('MutationTaster'))}")

    @pytest.mark.skipif(not os.path.exists(REAL_MAF),
                        reason="verification1 has not been run")
    def test_dbnsfp_matched_the_transcript_the_row_reports(self):
        """The positive proof the flag was added for.

        With transcript_match=1, `Ensembl_transcriptid` holds the transcript dbNSFP matched
        and `Feature` holds the one VEP selected from the user's list. Their agreeing is what
        makes "the predictors describe the transcript you asked for" a checkable statement
        rather than a promise.
        """
        h, r = tp.read_table(REAL_MAF)
        idx = {c: i for i, c in enumerate(h)}
        ok, bad = tp.transcript_list_agrees_with_feature(h, r, idx)
        assert ok > 0, "no row had a resolved dbNSFP transcript to compare"
        assert bad == 0, (
            f"{bad} row(s) carry dbNSFP predictors for a transcript other than the one the "
            f"row reports")

    @pytest.mark.skipif(not os.path.exists(REAL_MAF),
                        reason="verification1 has not been run")
    def test_hgvsc_names_the_reported_transcript(self):
        """The resolved case, asserted on real output rather than taken on trust."""
        h, r = tp.read_table(REAL_MAF)
        idx = {c: i for i, c in enumerate(h)}
        checked = 0
        for row in r:
            feature = tp.cell(row, idx, "Feature")
            hgvsc = tp.cell(row, idx, "HGVSc")
            if not feature or not hgvsc or ":" not in hgvsc:
                continue
            assert hgvsc.split(":")[0].split(".")[0] == feature.split(".")[0], (
                f"HGVSc {hgvsc} does not describe the reported transcript {feature}")
            checked += 1
        assert checked >= 20, f"only {checked} rows had both fields; too few to mean anything"
