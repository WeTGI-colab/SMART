"""SMART's own OncoKB client must send the consequence, as MafAnnotator does.

`oncokb2.0.py` sent `byProteinChange` the gene and the protein string and nothing else.
OncoKB then answers only about that literal string, and for a variant whose
interpretation depends on its consequence the answer is wrong in a way that looks like
an answer. Measured against the live API for TP53 `T125=`:

    alteration=T125=                                   -> hotspot False, Unknown, no Px
    alteration=T125=&consequence=splice_region_variant -> hotspot True,  Likely Oncogenic, LEVEL_Px1
    alteration=T125=&consequence=synonymous_variant    -> hotspot False, Unknown, no Px

The cost in the delivered output was `ONCOKB_HOTSPOT` False where OncoKB says True,
`ONCOKB_PROG_LVL` empty where OncoKB says LEVEL_Px1, and a variant summary reading "is
not annotated by OncoKB" -- all three written by this client, in the same row where
MafAnnotator's own columns had it right, because MafAnnotator sends the consequence.

No network: these tests assert the request built and the strings chosen.
"""

import os
import sys
from unittest.mock import MagicMock, patch

import pytest

from conftest import oncokb


class TestWhichConsequenceTermIsSent:
    """OncoKB does not understand a compound consequence, so one term must be chosen.

    `splice_region_variant%2Bsynonymous_variant` -- the form AnnotatorCore builds when it
    maps a Variant_Classification -- was measured returning Unknown. VEP orders the terms
    of one annotation by severity, most severe first, so the leading term is the one that
    carries the interpretation.
    """

    @pytest.mark.parametrize("raw,expected", [
        ("splice_region_variant&synonymous_variant", "splice_region_variant"),
        ("missense_variant", "missense_variant"),
        ("start_lost&feature_elongation", "start_lost"),
        ("stop_lost&feature_truncation", "stop_lost"),
        ("splice_donor_variant&intron_variant", "splice_donor_variant"),
        # a comma can appear when a caller has already joined annotations
        ("frameshift_variant,missense_variant", "frameshift_variant"),
        ("  upstream_gene_variant  ", "upstream_gene_variant"),
    ])
    def test_the_leading_term_is_taken(self, raw, expected):
        assert oncokb.most_severe_consequence(raw) == expected

    @pytest.mark.parametrize("raw", ["", None, "&", "   "])
    def test_nothing_usable_yields_nothing(self, raw):
        """An empty consequence must be omitted, not sent blank: the client's previous
        behaviour (send nothing) has to remain reachable."""
        assert oncokb.most_severe_consequence(raw) == ""

    def test_the_synonymous_half_is_never_what_gets_sent(self):
        """The specific trap. Sending the second term of TP53 T125T's consequence returns
        Unknown, which is what the client used to report."""
        assert oncokb.most_severe_consequence(
            "splice_region_variant&synonymous_variant") != "synonymous_variant"


class TestTheRequestCarriesIt:
    def _call(self, consequence, resp_json=None):
        resp = MagicMock()
        resp.status_code = 200
        resp.json.return_value = resp_json or {"oncogenic": "Likely Oncogenic"}
        with patch("requests.get", return_value=resp) as mock_get:
            oncokb.query_oncokb_mutation(
                gene="TP53", alteration="T125=", tumor_type="", token="t",
                cache={}, consequence=consequence)
        return mock_get.call_args

    def test_the_consequence_is_a_query_parameter(self):
        args = self._call("splice_region_variant")
        params = args.kwargs.get("params") or args[1]["params"]
        assert params["consequence"] == "splice_region_variant"
        assert params["alteration"] == "T125="
        assert params["hugoSymbol"] == "TP53"

    def test_an_empty_consequence_is_omitted_not_sent_blank(self):
        params = (lambda a: a.kwargs.get("params") or a[1]["params"])(self._call(""))
        assert "consequence" not in params

    def test_it_defaults_to_omitted_so_other_callers_are_unaffected(self):
        resp = MagicMock()
        resp.status_code = 200
        resp.json.return_value = {"oncogenic": "Oncogenic"}
        with patch("requests.get", return_value=resp) as mock_get:
            oncokb.query_oncokb_mutation(gene="BRAF", alteration="V600E",
                                         tumor_type="", token="t", cache={})
        params = mock_get.call_args.kwargs.get("params") or mock_get.call_args[1]["params"]
        assert "consequence" not in params


class TestTheCacheDistinguishesThem:
    """The same protein string with a different consequence is a different question."""

    def test_two_consequences_do_not_share_a_cache_entry(self):
        a = oncokb.mutation_cache_key("TP53", "T125=", "", "splice_region_variant")
        b = oncokb.mutation_cache_key("TP53", "T125=", "", "synonymous_variant")
        assert a != b

    def test_the_same_question_shares_one(self):
        a = oncokb.mutation_cache_key("TP53", "T125=", "", "splice_region_variant")
        b = oncokb.mutation_cache_key("TP53", "T125=", "", "splice_region_variant")
        assert a == b

    def test_a_cached_answer_is_not_reused_for_a_different_consequence(self):
        cache = {oncokb.mutation_cache_key("TP53", "T125=", "", "synonymous_variant"):
                 {"oncogenic": "Unknown"}}
        resp = MagicMock()
        resp.status_code = 200
        resp.json.return_value = {"oncogenic": "Likely Oncogenic"}
        with patch("requests.get", return_value=resp) as mock_get:
            out = oncokb.query_oncokb_mutation(
                gene="TP53", alteration="T125=", tumor_type="", token="t",
                cache=cache, consequence="splice_region_variant")
        mock_get.assert_called_once()
        assert out["oncogenic"] == "Likely Oncogenic"


class TestTheFallbackPathFindsTheRightAnnotation:
    """The fallback returns a gene and protein change but not the consequence, so it has
    to be located again -- by Feature, or it can describe a different transcript."""

    # Consequence | Feature, two fields
    CSQ = "missense_variant|NM_000001.1,splice_region_variant&synonymous_variant|NM_000546.6"

    def test_it_matches_the_selected_transcript_not_the_first(self):
        got = oncokb.fallback_consequence(self.CSQ, "NM_000546.6", 0, 1)
        assert got == "splice_region_variant"

    def test_the_version_is_ignored_when_matching(self):
        assert oncokb.fallback_consequence(self.CSQ, "NM_000546", 0, 1) == "splice_region_variant"

    def test_with_no_transcript_it_takes_veps_own_first_pick(self):
        assert oncokb.fallback_consequence(self.CSQ, "", 0, 1) == "missense_variant"

    def test_a_transcript_that_matches_nothing_sends_nothing(self):
        """Better no consequence than another transcript's: a wrong consequence is a
        wrong question, and OncoKB answers it."""
        assert oncokb.fallback_consequence(self.CSQ, "NM_999999", 0, 1) == ""

    def test_a_vep_line_without_the_field_sends_nothing(self):
        assert oncokb.fallback_consequence(self.CSQ, "NM_000546.6", -1, 1) == ""

    def test_an_empty_csq_string_sends_nothing(self):
        assert oncokb.fallback_consequence("", "NM_000546.6", 0, 1) == ""
