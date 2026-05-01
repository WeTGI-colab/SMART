#!/usr/bin/env python3
"""
test_multi_transcript.py — Unit + integration tests for multi-transcript output (v0.2.0).

Tests the core behaviour introduced in v0.2.0: when a variant overlaps more than
one transcript in the preferred whitelist, the pipeline produces one output row
per matching transcript, each with its own VEP annotation and OncoKB result.

Two levels of testing — no Docker, no reference files, no API token required:

  Unit tests (TestFindAllPreferredCsqMatches, TestFindFallbackCsq)
    Exercise find_all_preferred_csq_matches() and find_fallback_csq() directly
    using crafted CSQ strings with the real VEP field positions.

  Integration test (TestVcf2TableMultiTranscript)
    Writes a synthetic two-record VCF (simulating oncokb2.0.py output for an
    EGFR variant matched by two preferred transcripts), runs vcf2table.py's
    main() on it, and asserts each record produces a correctly populated row.

Usage:
    python3 tests/test_multi_transcript.py          # run all tests
    python3 tests/test_multi_transcript.py -v       # verbose output
"""

import sys
import os
import unittest
import tempfile
import argparse
import pandas as pd

# Make the scripts/ directory importable without installing anything
SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), "..", "scripts")
sys.path.insert(0, os.path.abspath(SCRIPTS_DIR))

import types as _types
import importlib.util as _ilu

# Stub heavy dependencies that are only available inside Docker
for _mod_name in ("cyvcf2", "requests"):
    if _mod_name not in sys.modules:
        sys.modules[_mod_name] = _types.ModuleType(_mod_name)

# Ensure cyvcf2.VCF / Writer exist as dummy classes so module-level imports succeed
_cyvcf2 = sys.modules["cyvcf2"]
for _cls in ("VCF", "Writer"):
    if not hasattr(_cyvcf2, _cls):
        setattr(_cyvcf2, _cls, type(_cls, (), {}))

_spec = _ilu.spec_from_file_location("oncokb2_0", os.path.join(SCRIPTS_DIR, "oncokb2.0.py"))
oncokb = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(oncokb)
find_all_preferred_csq_matches = oncokb.find_all_preferred_csq_matches
find_fallback_csq              = oncokb.find_fallback_csq

import vcf2table as v2t


# ---------------------------------------------------------------------------
# CSQ field positions matching the real VEP output used by this pipeline.
# Derived from the ##INFO=<ID=CSQ,...> header of a pipeline-annotated VCF.
# ---------------------------------------------------------------------------
N_CSQ_FIELDS   = 152   # total number of pipe-separated fields
IDX_SYMBOL     = 3
IDX_FEATURE    = 6
IDX_HGVSP      = 11
IDX_MANE       = 25
IDX_MANE_SEL   = 26    # MANE_SELECT   (contains NM_ id when applicable)
IDX_MANE_PLUS  = 27    # MANE_PLUS_CLINICAL


def make_csq(symbol="", feature="", hgvsp="", mane="", mane_select="", mane_plus="") -> str:
    """Build one CSQ annotation string with the correct number of fields."""
    fields = [""] * N_CSQ_FIELDS
    fields[IDX_SYMBOL]   = symbol
    fields[IDX_FEATURE]  = feature
    fields[IDX_HGVSP]    = hgvsp
    fields[IDX_MANE]     = mane
    fields[IDX_MANE_SEL] = mane_select
    fields[IDX_MANE_PLUS]= mane_plus
    return "|".join(fields)


class MockVariant:
    """Minimal stand-in for a cyvcf2 Variant — only used for log output."""
    CHROM = "chr7"
    POS   = 55191822
    REF   = "T"
    ALT   = ["G"]


# ===========================================================================
# Unit tests — find_all_preferred_csq_matches
# ===========================================================================

class TestFindAllPreferredCsqMatches(unittest.TestCase):

    def _run(self, csqs: str, preferred: set) -> list:
        return find_all_preferred_csq_matches(
            MockVariant(), csqs, preferred,
            IDX_SYMBOL, IDX_FEATURE, IDX_HGVSP,
            IDX_MANE_SEL, IDX_MANE_PLUS,
        )

    def test_two_preferred_transcripts_returns_two_matches(self):
        """Core feature: two preferred NM IDs → two entries, preserving order."""
        egfr_mane = make_csq(
            symbol="EGFR", feature="NM_005228.5",
            hgvsp="ENSP00000275493.2:p.Leu858Arg",
            mane="MANE_Select", mane_select="NM_005228.5",
        )
        egfr_alt = make_csq(
            symbol="EGFR", feature="NM_001346897.2",
            hgvsp="ENSP00000494659.1:p.Leu813Arg",
            mane_plus="NM_001346897.2",
        )
        csqs = f"{egfr_mane},{egfr_alt}"
        preferred = {"NM_005228", "NM_001346897"}

        matches = self._run(csqs, preferred)

        self.assertEqual(len(matches), 2)
        self.assertEqual(matches[0]["nm_id"], "NM_005228")
        self.assertEqual(matches[0]["gene"],  "EGFR")
        self.assertEqual(matches[0]["protein_change"], "L858R")
        self.assertEqual(matches[1]["nm_id"], "NM_001346897")
        self.assertEqual(matches[1]["protein_change"], "L813R")

    def test_one_preferred_transcript_returns_one_match(self):
        """Only one of two transcripts is preferred → single match."""
        egfr_mane = make_csq(
            symbol="EGFR", feature="NM_005228.5",
            hgvsp="ENSP:p.Leu858Arg", mane_select="NM_005228.5",
        )
        egfr_alt = make_csq(
            symbol="EGFR", feature="NM_001346897.2",
            hgvsp="ENSP:p.Leu813Arg", mane_plus="NM_001346897.2",
        )
        csqs = f"{egfr_mane},{egfr_alt}"

        matches = self._run(csqs, {"NM_005228"})

        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0]["nm_id"], "NM_005228")

    def test_no_preferred_match_returns_empty_list(self):
        """No preferred transcript in CSQ → empty list (fall through to Tier 2/3)."""
        egfr = make_csq(symbol="EGFR", feature="NM_005228.5",
                        hgvsp="ENSP:p.Leu858Arg", mane_select="NM_005228.5")

        matches = self._run(egfr, {"NM_999999"})

        self.assertEqual(matches, [])

    def test_empty_preferred_set_returns_empty_list(self):
        """No preferred transcript file provided → empty list."""
        egfr = make_csq(symbol="EGFR", mane_select="NM_005228.5")

        matches = self._run(egfr, set())

        self.assertEqual(matches, [])

    def test_duplicate_nm_ids_deduplicated(self):
        """Same NM ID appearing in both MANE_SELECT and MANE_PLUS → one entry."""
        ann = make_csq(
            symbol="GENE", mane_select="NM_000001.1", mane_plus="NM_000001.1",
        )
        matches = self._run(ann, {"NM_000001"})

        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0]["nm_id"], "NM_000001")

    def test_version_suffix_stripped(self):
        """NM_000077.5 in CSQ matches NM_000077 in preferred set."""
        cdkn2a = make_csq(
            symbol="CDKN2A", feature="NM_000077.5",
            hgvsp="ENSP:p.Ala68Val", mane_select="NM_000077.5",
        )
        matches = self._run(cdkn2a, {"NM_000077"})

        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0]["nm_id"], "NM_000077")

    def test_mane_plus_clinical_matched(self):
        """Transcript found via MANE_PLUS_CLINICAL field."""
        ann = make_csq(
            symbol="CDKN2A", feature="NM_058195.4",
            hgvsp="ENSP:p.Pro48Gly", mane_plus="NM_058195.4",
        )
        matches = self._run(ann, {"NM_058195"})

        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0]["nm_id"], "NM_058195")

    def test_cdkn2a_both_isoforms_in_preferred(self):
        """CDKN2A p16 (NM_000077) and p14ARF (NM_058195) → two rows."""
        p16 = make_csq(
            symbol="CDKN2A", feature="NM_000077.5",
            hgvsp="ENSP:p.Ala68Val", mane_select="NM_000077.5",
        )
        p14arf = make_csq(
            symbol="CDKN2A", feature="NM_058195.4",
            hgvsp="ENSP:p.Pro48Gly", mane_plus="NM_058195.4",
        )
        csqs = f"{p16},{p14arf}"
        preferred = {"NM_000077", "NM_058195"}

        matches = self._run(csqs, preferred)

        self.assertEqual(len(matches), 2)
        nms = {m["nm_id"] for m in matches}
        self.assertEqual(nms, {"NM_000077", "NM_058195"})


# ===========================================================================
# Unit tests — find_fallback_csq
# ===========================================================================

class TestFindFallbackCsq(unittest.TestCase):

    def _run(self, csqs: str):
        return find_fallback_csq(
            MockVariant(), csqs,
            IDX_SYMBOL, IDX_FEATURE, IDX_HGVSP,
            IDX_MANE_SEL, IDX_MANE_PLUS, IDX_MANE,
        )

    def test_picks_mane_select_over_first(self):
        first = make_csq(symbol="EGFR", feature="NM_999.1", mane="")
        mane  = make_csq(symbol="EGFR", feature="NM_005228.5", mane="MANE_Select")
        gene, protein, tx, reason = self._run(f"{first},{mane}")

        self.assertIn("MANE Select", reason)
        # Returns stripped NM ID (no version suffix) so vcf2table.py can match it
        self.assertEqual(tx, "NM_005228")

    def test_picks_mane_plus_when_no_mane_select(self):
        first = make_csq(symbol="EGFR", feature="NM_999.1")
        plus  = make_csq(symbol="EGFR", feature="NM_058195.4",
                         mane="MANE_Plus_Clinical")
        gene, protein, tx, reason = self._run(f"{first},{plus}")

        self.assertIn("MANE Plus Clinical", reason)

    def test_falls_back_to_first_when_no_mane(self):
        first  = make_csq(symbol="GENE_A", feature="NM_111.1")
        second = make_csq(symbol="GENE_B", feature="NM_222.2")
        gene, protein, tx, reason = self._run(f"{first},{second}")

        self.assertEqual(gene, "GENE_A")
        self.assertIn("Fallback", reason)

    def test_empty_csq_returns_none(self):
        gene, protein, tx, reason = self._run("")

        self.assertIsNone(gene)
        self.assertIsNone(protein)


# ===========================================================================
# Integration test — vcf2table.py with a synthetic two-record VCF
# ===========================================================================

# Minimal CSQ format for the synthetic test VCF — only the 6 fields we care
# about.  vcf2table.py uses the header to build the field list dynamically, so
# the data rows must match exactly.
_MINI_CSQ_FORMAT = "SYMBOL|Feature|HGVSp|MANE|MANE_SELECT|MANE_PLUS_CLINICAL"
_MINI_N = 6
_MINI_IDX = {
    "SYMBOL": 0, "Feature": 1, "HGVSp": 2,
    "MANE": 3, "MANE_SELECT": 4, "MANE_PLUS_CLINICAL": 5,
}


def _mini_csq(symbol="", feature="", hgvsp="", mane="", mane_select="", mane_plus="") -> str:
    return f"{symbol}|{feature}|{hgvsp}|{mane}|{mane_select}|{mane_plus}"


_VCF_HEADER = f"""\
##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP consequence. Format: {_MINI_CSQ_FORMAT}">
##INFO=<ID=ONCOKB_PREFERRED_TRANSCRIPT,Number=1,Type=String,Description="Preferred transcript NM ID">
##INFO=<ID=ONCOKB_QUERY_TYPE,Number=1,Type=String,Description="Query type">
##INFO=<ID=ONCOKB_ONCOGENIC,Number=1,Type=String,Description="Oncogenic classification">
##INFO=<ID=ONCOKB_SENS_LVL,Number=1,Type=String,Description="Sensitivity level">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""

# Two records for the same EGFR L858R variant position — one per preferred
# transcript, exactly as oncokb2.0.py would write them.
_EGFR_CSQ_BOTH = (
    _mini_csq("EGFR", "NM_005228.5",  "ENSP:p.Leu858Arg", "MANE_Select", "NM_005228.5",  "") + "," +
    _mini_csq("EGFR", "NM_001346897.2","ENSP:p.Leu813Arg", "",            "",              "NM_001346897.2")
)

_VCF_RECORDS = (
    # Record 1 — annotated with the MANE Select transcript (LEVEL_1)
    f"chr7\t55191822\t.\tT\tG\t.\tPASS\t"
    f"CSQ={_EGFR_CSQ_BOTH};"
    f"ONCOKB_PREFERRED_TRANSCRIPT=NM_005228;"
    f"ONCOKB_QUERY_TYPE=MUTATION;"
    f"ONCOKB_ONCOGENIC=Oncogenic;"
    f"ONCOKB_SENS_LVL=LEVEL_1\tGT\t0/1\n"
    # Record 2 — annotated with the alternative transcript (Unknown)
    f"chr7\t55191822\t.\tT\tG\t.\tPASS\t"
    f"CSQ={_EGFR_CSQ_BOTH};"
    f"ONCOKB_PREFERRED_TRANSCRIPT=NM_001346897;"
    f"ONCOKB_QUERY_TYPE=MUTATION;"
    f"ONCOKB_ONCOGENIC=Unknown;"
    f"ONCOKB_SENS_LVL=.\tGT\t0/1\n"
)


class TestVcf2TableMultiTranscript(unittest.TestCase):
    """
    End-to-end test of vcf2table.py with a synthetic two-record VCF.

    Verifies that:
    1. Each input record produces exactly one output row.
    2. The ONCOKB_PREFERRED_TRANSCRIPT tag correctly pins the CSQ selection —
       row 1 gets the MANE Select transcript, row 2 gets the alternative.
    3. The OncoKB annotations from the INFO field (Oncogenic, SENS_LVL) appear
       on the correct row.
    """

    def setUp(self):
        # Write synthetic VCF
        self.vcf_file = tempfile.NamedTemporaryFile(
            mode="w", suffix=".vcf", delete=False
        )
        self.vcf_file.write(_VCF_HEADER + _VCF_RECORDS)
        self.vcf_file.flush()

        # Transcript whitelist — both EGFR NM IDs
        self.tx_file = tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        )
        self.tx_file.write("NM_005228.5\nNM_001346897.2\n")
        self.tx_file.flush()

        self.out_file = tempfile.NamedTemporaryFile(
            suffix=".csv", delete=False
        )
        self.out_file.close()

        args = argparse.Namespace(
            vcf=self.vcf_file.name,
            transcripts=self.tx_file.name,
            output=self.out_file.name,
            debug=False,
        )
        df = v2t.main(args)
        self.df = df

    def tearDown(self):
        for f in (self.vcf_file, self.tx_file):
            try:
                os.unlink(f.name)
            except OSError:
                pass
        try:
            os.unlink(self.out_file.name)
        except OSError:
            pass

    def test_two_input_records_produce_two_output_rows(self):
        self.assertEqual(len(self.df), 2,
            f"Expected 2 rows (one per preferred transcript), got {len(self.df)}")

    def test_row1_uses_mane_select_transcript(self):
        row = self.df.iloc[0]
        self.assertEqual(row["NM_Transcript"], "NM_005228",
            "Row 1 should be pinned to the MANE Select transcript NM_005228")

    def test_row2_uses_alternative_transcript(self):
        row = self.df.iloc[1]
        self.assertEqual(row["NM_Transcript"], "NM_001346897",
            "Row 2 should be pinned to the alternative transcript NM_001346897")

    def test_row1_vep_fields_match_mane_select_csq(self):
        row = self.df.iloc[0]
        self.assertEqual(row.get("SYMBOL", ""), "EGFR")
        self.assertEqual(row.get("Feature", ""), "NM_005228.5")
        self.assertIn("Leu858Arg", row.get("HGVSp", ""))

    def test_row2_vep_fields_match_alternative_csq(self):
        row = self.df.iloc[1]
        self.assertEqual(row.get("SYMBOL", ""), "EGFR")
        self.assertEqual(row.get("Feature", ""), "NM_001346897.2")
        self.assertIn("Leu813Arg", row.get("HGVSp", ""))

    def test_row1_oncokb_annotation_is_oncogenic(self):
        row = self.df.iloc[0]
        self.assertEqual(row.get("ONCOKB_ONCOGENIC", ""), "Oncogenic")
        self.assertEqual(row.get("ONCOKB_SENS_LVL", ""), "LEVEL_1")

    def test_row2_oncokb_annotation_is_unknown(self):
        row = self.df.iloc[1]
        self.assertEqual(row.get("ONCOKB_ONCOGENIC", ""), "Unknown")
        # LEVEL_1 must NOT appear on the second row
        self.assertNotEqual(row.get("ONCOKB_SENS_LVL", ""), "LEVEL_1")

    def test_both_rows_same_genomic_position(self):
        """Both rows originate from the same variant — coordinates are identical."""
        self.assertEqual(self.df.iloc[0]["POS"], self.df.iloc[1]["POS"])
        self.assertEqual(self.df.iloc[0]["REF"], self.df.iloc[1]["REF"])
        self.assertEqual(self.df.iloc[0]["ALT"], self.df.iloc[1]["ALT"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
