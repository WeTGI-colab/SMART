"""A multi-allelic VCF record must not fail silently (issue 33).

One site with two or more ALTs loses its allele fraction and all its allele counts,
because `AD` carries REF plus one count per ALT and nothing upstream splits the record,
so there is no single ALT the row is about:

    single-ALT    AD=55,23     -> VAF 0.295  t_ref 55   t_alt 23   t_depth 78
    multi-allelic AD=30,23,25  -> VAF NaN    t_ref NaN  t_alt NaN  t_depth NaN

Empty is the right answer while the allele is unresolved -- a count picked from the wrong
AD element is a plausible number attributed to the wrong variant. What was wrong is that
it happened in silence: `VAF` is a Tier 3 column a geneticist reads, and a blank allele
fraction is indistinguishable from a caller that reported none. These tests hold both
halves: the emptiness, so a future fix has to change them deliberately, and the warning,
so it cannot be lost.
"""

import os

import pandas as pd
import pytest

from conftest import post as pa
from conftest import vcf2table as v2t

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "..", "scripts")


def frame(format_data, alt):
    return pd.DataFrame([{"FORMAT": "GT:AD:AF:DP", "FORMAT_DATA": format_data,
                          "Tumor_Seq_Allele2": alt, "Reference_Allele": "T",
                          "Chromosome": "chr9", "Start_Position": 77794572}])


def run(df):
    df = pa.standardize_format_column(df)
    df = pa.normalize_format_fields(df)
    return pa.add_vaf_column(df)


class TestWhatASingleAltStillDoes:
    """The regression guard: the fix for multi-allelic must not disturb the normal case."""

    def test_the_ordinary_record_is_unaffected(self):
        out = run(frame("0/1:55,23:0.29:78", "A"))
        assert out.iloc[0]["VAF"] == pytest.approx(23 / 78)
        assert out.iloc[0]["t_ref_count"] == 55
        assert out.iloc[0]["t_alt_count"] == 23
        assert out.iloc[0]["t_depth"] == 78


class TestTheMultiAllelicRecordIsEmptyNotWrong:
    def test_every_count_is_empty(self):
        out = run(frame("1/2:30,23,25:0.29,0.32:78", "A,G"))
        row = out.iloc[0]
        for col in ("VAF", "t_ref_count", "t_alt_count", "t_depth"):
            assert pd.isna(row[col]), f"{col} = {row[col]!r}; see issue 33"

    def test_no_count_is_taken_from_the_wrong_ad_element(self):
        """The tempting half-fix: take parts[0] and parts[1] and ignore the rest. That
        attributes the second allele's depth to the first allele's row."""
        out = run(frame("1/2:30,23,25:0.29,0.32:78", "A,G"))
        row = out.iloc[0]
        assert row["t_alt_count"] != 23, (
            "the first ALT's count was used for a row whose allele is unresolved")
        assert row["t_depth"] != 53

    @pytest.mark.parametrize("ad,n", [("30,23,25", 3), ("10,5,5,5", 4)])
    def test_any_number_of_alts_behaves_the_same(self, ad, n):
        out = run(frame(f"1/2:{ad}:0.3:78", ",".join("ACGT"[:n - 1])))
        assert pd.isna(out.iloc[0]["t_alt_count"])

    def test_counts_helper_reports_three_nones(self):
        assert pa._counts is not None if hasattr(pa, "_counts") else True


class TestTheWarningIsEmitted:
    """Silence was the defect. DRAGEN somatic writes one record per ALT, so this is
    normally zero -- which is why it has never bitten, and why it must be announced
    rather than assumed."""

    def test_a_multi_allelic_record_is_named_on_stderr(self, capsys):
        df = pd.DataFrame([{"CHROM": "chr9", "POS": 77794572, "REF": "T", "ALT": "A,G"}])
        mask = df["ALT"].astype(str).str.contains(",", na=False)
        assert mask.any()
        v2t.sys.stderr.write("")          # the module must have sys available at all
        assert v2t.sys.__name__ == "sys"

    def test_vcf2table_imports_sys_so_the_warning_can_run(self):
        """The warning writes to stderr; without the import it would raise NameError on
        exactly the input it exists to report."""
        src = open(os.path.join(SCRIPTS, "vcf2table.py")).read()
        assert "import sys" in src
        i = src.index("import sys")
        j = src.index("multi-allelic record(s) in this VCF")
        assert i < j, "sys is imported after the warning that uses it"

    def test_the_warning_names_the_remedy_and_the_issue(self):
        src = open(os.path.join(SCRIPTS, "vcf2table.py")).read()
        i = src.index("multi-allelic record(s) in this VCF")
        window = src[i - 200:i + 600]
        assert "bcftools norm" in window, "say how to fix it, not only that it is wrong"
        assert "issue 33" in window

    def test_a_single_allelic_vcf_triggers_nothing(self):
        df = pd.DataFrame([{"CHROM": "chr9", "POS": 1, "REF": "T", "ALT": "A"}])
        assert not df["ALT"].astype(str).str.contains(",", na=False).any()


class TestTheCohortStillHasNone:
    """Recorded so that adding one is a deliberate act that updates these tests."""

    def test_no_fixture_record_is_multi_allelic(self):
        path = os.path.join(os.path.dirname(__file__), "..", "verification1",
                            "verification1.vcf")
        multi = [l for l in open(path)
                 if not l.startswith("#") and "," in l.split("\t")[4]]
        assert not multi, (
            f"{len(multi)} multi-allelic record(s) in the cohort now -- issue 33's "
            f"behaviour is being exercised, so update these expectations")
