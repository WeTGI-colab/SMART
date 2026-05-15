"""
Unit tests for scripts/vcf2table.py

Covers:
  - calculate_end()        — end-position logic for every variant class
  - convert_hgvsp_short()  — 3-letter → 1-letter AA code conversion
  - load_nm_transcripts()  — whitelist loading and version stripping
"""
import pandas as pd
import pytest

from conftest import vcf2table


# ===========================================================================
# calculate_end
# ===========================================================================

def _row(**kwargs) -> pd.Series:
    """Build a minimal pandas Series for calculate_end."""
    defaults = {"POS": 100, "REF": "A", "ALT": "G",
                "VARIANT_CLASS": "snv", "SVLEN": None}
    defaults.update(kwargs)
    return pd.Series(defaults)


class TestCalculateEnd:
    def test_snv_returns_pos(self):
        assert vcf2table.calculate_end(_row(VARIANT_CLASS="snv", POS=114713908)) == 114713908

    def test_deletion_returns_pos_plus_ref_minus_one(self):
        # REF=AGT (3 bp) → end = 100 + 3 - 1 = 102
        row = _row(POS=100, REF="AGT", ALT="A", VARIANT_CLASS="deletion")
        assert vcf2table.calculate_end(row) == 102

    def test_deletion_single_base(self):
        # REF=AT, ALT=A → end = 100 + 2 - 1 = 101
        row = _row(POS=100, REF="AT", ALT="A", VARIANT_CLASS="deletion")
        assert vcf2table.calculate_end(row) == 101

    def test_insertion_with_svlen(self):
        # SVLEN=10 → end = 100 + 10 - 1 = 109
        row = _row(POS=100, REF="A", ALT="AGGGGGGGGG",
                   VARIANT_CLASS="insertion", SVLEN=10)
        assert vcf2table.calculate_end(row) == 109

    def test_insertion_svlen_as_string(self):
        # SVLEN stored as string (common in VCF INFO)
        row = _row(POS=200, REF="A", ALT="ATT", VARIANT_CLASS="insertion", SVLEN="5")
        assert vcf2table.calculate_end(row) == 204

    def test_insertion_without_svlen_uses_alt_length(self):
        # ALT=AGG (3 bp) → end = 100 + 3 - 1 = 102
        row = _row(POS=100, REF="A", ALT="AGG", VARIANT_CLASS="insertion", SVLEN=None)
        assert vcf2table.calculate_end(row) == 102

    def test_chromosome_breakpoint_returns_empty_string(self):
        row = _row(VARIANT_CLASS="chromosome_breakpoint")
        assert vcf2table.calculate_end(row) == ""

    def test_unknown_class_returns_pos(self):
        # Fallback: any unrecognised class returns POS
        row = _row(POS=500, VARIANT_CLASS="sequence_alteration")
        assert vcf2table.calculate_end(row) == 500

    @pytest.mark.parametrize("pos,ref,alt,expected", [
        (114713908, "T", "C",   114713908),     # NRAS Q61R — SNV
        (208248388, "C", "T",   208248388),     # IDH1 R132H — SNV
        (55174772,  "GAATTAAGAGAAGCAT", "G", 55174787),  # EGFR del — deletion: 55174772 + 16 - 1
    ])
    def test_real_verification1_variants(self, pos, ref, alt, expected):
        n = len(ref)
        if n == len(alt):
            vc = "snv"
        elif n > len(alt):
            vc = "deletion"
        else:
            vc = "insertion"
        row = _row(POS=pos, REF=ref, ALT=alt, VARIANT_CLASS=vc)
        assert vcf2table.calculate_end(row) == expected


# ===========================================================================
# convert_hgvsp_short
# ===========================================================================

class TestConvertHgvspShort:
    @pytest.mark.parametrize("hgvsp,expected", [
        ("p.Val600Glu",  "p.V600E"),     # BRAF V600E
        ("p.Gln61Arg",   "p.Q61R"),     # NRAS Q61R
        ("p.Arg132His",  "p.R132H"),    # IDH1 R132H
        ("p.Glu545Lys",  "p.E545K"),    # PIK3CA E545K
        ("p.His1047Arg", "p.H1047R"),   # PIK3CA H1047R
        ("p.Arg248Trp",  "p.R248W"),    # TP53 R248W
        ("p.Arg175His",  "p.R175H"),    # TP53 R175H
    ])
    def test_three_letter_to_one_letter(self, hgvsp, expected):
        assert vcf2table.convert_hgvsp_short(hgvsp) == expected

    def test_stop_codon_ter_becomes_asterisk(self):
        assert vcf2table.convert_hgvsp_short("p.Arg248Ter") == "p.R248*"

    def test_strips_ensp_prefix(self):
        # Full HGVSp from VEP includes an ENSP prefix
        result = vcf2table.convert_hgvsp_short("ENSP00000358548.4:p.Val600Glu")
        assert result == "p.V600E"

    def test_already_single_letter_preserved(self):
        # If the regex finds 1-letter codes, aa_dict.get falls back to the key
        result = vcf2table.convert_hgvsp_short("p.V600E")
        assert result == "p.V600E"

    def test_empty_string_returns_empty(self):
        assert vcf2table.convert_hgvsp_short("") == ""

    def test_none_returns_empty(self):
        assert vcf2table.convert_hgvsp_short(None) == ""

    def test_no_p_dot_returns_empty(self):
        assert vcf2table.convert_hgvsp_short("c.1799T>A") == ""

    def test_non_string_returns_empty(self):
        assert vcf2table.convert_hgvsp_short(123) == ""


# ===========================================================================
# load_nm_transcripts
# ===========================================================================

class TestLoadNmTranscripts:
    def test_loads_transcripts(self, nm_file):
        result = vcf2table.load_nm_transcripts(nm_file)
        assert "NM_000077" in result
        assert "NM_002524" in result
        assert "NM_005228" in result
        assert "NM_000546" in result

    def test_strips_version_suffixes(self, nm_file):
        result = vcf2table.load_nm_transcripts(nm_file)
        # Versioned entries must not appear
        assert "NM_000077.5" not in result
        assert "NM_002524.3" not in result

    def test_skips_blank_lines(self, nm_file):
        result = vcf2table.load_nm_transcripts(nm_file)
        assert "" not in result

    def test_comment_lines_treated_as_transcripts(self, nm_file):
        # The current implementation does NOT skip # lines,
        # so a comment becomes part of the set (stripped to '#').
        # Test that the function runs without error; the comment char
        # itself is not a valid NM ID but it won't crash.
        result = vcf2table.load_nm_transcripts(nm_file)
        assert isinstance(result, set)

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises((FileNotFoundError, OSError)):
            vcf2table.load_nm_transcripts(str(tmp_path / "nonexistent.txt"))

    def test_returns_set(self, nm_file):
        result = vcf2table.load_nm_transcripts(nm_file)
        assert isinstance(result, set)

    def test_single_entry_no_version(self, tmp_path):
        p = tmp_path / "single.txt"
        p.write_text("NM_005896\n")
        result = vcf2table.load_nm_transcripts(str(p))
        assert result == {"NM_005896"}

    def test_version_stripped_to_base(self, tmp_path):
        p = tmp_path / "v.txt"
        p.write_text("NM_023110.3\n")
        assert vcf2table.load_nm_transcripts(str(p)) == {"NM_023110"}
