"""
Unit tests for scripts/post_analysis.py

Covers:
  - load_field_config()           — YAML loading and validation
  - resolve_field_meta()          — exact + pattern field lookup
  - standardize_format_column()   — FORMAT column renaming
  - normalize_format_fields()     — GT/AD/VAF expansion
  - add_vaf_column()              — VAF calculation from allelic depth
  - clean_column_values()         — value sanitisation
  - calculate_end_position()      — end-coordinate logic
  - drop_columns()                — column removal
"""
import pandas as pd
import pytest

from conftest import post


# ===========================================================================
# load_field_config
# ===========================================================================

class TestLoadFieldConfig:
    def test_returns_fields_and_patterns(self, minimal_yaml):
        fields, patterns = post.load_field_config(minimal_yaml)
        assert isinstance(fields, dict)
        assert isinstance(patterns, list)

    def test_fields_contain_expected_keys(self, minimal_yaml):
        fields, _ = post.load_field_config(minimal_yaml)
        assert "SYMBOL" in fields
        assert "ONCOKB_ONCOGENIC" in fields

    def test_patterns_contain_wildcard_entries(self, minimal_yaml):
        _, patterns = post.load_field_config(minimal_yaml)
        assert len(patterns) >= 1
        assert any("ONCOKB_TX_*" in str(p) for p in patterns)

    def test_missing_fields_key_raises_value_error(self, tmp_path):
        p = tmp_path / "bad.yaml"
        p.write_text("other_key:\n  foo: bar\n")
        with pytest.raises(ValueError, match="fields"):
            post.load_field_config(str(p))

    def test_nonexistent_file_raises(self, tmp_path):
        with pytest.raises((FileNotFoundError, OSError)):
            post.load_field_config(str(tmp_path / "missing.yaml"))

    def test_field_metadata_structure(self, minimal_yaml):
        fields, _ = post.load_field_config(minimal_yaml)
        symbol_meta = fields["SYMBOL"]
        assert "description" in symbol_meta
        assert "source" in symbol_meta
        assert "version" in symbol_meta


# ===========================================================================
# resolve_field_meta
# ===========================================================================

class TestResolveFieldMeta:
    def setup_method(self):
        self.fields = {
            "SYMBOL": {"description": "Gene symbol", "source": "VEP", "tier": 3},
            "ONCOKB_ONCOGENIC": {"description": "Oncogenicity", "source": "OncoKB", "tier": 3},
        }
        self.patterns = [
            {"pattern": "ONCOKB_TX_*", "description": "Treatment", "source": "OncoKB", "tier": 2},
            {"pattern": "ONCOKB_DIAG_*", "description": "Diagnostic", "source": "OncoKB", "tier": 2},
        ]

    def test_exact_match(self):
        meta = post.resolve_field_meta("SYMBOL", self.fields, self.patterns)
        assert meta["source"] == "VEP"

    def test_exact_match_takes_priority_over_pattern(self):
        meta = post.resolve_field_meta("ONCOKB_ONCOGENIC", self.fields, self.patterns)
        assert meta["description"] == "Oncogenicity"

    def test_pattern_match_oncokb_tx(self):
        meta = post.resolve_field_meta("ONCOKB_TX_0_level", self.fields, self.patterns)
        assert meta["source"] == "OncoKB"
        assert meta["tier"] == 2

    def test_pattern_match_oncokb_diag(self):
        meta = post.resolve_field_meta("ONCOKB_DIAG_1_tumorType", self.fields, self.patterns)
        assert meta["description"] == "Diagnostic"

    def test_no_match_returns_empty_dict(self):
        meta = post.resolve_field_meta("COMPLETELY_UNKNOWN_FIELD", self.fields, self.patterns)
        assert meta == {}


# ===========================================================================
# standardize_format_column
# ===========================================================================

class TestStandardizeFormatColumn:
    def test_renames_column_after_format(self):
        df = pd.DataFrame({"FORMAT": ["GT:AD"], "SAMPLE": ["0/1:50,25"]})
        result = post.standardize_format_column(df)
        assert "FORMAT_DATA" in result.columns
        assert "SAMPLE" not in result.columns

    def test_no_format_column_unchanged(self):
        df = pd.DataFrame({"SYMBOL": ["BRAF"], "HGVSp": ["p.V600E"]})
        result = post.standardize_format_column(df)
        assert list(result.columns) == ["SYMBOL", "HGVSp"]

    def test_format_last_column_no_crash(self):
        df = pd.DataFrame({"SYMBOL": ["BRAF"], "FORMAT": ["GT"]})
        result = post.standardize_format_column(df)
        assert "FORMAT" in result.columns


# ===========================================================================
# normalize_format_fields
# ===========================================================================

class TestNormalizeFormatFields:
    def test_expands_gt_ad_vaf(self):
        df = pd.DataFrame({
            "FORMAT":      ["GT:AD:VAF"],
            "FORMAT_DATA": ["0/1:50,25:0.33"],
        })
        result = post.normalize_format_fields(df)
        assert "GT" in result.columns
        assert "AD" in result.columns
        assert "VAF" in result.columns
        assert result["GT"].iloc[0] == "0/1"
        assert result["AD"].iloc[0] == "50,25"

    def test_missing_format_columns_returns_unchanged(self):
        df = pd.DataFrame({"SYMBOL": ["NRAS"]})
        result = post.normalize_format_fields(df)
        assert "GT" not in result.columns

    def test_fewer_values_than_keys_handled(self):
        df = pd.DataFrame({
            "FORMAT":      ["GT:AD:VAF"],
            "FORMAT_DATA": ["0/1:50,25"],   # VAF missing
        })
        result = post.normalize_format_fields(df)
        assert "GT" in result.columns
        assert "AD" in result.columns
        # VAF key exists but value is None
        assert result["VAF"].iloc[0] is None

    def test_existing_column_not_overwritten(self):
        df = pd.DataFrame({
            "FORMAT":      ["GT:AD"],
            "FORMAT_DATA": ["0/1:10,5"],
            "GT":          ["existing"],
        })
        result = post.normalize_format_fields(df)
        # GT was already present — should not be overwritten
        assert result["GT"].iloc[0] == "existing"


# ===========================================================================
# add_vaf_column
# ===========================================================================

class TestAddVafColumn:
    def test_normal_calculation(self):
        df = pd.DataFrame({"AD": ["50,25"]})
        result = post.add_vaf_column(df)
        expected = 25 / 75
        assert abs(result["VAF"].iloc[0] - expected) < 1e-9

    def test_zero_alt_reads(self):
        df = pd.DataFrame({"AD": ["100,0"]})
        result = post.add_vaf_column(df)
        assert result["VAF"].iloc[0] == 0.0

    def test_zero_total_returns_zero(self):
        df = pd.DataFrame({"AD": ["0,0"]})
        result = post.add_vaf_column(df)
        assert result["VAF"].iloc[0] == 0

    def test_invalid_format_returns_none(self):
        df = pd.DataFrame({"AD": ["garbage"]})
        result = post.add_vaf_column(df)
        assert result["VAF"].iloc[0] is None

    def test_missing_ad_column_unchanged(self):
        df = pd.DataFrame({"SYMBOL": ["BRAF"]})
        result = post.add_vaf_column(df)
        assert "VAF" not in result.columns

    @pytest.mark.parametrize("ad,expected", [
        ("62,30",   30 / 92),    # IDH1 R132H from verification1
        ("55,23",   23 / 78),    # NRAS Q61R
        ("100,100", 0.5),
    ])
    def test_real_vaf_values(self, ad, expected):
        df = pd.DataFrame({"AD": [ad]})
        result = post.add_vaf_column(df)
        assert abs(result["VAF"].iloc[0] - expected) < 1e-9

    def test_multiple_rows(self):
        df = pd.DataFrame({"AD": ["100,50", "80,20", "60,60"]})
        result = post.add_vaf_column(df)
        assert len(result["VAF"]) == 3
        assert abs(result["VAF"].iloc[0] - 50 / 150) < 1e-9
        assert abs(result["VAF"].iloc[1] - 20 / 100) < 1e-9
        assert abs(result["VAF"].iloc[2] - 0.5) < 1e-9


# ===========================================================================
# clean_column_values
# ===========================================================================

class TestCleanColumnValues:
    def test_empty_brackets_replaced(self):
        df = pd.DataFrame({"FIELD": ["[]"]})
        result = post.clean_column_values(df)
        assert result["FIELD"].iloc[0] == ""

    def test_url_encoded_equals_replaced(self):
        # df.replace with regex=False matches the whole cell value, not substrings.
        # Only a cell that IS exactly "%3D" gets replaced with "=".
        df = pd.DataFrame({"FIELD": ["%3D"]})
        result = post.clean_column_values(df)
        assert result["FIELD"].iloc[0] == "="

    def test_hgvsc_transcript_prefix_stripped(self):
        df = pd.DataFrame({"HGVSc": ["NM_005228.5:c.2573T>G"]})
        result = post.clean_column_values(df)
        assert result["HGVSc"].iloc[0] == "c.2573T>G"

    def test_hgvsc_without_colon_unchanged(self):
        df = pd.DataFrame({"HGVSc": ["c.2573T>G"]})
        result = post.clean_column_values(df)
        assert result["HGVSc"].iloc[0] == "c.2573T>G"

    def test_normal_values_unchanged(self):
        df = pd.DataFrame({"SYMBOL": ["BRAF"], "ONCOGENIC": ["Oncogenic"]})
        result = post.clean_column_values(df)
        assert result["SYMBOL"].iloc[0] == "BRAF"
        assert result["ONCOGENIC"].iloc[0] == "Oncogenic"


# ===========================================================================
# calculate_end_position
# ===========================================================================

class TestCalculateEndPosition:
    def _df(self, **kwargs) -> pd.DataFrame:
        defaults = {
            "Start_Position": 1000,
            "ID": ".",
            "REF": "A",
            "ALT": "G",
            "DUPSVLEN": "",
        }
        defaults.update(kwargs)
        return pd.DataFrame([defaults])

    def test_snv_end_equals_start(self):
        df = self._df(REF="A", ALT="G")    # len(REF)==len(ALT)==1
        result = post.calculate_end_position(df)
        assert result["End_Position"].iloc[0] == 1000

    def test_deletion_end_past_start(self):
        # REF=AG (2), ALT=A (1) → max(2,1)=2 → end = 1000 + 2 - 1 = 1001
        df = self._df(REF="AG", ALT="A")
        result = post.calculate_end_position(df)
        assert result["End_Position"].iloc[0] == 1001

    def test_manta_bnd_returns_none(self):
        df = self._df(ID="MantaBND:1:1000:2000")
        result = post.calculate_end_position(df)
        assert result["End_Position"].iloc[0] is None

    def test_manta_ins_with_symbolic_alt_uses_dupsvlen(self):
        df = self._df(ID="MantaINS:1:1000:2000", ALT="<INS>", DUPSVLEN=10)
        result = post.calculate_end_position(df)
        assert result["End_Position"].iloc[0] == 1010

    def test_manta_ins_with_sequence_alt_uses_alt_length(self):
        df = self._df(ID="MantaINS:1:1000:2000", ALT="ACTG", DUPSVLEN="")
        result = post.calculate_end_position(df)
        assert result["End_Position"].iloc[0] == 1004   # 1000 + len("ACTG")

    def test_manta_dup_uses_dupsvlen(self):
        df = self._df(ID="MantaDUP:1:1000:1500", DUPSVLEN=500)
        result = post.calculate_end_position(df)
        assert result["End_Position"].iloc[0] == 1500

    def test_manta_dup_missing_svlen_returns_none(self):
        df = self._df(ID="MantaDUP:1:1000:1500", DUPSVLEN="")
        result = post.calculate_end_position(df)
        assert result["End_Position"].iloc[0] is None

    def test_end_position_column_created(self):
        df = self._df()
        result = post.calculate_end_position(df)
        assert "End_Position" in result.columns


# ===========================================================================
# drop_columns
# ===========================================================================

class TestDropColumns:
    def test_drops_existing_columns(self):
        df = pd.DataFrame({"A": [1], "B": [2], "C": [3]})
        result = post.drop_columns(df, ["A", "C"])
        assert list(result.columns) == ["B"]

    def test_ignores_nonexistent_columns(self):
        df = pd.DataFrame({"A": [1]})
        result = post.drop_columns(df, ["MISSING", "ALSO_MISSING"])
        assert list(result.columns) == ["A"]

    def test_empty_list_unchanged(self):
        df = pd.DataFrame({"A": [1], "B": [2]})
        result = post.drop_columns(df, [])
        assert list(result.columns) == ["A", "B"]

    def test_returns_dataframe(self):
        df = pd.DataFrame({"A": [1]})
        result = post.drop_columns(df, [])
        assert isinstance(result, pd.DataFrame)
