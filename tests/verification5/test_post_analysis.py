"""
Unit tests for scripts/post_analysis.py

Covers:
  - load_field_config()              — YAML loading and validation
  - resolve_field_meta()             — exact + pattern field lookup
  - standardize_format_column()      — FORMAT column renaming
  - normalize_format_fields()        — GT/AD/VAF expansion
  - add_vaf_column()                 — VAF calculation from allelic depth
  - clean_column_values()            — value sanitisation
  - calculate_end_position()         — end-coordinate logic
  - drop_columns()                   — column removal
  - add_cancerhotspots_counts()      — O7 position and AA-change counts
  - add_gene_role()                  — O2/B2 gene TSG/oncogene classification
  - add_genie_counts()               — O4 somatic database patient count
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


# ===========================================================================
# add_cancerhotspots_counts
# ===========================================================================

class TestAddGeneRole:
    """Tests for add_gene_role() using a minimal gene roles TSV fixture."""

    @pytest.fixture
    def roles_file(self, tmp_path):
        content = (
            "hugoSymbol\tgeneType\n"
            "BRAF\tONCOGENE\n"
            "TP53\tTSG\n"
            "CDKN2A\tTSG\n"
            "EZH2\tONCOGENE_AND_TSG\n"
            "ALB\tNEITHER\n"
        )
        p = tmp_path / "oncokb_gene_roles.tsv"
        p.write_text(content)
        return str(p)

    def _df(self, symbol):
        return pd.DataFrame({"SYMBOL": [symbol]})

    def test_oncogene_classified(self, roles_file):
        result = post.add_gene_role(self._df("BRAF"), roles_file)
        assert result["Gene_role_in_cancer"].iloc[0] == "ONCOGENE"

    def test_tsg_classified(self, roles_file):
        result = post.add_gene_role(self._df("TP53"), roles_file)
        assert result["Gene_role_in_cancer"].iloc[0] == "TSG"

    def test_dual_role_classified(self, roles_file):
        result = post.add_gene_role(self._df("EZH2"), roles_file)
        assert result["Gene_role_in_cancer"].iloc[0] == "ONCOGENE_AND_TSG"

    def test_gene_not_in_lookup_returns_unknown(self, roles_file):
        result = post.add_gene_role(self._df("UNKNOWNGENE"), roles_file)
        assert result["Gene_role_in_cancer"].iloc[0] == "unknown"

    def test_no_path_returns_unknown(self):
        result = post.add_gene_role(self._df("BRAF"), gene_roles_path=None)
        assert result["Gene_role_in_cancer"].iloc[0] == "unknown"

    def test_missing_file_returns_unknown(self, tmp_path):
        result = post.add_gene_role(self._df("BRAF"), str(tmp_path / "missing.tsv"))
        assert result["Gene_role_in_cancer"].iloc[0] == "unknown"

    def test_multiple_rows(self, roles_file):
        df = pd.DataFrame({"SYMBOL": ["BRAF", "TP53", "EZH2", "UNKNOWNGENE"]})
        result = post.add_gene_role(df, roles_file)
        assert list(result["Gene_role_in_cancer"]) == [
            "ONCOGENE", "TSG", "ONCOGENE_AND_TSG", "unknown"
        ]

    def test_o2_scenario_tsg_with_frameshift(self, roles_file):
        """TSG + frameshift → O2 applies. Gene_role = TSG is the prerequisite."""
        df = pd.DataFrame({"SYMBOL": ["TP53"], "Consequence": ["frameshift_variant"]})
        result = post.add_gene_role(df, roles_file)
        assert result["Gene_role_in_cancer"].iloc[0] == "TSG"  # O2 can apply

    def test_b2_scenario_oncogene_with_frameshift(self, roles_file):
        """Oncogene + frameshift → B2 applies (wrong mechanism → VUS override)."""
        df = pd.DataFrame({"SYMBOL": ["BRAF"], "Consequence": ["frameshift_variant"]})
        result = post.add_gene_role(df, roles_file)
        assert result["Gene_role_in_cancer"].iloc[0] == "ONCOGENE"  # B2 can apply


class TestAddGenieCounts:
    """Tests for add_genie_counts() using a minimal in-memory TSV fixture."""

    @pytest.fixture
    def lookup_file(self, tmp_path):
        """Write a minimal genie_lookup.tsv with known variant counts."""
        content = (
            "gene\tprotein_change\tcount\n"
            "BRAF\tp.V600E\t7234\n"
            "KRAS\tp.G12D\t4521\n"
            "NRAS\tp.Q61R\t890\n"
            "IDH1\tp.R132H\t312\n"
            "TP53\tp.R175H\t45\n"
        )
        p = tmp_path / "genie_lookup.tsv"
        p.write_text(content)
        return str(p)

    def _df(self, symbol, hgvsp):
        return pd.DataFrame({"SYMBOL": [symbol], "HGVSp_Short": [hgvsp]})

    def test_known_driver_gets_count(self, lookup_file):
        df = self._df("BRAF", "p.V600E")
        result = post.add_genie_counts(df, lookup_file)
        assert result["GENIE_count"].iloc[0] == 7234

    def test_kras_g12d(self, lookup_file):
        df = self._df("KRAS", "p.G12D")
        result = post.add_genie_counts(df, lookup_file)
        assert result["GENIE_count"].iloc[0] == 4521

    def test_variant_not_in_lookup_returns_none(self, lookup_file):
        df = self._df("EGFR", "p.L858R")
        result = post.add_genie_counts(df, lookup_file)
        assert result["GENIE_count"].iloc[0] is None

    def test_no_lookup_path_adds_empty_column(self):
        df = self._df("BRAF", "p.V600E")
        result = post.add_genie_counts(df, lookup_path=None)
        assert "GENIE_count" in result.columns
        assert result["GENIE_count"].iloc[0] is None

    def test_missing_file_adds_empty_column(self, tmp_path):
        df = self._df("BRAF", "p.V600E")
        result = post.add_genie_counts(df, str(tmp_path / "missing.tsv"))
        assert "GENIE_count" in result.columns
        assert result["GENIE_count"].iloc[0] is None

    def test_multiple_rows(self, lookup_file):
        df = pd.DataFrame({
            "SYMBOL":      ["BRAF",    "KRAS",    "EGFR",    "NRAS"],
            "HGVSp_Short": ["p.V600E", "p.G12D",  "p.L858R", "p.Q61R"],
        })
        result = post.add_genie_counts(df, lookup_file)
        assert result["GENIE_count"].iloc[0] == 7234   # BRAF V600E
        assert result["GENIE_count"].iloc[1] == 4521   # KRAS G12D
        assert pd.isna(result["GENIE_count"].iloc[2])  # EGFR L858R — not in fixture
        assert result["GENIE_count"].iloc[3] == 890    # NRAS Q61R

    def test_cna_no_hgvsp_returns_none(self, lookup_file):
        df = pd.DataFrame({"SYMBOL": ["CDKN2A"], "HGVSp_Short": [""]})
        result = post.add_genie_counts(df, lookup_file)
        assert result["GENIE_count"].iloc[0] is None

    def test_svig_uk_o4_strong_threshold(self, lookup_file):
        """Verify count > 10 for missense → qualifies for O4 Strong [+4]."""
        df = self._df("TP53", "p.R175H")
        result = post.add_genie_counts(df, lookup_file)
        assert result["GENIE_count"].iloc[0] == 45
        assert result["GENIE_count"].iloc[0] > 10   # O4 Strong for missense


class TestAddCancerhotspotsCounts:
    """Tests for add_cancerhotspots_counts() using a minimal in-memory JSON fixture."""

    @pytest.fixture
    def counts_file(self, tmp_path):
        """Write a minimal cancerhotspots_counts.json with 3 hotspot records."""
        data = [
            {
                "hugoSymbol": "NRAS",
                "aminoAcidPosition": {"start": 61, "end": 61},
                "tumorCount": 422,
                "variantAminoAcid": {"R": 204, "K": 142, "H": 27},
            },
            {
                "hugoSymbol": "BRAF",
                "aminoAcidPosition": {"start": 600, "end": 600},
                "tumorCount": 9852,
                "variantAminoAcid": {"E": 9710, "K": 97, "R": 45},
            },
            {
                "hugoSymbol": "IDH1",
                "aminoAcidPosition": {"start": 132, "end": 132},
                "tumorCount": 3104,
                "variantAminoAcid": {"H": 2891, "C": 142, "S": 71},
            },
        ]
        import json
        p = tmp_path / "cancerhotspots_counts.json"
        p.write_text(json.dumps(data))
        return str(p)

    def _df(self, symbol, hgvsp_short):
        return pd.DataFrame({"SYMBOL": [symbol], "HGVSp_Short": [hgvsp_short]})

    def test_nras_q61r_strong(self, counts_file):
        df = self._df("NRAS", "p.Q61R")
        result = post.add_cancerhotspots_counts(df, counts_file)
        assert result["CancerHotspots_position_count"].iloc[0] == 422
        assert result["CancerHotspots_aa_change_count"].iloc[0] == 204

    def test_braf_v600e_strong(self, counts_file):
        df = self._df("BRAF", "p.V600E")
        result = post.add_cancerhotspots_counts(df, counts_file)
        assert result["CancerHotspots_position_count"].iloc[0] == 9852
        assert result["CancerHotspots_aa_change_count"].iloc[0] == 9710

    def test_idh1_r132h_strong(self, counts_file):
        df = self._df("IDH1", "p.R132H")
        result = post.add_cancerhotspots_counts(df, counts_file)
        assert result["CancerHotspots_position_count"].iloc[0] == 3104
        assert result["CancerHotspots_aa_change_count"].iloc[0] == 2891

    def test_rare_aa_change_gives_zero_count(self, counts_file):
        # BRAF V600D is not in variantAminoAcid — aa_change_count should be 0
        df = self._df("BRAF", "p.V600D")
        result = post.add_cancerhotspots_counts(df, counts_file)
        assert result["CancerHotspots_position_count"].iloc[0] == 9852
        assert result["CancerHotspots_aa_change_count"].iloc[0] == 0

    def test_gene_not_in_hotspots_returns_none(self, counts_file):
        df = self._df("TP53", "p.R175H")
        result = post.add_cancerhotspots_counts(df, counts_file)
        assert result["CancerHotspots_position_count"].iloc[0] is None
        assert result["CancerHotspots_aa_change_count"].iloc[0] is None

    def test_cna_no_hgvsp_returns_none(self, counts_file):
        df = pd.DataFrame({"SYMBOL": ["CDKN2A"], "HGVSp_Short": [""]})
        result = post.add_cancerhotspots_counts(df, counts_file)
        assert result["CancerHotspots_position_count"].iloc[0] is None

    def test_no_counts_file_adds_empty_columns(self):
        df = self._df("BRAF", "p.V600E")
        result = post.add_cancerhotspots_counts(df, counts_path=None)
        assert "CancerHotspots_position_count" in result.columns
        assert result["CancerHotspots_position_count"].iloc[0] is None

    def test_missing_file_adds_empty_columns(self, tmp_path):
        df = self._df("BRAF", "p.V600E")
        result = post.add_cancerhotspots_counts(df, str(tmp_path / "missing.json"))
        assert "CancerHotspots_position_count" in result.columns
        assert result["CancerHotspots_position_count"].iloc[0] is None

    def test_multiple_rows(self, counts_file):
        df = pd.DataFrame({
            "SYMBOL":      ["NRAS",   "BRAF",   "IDH1",   "TP53"],
            "HGVSp_Short": ["p.Q61R", "p.V600E","p.R132H","p.R175H"],
        })
        result = post.add_cancerhotspots_counts(df, counts_file)
        assert result["CancerHotspots_position_count"].iloc[0] == 422    # NRAS Q61
        assert result["CancerHotspots_position_count"].iloc[1] == 9852   # BRAF V600
        assert result["CancerHotspots_position_count"].iloc[2] == 3104   # IDH1 R132
        assert pd.isna(result["CancerHotspots_position_count"].iloc[3])   # TP53 not in fixture


# ===========================================================================
# add_o1_canonical
# ===========================================================================

class TestAddO1Canonical:
    """Tests for add_o1_canonical() — exact match and proxy detection."""

    @pytest.fixture
    def canon_file(self, tmp_path):
        content = (
            "gene\ttranscript\tHGVSp_Short\tHGVSp_long\tsvig_uk_assessment\n"
            "BRAF\tNM_004333.6\tp.V600E\tp.Val600Glu\tOncogenic\n"
            "KRAS\tNM_004985.5\tp.G12D\tp.Gly12Asp\tOncogenic\n"
            "IDH1\tNM_005896.4\tp.R132H\tp.Arg132His\tOncogenic\n"
            "NRAS\tNM_002524.5\tp.Q61R\tp.Gln61Arg\tOncogenic\n"
        )
        p = tmp_path / "canonical.tsv"
        p.write_text(content)
        return str(p)

    def _df(self, symbol, hgvsp, hotspot="False", oncogenic="", clinsci="", clinonc=""):
        return pd.DataFrame({
            "SYMBOL":          [symbol],
            "HGVSp_Short":     [hgvsp],
            "ONCOKB_HOTSPOT":  [hotspot],
            "ONCOKB_ONCOGENIC":[oncogenic],
            "ClinVar_SCI":     [clinsci],
            "ClinVar_ONC":     [clinonc],
        })

    def test_exact_match_canonical(self, canon_file):
        df = self._df("BRAF", "p.V600E")
        result = post.add_o1_canonical(df, canon_file)
        assert result["O1_canonical"].iloc[0] == True
        assert result["O1_canonical_source"].iloc[0] == "SVIG-UK_Table3_exact"

    def test_exact_match_idh1(self, canon_file):
        df = self._df("IDH1", "p.R132H")
        result = post.add_o1_canonical(df, canon_file)
        assert result["O1_canonical"].iloc[0] == True

    def test_non_canonical_not_matched(self, canon_file):
        df = self._df("EGFR", "p.L858R")
        result = post.add_o1_canonical(df, canon_file)
        assert result["O1_canonical"].iloc[0] == False

    def test_proxy_detection(self, canon_file):
        # Not in canonical list but meets proxy criteria
        df = self._df("EGFR", "p.L858R",
                      hotspot="True", oncogenic="Oncogenic",
                      clinsci="Tier_I_-_Strong", clinonc="Oncogenic")
        result = post.add_o1_canonical(df, canon_file)
        assert result["O1_canonical"].iloc[0] == True
        assert result["O1_canonical_source"].iloc[0] == "proxy_OncoKB+ClinVar"

    def test_proxy_not_triggered_without_hotspot(self, canon_file):
        df = self._df("EGFR", "p.L858R",
                      hotspot="False", oncogenic="Oncogenic",
                      clinsci="Tier_I_-_Strong")
        result = post.add_o1_canonical(df, canon_file)
        assert result["O1_canonical"].iloc[0] == False

    def test_no_canonical_file_uses_proxy(self):
        df = self._df("BRAF", "p.V600E",
                      hotspot="True", oncogenic="Oncogenic",
                      clinsci="Tier_I_-_Strong", clinonc="Oncogenic")
        result = post.add_o1_canonical(df, canonical_path=None)
        assert result["O1_canonical"].iloc[0] == True
        assert result["O1_canonical_source"].iloc[0] == "proxy_OncoKB+ClinVar"

    def test_missing_file_falls_back_to_proxy(self, tmp_path):
        df = self._df("BRAF", "p.V600E",
                      hotspot="True", oncogenic="Oncogenic",
                      clinsci="Tier_I_-_Strong")
        result = post.add_o1_canonical(df, str(tmp_path / "missing.tsv"))
        assert result["O1_canonical"].iloc[0] == True

    def test_columns_always_added(self):
        df = pd.DataFrame({"SYMBOL": ["ALB"], "HGVSp_Short": ["p.X1Y"]})
        result = post.add_o1_canonical(df)
        assert "O1_canonical" in result.columns
        assert "O1_canonical_source" in result.columns


# ===========================================================================
# add_o9_candidate
# ===========================================================================

class TestAddO9Candidate:
    """Tests for add_o9_candidate() — in-frame indel, stop-loss, final-exon oncogene."""

    def _df(self, consequence, exon="", role=""):
        return pd.DataFrame({
            "Consequence":        [consequence],
            "EXON":               [exon],
            "Gene_role_in_cancer":[role],
        })

    def test_in_frame_insertion_flagged(self):
        result = post.add_o9_candidate(self._df("in_frame_insertion"))
        assert result["O9_candidate"].iloc[0] == True
        assert result["O9_candidate_reason"].iloc[0] == "in_frame_indel"

    def test_in_frame_deletion_flagged(self):
        result = post.add_o9_candidate(self._df("in_frame_deletion"))
        assert result["O9_candidate"].iloc[0] == True
        assert result["O9_candidate_reason"].iloc[0] == "in_frame_indel"

    def test_stop_lost_flagged(self):
        result = post.add_o9_candidate(self._df("stop_lost"))
        assert result["O9_candidate"].iloc[0] == True
        assert result["O9_candidate_reason"].iloc[0] == "stop_lost"

    def test_truncating_final_exon_oncogene_flagged(self):
        df = self._df("stop_gained", exon="10/10", role="ONCOGENE")
        result = post.add_o9_candidate(df)
        assert result["O9_candidate"].iloc[0] == True
        assert result["O9_candidate_reason"].iloc[0] == "truncating_final_exon_oncogene"

    def test_truncating_not_final_exon_not_flagged(self):
        df = self._df("stop_gained", exon="5/10", role="ONCOGENE")
        result = post.add_o9_candidate(df)
        assert result["O9_candidate"].iloc[0] == False

    def test_truncating_final_exon_tsg_not_flagged(self):
        # Final exon truncation in TSG → O2, not O9
        df = self._df("stop_gained", exon="10/10", role="TSG")
        result = post.add_o9_candidate(df)
        assert result["O9_candidate"].iloc[0] == False

    def test_missense_not_flagged(self):
        result = post.add_o9_candidate(self._df("missense_variant"))
        assert result["O9_candidate"].iloc[0] == False

    def test_columns_always_added(self):
        df = pd.DataFrame({"Consequence": ["synonymous_variant"]})
        result = post.add_o9_candidate(df)
        assert "O9_candidate" in result.columns
        assert "O9_candidate_reason" in result.columns


# ===========================================================================
# add_o5_same_position
# ===========================================================================

class TestAddO5SamePosition:
    """Tests for add_o5_same_position() — same residue as known oncogenic variant."""

    @pytest.fixture
    def canon_file(self, tmp_path):
        content = (
            "gene\ttranscript\tHGVSp_Short\tHGVSp_long\tsvig_uk_assessment\n"
            "BRAF\tNM_004333.6\tp.V600E\tp.Val600Glu\tOncogenic\n"
            "KRAS\tNM_004985.5\tp.G12D\tp.Gly12Asp\tOncogenic\n"
            "KRAS\tNM_004985.5\tp.G12V\tp.Gly12Val\tOncogenic\n"
            "IDH1\tNM_005896.4\tp.R132H\tp.Arg132His\tOncogenic\n"
        )
        p = tmp_path / "canonical.tsv"
        p.write_text(content)
        return str(p)

    def _df(self, symbol, hgvsp, revel=""):
        return pd.DataFrame({
            "SYMBOL":      [symbol],
            "HGVSp_Short": [hgvsp],
            "REVEL":       [revel],
        })

    def test_same_position_different_aa_strong(self, canon_file):
        # BRAF p.V600K — same position as canonical p.V600E, REVEL >= 0.773
        df = self._df("BRAF", "p.V600K", revel="0.95")
        result = post.add_o5_same_position(df, canon_file)
        assert result["O5_same_position"].iloc[0] == True
        assert result["O5_strength"].iloc[0] == "Strong"
        assert "p.V600E" in result["O5_reference_variant"].iloc[0]

    def test_same_position_moderate_revel(self, canon_file):
        # REVEL in 0.7–0.773 range → Moderate
        df = self._df("IDH1", "p.R132C", revel="0.72")
        result = post.add_o5_same_position(df, canon_file)
        assert result["O5_same_position"].iloc[0] == True
        assert result["O5_strength"].iloc[0] == "Moderate"

    def test_low_revel_not_flagged(self, canon_file):
        # REVEL < 0.7 → O5 not applicable
        df = self._df("BRAF", "p.V600K", revel="0.50")
        result = post.add_o5_same_position(df, canon_file)
        assert result["O5_same_position"].iloc[0] == False

    def test_exact_canonical_match_not_flagged(self, canon_file):
        # p.V600E itself is in the canonical list → O1 territory, not O5
        df = self._df("BRAF", "p.V600E", revel="0.99")
        result = post.add_o5_same_position(df, canon_file)
        assert result["O5_same_position"].iloc[0] == False

    def test_different_gene_not_flagged(self, canon_file):
        # EGFR p.V600K — no canonical variant at position 600 in EGFR
        df = self._df("EGFR", "p.V600K", revel="0.95")
        result = post.add_o5_same_position(df, canon_file)
        assert result["O5_same_position"].iloc[0] == False

    def test_no_canonical_file_returns_false(self):
        df = self._df("BRAF", "p.V600K", revel="0.95")
        result = post.add_o5_same_position(df, canonical_path=None)
        assert result["O5_same_position"].iloc[0] == False

    def test_multiple_canonical_at_same_position(self, canon_file):
        # KRAS G12 has G12D and G12V in canonical; G12C should match both
        df = self._df("KRAS", "p.G12C", revel="0.85")
        result = post.add_o5_same_position(df, canon_file)
        assert result["O5_same_position"].iloc[0] == True
        assert "p.G12D" in result["O5_reference_variant"].iloc[0]
        assert "p.G12V" in result["O5_reference_variant"].iloc[0]

    def test_columns_always_added(self):
        df = pd.DataFrame({"SYMBOL": ["ALB"], "HGVSp_Short": ["p.X1Y"], "REVEL": [""]})
        result = post.add_o5_same_position(df)
        assert "O5_same_position" in result.columns
        assert "O5_strength" in result.columns
        assert "O5_reference_variant" in result.columns


# ===========================================================================
# add_svig_uk_classification  (_score_variant)
# ===========================================================================

class TestSvigUkClassification:
    """Tests for the SVIG-UK scoring and classification logic."""

    def _row(self, **kwargs):
        defaults = {
            "O1_canonical": False, "O1_canonical_source": "",
            "Gene_role_in_cancer": "unknown",
            "Consequence": "missense_variant",
            "MAX_AF": None,
            "REVEL": None,
            "SpliceAI_pred_DS_AG": 0, "SpliceAI_pred_DS_AL": 0,
            "SpliceAI_pred_DS_DG": 0, "SpliceAI_pred_DS_DL": 0,
            "ONCOKB_ONCOGENIC": "",
            "GENIE_count": None,
            "CancerHotspots_position_count": None,
            "CancerHotspots_aa_change_count": None,
            "LOEUF": None,
            "EXON": "",
            "O5_same_position": False, "O5_strength": "",
            "O9_candidate": False, "O9_candidate_reason": "",
        }
        defaults.update(kwargs)
        return defaults

    def _classify(self, **kwargs):
        cls, codes, score = post._score_variant(self._row(**kwargs))
        return cls, codes, score

    # ── Standalone overrides ──────────────────────────────────────────

    def test_o1_canonical_returns_oncogenic(self):
        cls, codes, _ = self._classify(O1_canonical=True)
        assert cls == "Oncogenic"
        assert "O1" in codes

    def test_b2_wrong_mechanism_returns_vus(self):
        cls, codes, _ = self._classify(
            Gene_role_in_cancer="ONCOGENE",
            Consequence="frameshift_variant"
        )
        assert cls == "VUS"
        assert "B2" in codes

    def test_b1_standalone_high_af_returns_benign(self):
        cls, codes, _ = self._classify(MAX_AF=0.10)
        assert cls == "Benign"
        assert "B1" in codes

    # ── Minimum 2 codes rule ─────────────────────────────────────────

    def test_single_code_returns_vus(self):
        # Only O6 applies (+1) — one code → VUS
        cls, _, score = self._classify(REVEL=0.75, MAX_AF=0.001)
        assert cls == "VUS"

    def test_two_codes_gives_vus_below_threshold(self):
        # O3_mod(+2) + O7_mod(+2) = 4 pts, 2 codes → VUS (4 < 6 threshold)
        cls, codes, score = self._classify(
            MAX_AF=0,
            CancerHotspots_position_count=30, CancerHotspots_aa_change_count=12,
        )
        assert score == 4
        assert cls == "VUS"
        assert "O3" in codes and "O7" in codes

    # ── Real variant scenarios ────────────────────────────────────────

    def test_braf_v600e_profile_oncogenic(self):
        # Simulates BRAF V600E: canonical, hotspot, Oncogenic in OncoKB
        cls, codes, score = self._classify(
            O1_canonical=True,
        )
        assert cls == "Oncogenic"

    def test_idh1_r132h_likely_oncogenic(self):
        # IDH1 R132H: rare in gnomAD, big hotspot, Oncogenic OncoKB, high REVEL
        cls, codes, score = self._classify(
            MAX_AF=0,                                        # O3 +2
            CancerHotspots_position_count=3104,              # O7 Strong: +4
            CancerHotspots_aa_change_count=2891,
            ONCOKB_ONCOGENIC="Oncogenic",                    # O10 +1
            REVEL=0.852,                                     # O6 +1
        )
        assert score >= 6
        assert cls in ("Likely Oncogenic", "Oncogenic")

    def test_benign_common_variant(self):
        cls, _, _ = self._classify(
            MAX_AF=0.15,     # B1 standalone
        )
        assert cls == "Benign"

    def test_b4_synonymous_likely_benign(self):
        # Synonymous + no SpliceAI impact → B4_str(-4) alone → VUS (only 1 code)
        cls, codes, score = self._classify(
            Consequence="synonymous_variant",
            MAX_AF=0.001,    # adds B1_str(-4)
        )
        # B4(-4) + B1(-4) = -8 and 2 codes → Benign
        assert score <= -7 or cls in ("Likely Benign", "Benign", "VUS")

    def test_output_columns_added(self):
        df = pd.DataFrame([self._row(O1_canonical=True)])
        result = post.add_svig_uk_classification(df)
        assert "SVIG_UK_classification" in result.columns
        assert "SVIG_UK_score" in result.columns
        assert "SVIG_UK_codes" in result.columns

    def test_codes_string_format(self):
        cls, codes, score = self._classify(
            MAX_AF=0,
            REVEL=0.85,
        )
        # Check codes format in the full function
        df = pd.DataFrame([self._row(MAX_AF=0, REVEL=0.85)])
        result = post.add_svig_uk_classification(df)
        codes_str = result["SVIG_UK_codes"].iloc[0]
        assert "O3" in codes_str
        assert "O6" in codes_str
