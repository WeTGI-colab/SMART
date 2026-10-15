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


class TestStripEmbeddedWhitespace:
    def test_newline_replaced_with_space(self):
        df = pd.DataFrame({"ONCOKB_EFFECT_DESC": ["damaging.\nIn silico predicts (PMID: 1)."]})
        result = post.strip_embedded_whitespace(df)
        assert "\n" not in result["ONCOKB_EFFECT_DESC"].iloc[0]
        assert result["ONCOKB_EFFECT_DESC"].iloc[0] == "damaging. In silico predicts (PMID: 1)."

    def test_tab_and_cr_replaced(self):
        df = pd.DataFrame({"FIELD": ["a\tb\r\nc"]})
        result = post.strip_embedded_whitespace(df)
        assert result["FIELD"].iloc[0] == "a b c"

    def test_normal_values_unchanged(self):
        df = pd.DataFrame({"SYMBOL": ["BRAF"], "HGVSp_Short": ["p.V600E"]})
        result = post.strip_embedded_whitespace(df)
        assert result["SYMBOL"].iloc[0] == "BRAF"
        assert result["HGVSp_Short"].iloc[0] == "p.V600E"

    def test_row_does_not_split_when_written(self):
        # The actual bug: an embedded newline splits a MAF record into a phantom
        # row with a blank Tumor_Sample_Barcode. After sanitising, N records must
        # still write as exactly N physical lines.
        import io
        df = pd.DataFrame({
            "Tumor_Sample_Barcode": ["S1", "S2"],
            "ONCOKB_EFFECT_DESC": ["ok", "text.\nIn silico (PMID: 1).\t['1']"],
        })
        clean = post.strip_embedded_whitespace(df)
        buf = io.StringIO()
        clean.to_csv(buf, sep="\t", index=False, header=False)
        assert buf.getvalue().rstrip("\n").count("\n") + 1 == 2


# ===========================================================================
# calculate_end_position
# ===========================================================================

class TestCalculateEndPosition:
    """End_Position: the last reference base the record spans.

    These tests used to assert the previous behaviour, which produced an **empty**
    End_Position on 44 of 44 rows of every real run -- the Tier 1 MAF the laboratory reads,
    copy-number calls included. Three things were wrong and all three are asserted against
    here now:

      * an end was computed only when `ID == "."`, so a caller that names its records got
        nothing. Manta names its records; so does any VCF carrying an rsID.
      * `DUPSVLEN` was read for the Manta rows. No caller writes it -- Manta and DRAGEN-SV
        write `SVLEN` -- and it is not a column of this table, so those rows were blanked
        even though `END` sat correct in the same row.
      * an insertion's end was `start + len(ALT)`, which overstates the reference span.

    The convention is `start + len(REF) - 1`, which is what `oncokb_maf_query.py` already
    computes for the OncoKB genomic query and is the VCF definition: REF carries the anchor
    base, so its length is the number of reference bases covered.
    """

    def _df(self, **kwargs) -> pd.DataFrame:
        defaults = {"Start_Position": 1000, "ID": ".", "REF": "A", "ALT": "G"}
        defaults.update(kwargs)
        return pd.DataFrame([defaults])

    def _end(self, **kwargs):
        v = post.calculate_end_position(self._df(**kwargs))["End_Position"].iloc[0]
        return None if v is None or v != v else int(v)

    def test_snv_end_equals_start(self):
        assert self._end(REF="A", ALT="G") == 1000

    def test_a_named_record_is_computed_too(self):
        """The defect that emptied the column: an ID other than '.' was not a reason to stop."""
        assert self._end(ID="NRAS_Q61R", REF="A", ALT="G") == 1000
        assert self._end(ID="rs113488022", REF="A", ALT="G") == 1000

    def test_deletion_spans_the_deleted_bases(self):
        assert self._end(REF="AG", ALT="A") == 1001

    def test_insertion_spans_only_the_anchor_base(self):
        """`A>ACTG` inserts after position 1000; the record covers one reference base.

        The previous behaviour returned 1004, overstating the span by len(ALT)-1.
        """
        assert self._end(REF="A", ALT="ACTG") == 1000

    def test_mnv_spans_its_own_length(self):
        assert self._end(REF="CA", ALT="GT") == 1001

    def test_the_callers_end_wins_for_a_symbolic_alt(self):
        """A copy-number call has no sequence to measure, so END is the only authority.

        This is the real MET amplification from the verification cohort, whose segment end
        was present and discarded.
        """
        assert self._end(ID="MantaDUP:MET_AMP", Start_Position=116672196,
                         REF="N", ALT="<DUP>", END="116798037") == 116798037

    def test_svlen_is_used_when_there_is_no_end(self):
        """SVLEN is what Manta and DRAGEN-SV actually write, and it may be negative."""
        assert self._end(ID="MantaDEL:1", REF="N", ALT="<DEL>", SVLEN="-500") == 1500
        assert self._end(ID="MantaDUP:1", REF="N", ALT="<DUP>", SVLEN="500") == 1500

    def test_svlen_is_not_applied_to_a_sequence_alt(self):
        """A record with real alleles has its span in REF; SVLEN would double-count."""
        assert self._end(ID="x", REF="AG", ALT="A", SVLEN="500") == 1001

    def test_a_breakend_has_no_span(self):
        """A junction, not an interval. Its REF is a placeholder, which would otherwise
        give end == start and read as a point event."""
        assert self._end(ID="MantaBND:1", REF="N", ALT="N[chr9:133729451[") is None
        assert self._end(ID="MantaBND:1", REF="N", ALT="]chr9:133729451]N") is None

    def test_a_placeholder_ref_gives_no_span(self):
        assert self._end(ID="x", REF=".", ALT="<DUP>") is None

    @pytest.mark.parametrize("end", ["", ".", "nan", "None", None])
    def test_an_unusable_end_falls_through_to_the_reference_span(self, end):
        assert self._end(ID="x", REF="AG", ALT="A", END=end) == 1001

    def test_an_end_before_the_start_is_not_trusted(self):
        """A malformed END must not produce an interval that runs backwards."""
        assert self._end(ID="x", Start_Position=1000, REF="AG", ALT="A", END="900") == 1001

    def test_a_non_numeric_start_gives_nothing(self):
        assert self._end(Start_Position="chr1", REF="A", ALT="G") is None


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


class TestOncoKBTreatmentsShapes:
    """_apply_oncokb_fixes parses ONCOKB_treatments and must survive every shape it meets.

    json.loads and ast.literal_eval can both succeed and still not return a list of dicts: a
    quoted string parses to a str, and iterating a str yields characters, so `t.get(...)` raised
    AttributeError and aborted post_analysis for the whole run. That happened on a CNA row whose
    treatments had already been flattened to drug names. The parse is now validated rather than
    trusted, and these pin it.
    """

    @staticmethod
    def _run(value):
        df = pd.DataFrame([{
            # ONCOKB_QUERY_TYPE == "CNA" is what selects the rows _extract_levels runs on.
            # Without it the mask is empty and these tests would pass without executing the
            # code they exist to cover.
            "ONCOKB_QUERY_TYPE": "CNA",
            "Variant_Type": "INS", "Tumor_Seq_Allele2": "<DUP>",
            "ONCOKB_treatments": value,
            "LEVEL_1": "", "LEVEL_2": "", "LEVEL_3A": "",
            "HIGHEST_SENSITIVE_LEVEL": "", "ONCOKB_SENS_LVL": "",
            "HIGHEST_RESISTANCE_LEVEL": "", "ONCOKB_RESIST_LVL": "",
            "VARIANT_IN_ONCOKB": "", "ONCOKB_variantExist": "",
        }])
        out, _, _ = post._apply_oncokb_fixes(df)
        return out

    @pytest.mark.parametrize("value", [
        "'Erlotinib,Gefitinib,Osimertinib'",   # parses to a str — the case that aborted a run
        '"a string"',                           # json-quoted string
        "Erlotinib,Gefitinib",                  # bare text, parses to nothing
        "[]", "", "nan", "{}",                  # empty forms
        "[1, 2, 3]",                            # a list, but not of dicts
        "['a', 'b']",                           # a list of strings
        "[{'level': 'LEVEL_1'}]",               # a dict with no drugs
        "[{'level': 'LEVEL_1', 'drugs': 'not a list'}]",
        "[{'level': 'LEVEL_1', 'drugs': ['not a dict']}]",
    ])
    def test_no_shape_raises(self, value):
        self._run(value)      # the assertion is that it returns at all

    def test_a_well_formed_value_still_extracts_the_drug(self):
        out = self._run("[{'level': 'LEVEL_1', 'drugs': [{'drugName': 'Osimertinib'}]}]")
        assert out.at[0, "LEVEL_1"] == "Osimertinib"

    def test_a_double_encoded_value_is_parsed_through(self):
        """CDK4 and ERBB2 arrive as a Python repr wrapped in a JSON string.

        One parse yields the inner string, not the list. Iterating that string raised
        AttributeError; stopping at it would silently lose every treatment level for those
        rows, which is the state the suite caught as 4 empty LEVEL columns.
        """
        out = self._run('"[{\'level\': \'LEVEL_1\', \'drugs\': [{\'drugName\': \'Osimertinib\'}]}]"')
        assert out.at[0, "LEVEL_1"] == "Osimertinib"

    def test_gene_in_oncokb_is_taken_from_the_gene_query_for_cna_rows(self):
        """MafAnnotator cannot annotate a CNA, so it leaves GENE_IN_ONCOKB False.

        That reads as "OncoKB does not know this gene" for ERBB2, MET, CDK4 and CDKN2A.
        oncokb2.0.py asked separately and recorded the answer in ONCOKB_geneExist.
        """
        import pandas as _pd
        df = _pd.DataFrame([{
            "ONCOKB_QUERY_TYPE": "CNA", "Variant_Type": "INS",
            "Tumor_Seq_Allele2": "<DUP>", "ONCOKB_treatments": "[]",
            "GENE_IN_ONCOKB": "False", "ONCOKB_geneExist": "True",
            "VARIANT_IN_ONCOKB": "", "ONCOKB_variantExist": "True",
        }])
        out, _, _ = post._apply_oncokb_fixes(df)
        assert out.at[0, "GENE_IN_ONCOKB"] == "True"
        assert out.at[0, "VARIANT_IN_ONCOKB"] == "True"

    def test_duplicate_drugs_are_not_repeated(self):
        out = self._run("[{'level': 'LEVEL_1', 'drugs': [{'drugName': 'X'}, "
                        "{'drugName': 'X'}, {'drugName': 'Y'}]}]")
        assert out.at[0, "LEVEL_1"] == "X,Y"


class TestCNALevelsFromStructuredColumns:
    """CNA treatment levels come from ONCOKB_TX_{n}_level / _drugs, not from the blob.

    ONCOKB_treatments holds the API response stringified, and its free text carries unescaped
    apostrophes that close a Python literal early, so it cannot be parsed back for most rows.
    oncokb2.0.py has already written the same information out per entry from the real JSON.
    Reading those columns needs no parsing, which is why this works where re-parsing could not.
    """

    def _row(self, **extra):
        base = {
            "ONCOKB_QUERY_TYPE": "CNA", "Variant_Type": "INS",
            "Tumor_Seq_Allele2": "<DUP>", "ONCOKB_treatments": "[]",
            "LEVEL_1": "", "LEVEL_2": "", "LEVEL_3A": "", "LEVEL_4": "",
            "HIGHEST_SENSITIVE_LEVEL": "", "ONCOKB_SENS_LVL": "",
            "HIGHEST_RESISTANCE_LEVEL": "", "ONCOKB_RESIST_LVL": "",
            "VARIANT_IN_ONCOKB": "", "ONCOKB_variantExist": "",
            "GENE_IN_ONCOKB": "False", "ONCOKB_geneExist": "True",
        }
        base.update(extra)
        out, _, _ = post._apply_oncokb_fixes(pd.DataFrame([base]))
        return out.iloc[0]

    def test_a_single_entry_fills_its_level(self):
        r = self._row(**{"ONCOKB_TX_0_level": "LEVEL_1",
                         "ONCOKB_TX_0_drugs": "Ado-Trastuzumab Emtansine"})
        assert r["LEVEL_1"] == "Ado-Trastuzumab Emtansine"

    def test_several_entries_at_one_level_are_joined(self):
        r = self._row(**{"ONCOKB_TX_0_level": "LEVEL_4", "ONCOKB_TX_0_drugs": "Palbociclib",
                         "ONCOKB_TX_1_level": "LEVEL_4", "ONCOKB_TX_1_drugs": "Abemaciclib"})
        assert r["LEVEL_4"] == "Palbociclib,Abemaciclib"

    def test_a_combination_is_split_into_its_drugs(self):
        # OncoKB writes a combination as "A + B"; the levels list each drug
        r = self._row(**{"ONCOKB_TX_0_level": "LEVEL_1",
                         "ONCOKB_TX_0_drugs": "Lapatinib + Capecitabine"})
        assert r["LEVEL_1"] == "Lapatinib,Capecitabine"

    def test_a_drug_is_not_repeated_across_entries(self):
        r = self._row(**{"ONCOKB_TX_0_level": "LEVEL_2", "ONCOKB_TX_0_drugs": "Trastuzumab",
                         "ONCOKB_TX_1_level": "LEVEL_2", "ONCOKB_TX_1_drugs": "Trastuzumab"})
        assert r["LEVEL_2"] == "Trastuzumab"

    def test_entries_go_to_their_own_levels(self):
        r = self._row(**{"ONCOKB_TX_0_level": "LEVEL_1", "ONCOKB_TX_0_drugs": "A",
                         "ONCOKB_TX_1_level": "LEVEL_3A", "ONCOKB_TX_1_drugs": "B"})
        assert r["LEVEL_1"] == "A" and r["LEVEL_3A"] == "B"

    def test_an_unknown_level_is_ignored_rather_than_crashing(self):
        r = self._row(**{"ONCOKB_TX_0_level": "LEVEL_NONSENSE", "ONCOKB_TX_0_drugs": "A"})
        assert r["LEVEL_1"] == ""

    def test_an_empty_drug_entry_is_skipped(self):
        for empty in ("", "nan", "."):
            r = self._row(**{"ONCOKB_TX_0_level": "LEVEL_1", "ONCOKB_TX_0_drugs": empty})
            assert r["LEVEL_1"] == ""

    def test_an_existing_value_is_not_overwritten(self):
        # MafAnnotator's own answer wins where it gave one; this only fills gaps
        r = self._row(LEVEL_1="Already there",
                      **{"ONCOKB_TX_0_level": "LEVEL_1", "ONCOKB_TX_0_drugs": "Other"})
        assert r["LEVEL_1"] == "Already there"

    def test_nan_cells_are_treated_as_empty(self):
        """post_analysis reads the MAF with dtype=str, so an empty cell arrives as NaN.

        NaN is truthy in Python and str(NaN) is "nan", so `v or ""` keeps the NaN and a naive
        emptiness test concludes the cell is already filled. That silently skipped every row
        this was meant to fill: the unit tests passed because they built the frame from strings,
        while the pipeline produced NaN and the levels stayed empty in the real output.
        """
        import numpy as _np
        import pandas as _pd
        df = _pd.DataFrame([{
            "ONCOKB_QUERY_TYPE": "CNA", "Variant_Type": "INS",
            "Tumor_Seq_Allele2": "<DUP>", "ONCOKB_treatments": "[]",
            "LEVEL_1": _np.nan, "LEVEL_3A": _np.nan,
            "GENE_IN_ONCOKB": "False", "ONCOKB_geneExist": "True",
            "VARIANT_IN_ONCOKB": _np.nan, "ONCOKB_variantExist": "True",
            "ONCOKB_TX_0_level": "LEVEL_1",
            "ONCOKB_TX_0_drugs": "Ado-Trastuzumab Emtansine",
            "ONCOKB_TX_1_level": "LEVEL_3A",
            "ONCOKB_TX_1_drugs": "Trastuzumab Deruxtecan",
        }])
        out, _, _ = post._apply_oncokb_fixes(df)
        assert out.at[0, "LEVEL_1"] == "Ado-Trastuzumab Emtansine"
        assert out.at[0, "LEVEL_3A"] == "Trastuzumab Deruxtecan"

    def test_no_tx_columns_at_all_is_survivable(self):
        r = self._row()
        assert r["LEVEL_1"] == ""
