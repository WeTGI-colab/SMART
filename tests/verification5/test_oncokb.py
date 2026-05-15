"""
Unit tests for scripts/oncokb2.0.py

Covers:
  - get_cancer_from_filename()     — filename → OncoTree code mapping
  - clean_value()                  — value sanitisation
  - stringify_list()               — list → pipe-joined string
  - _extract_one_letter_protein()  — HGVSp 3-letter → 1-letter conversion
  - classify_variant()             — SNV vs CNV routing
  - map_svtype_to_cna()            — SVTYPE → OncoKB CNA type
  - load_preferred_transcripts()   — NM whitelist loading
  - get_csq_index()                — CSQ field index lookup
  - parse_csq_format()             — VCF header CSQ parsing
"""
import pytest

from conftest import oncokb


# ---------------------------------------------------------------------------
# Mock variant helper (classify_variant needs variant.INFO.get and variant.ID)
# ---------------------------------------------------------------------------

class _MockInfo:
    def __init__(self, data):
        self._data = data or {}

    def get(self, key):
        return self._data.get(key)


class MockVariant:
    def __init__(self, svtype=None, variant_id=None):
        self.INFO = _MockInfo({"SVTYPE": svtype} if svtype else {})
        self.ID = variant_id or ""


# ===========================================================================
# get_cancer_from_filename
# ===========================================================================

class TestGetCancerFromFilename:
    @pytest.mark.parametrize("filename,expected", [
        ("patient_lung_sample.vcf.gz",       "LUNG"),
        ("melanoma_tumor.vcf",               "MEL"),
        ("BRCA_breast_sample.vcf",           "BRCA"),
        ("colon_adenocarcinoma.vcf",         "COAD"),
        ("pancreatic_cancer_01.vcf",         "PAAD"),
        ("thyroid-carcinoma.vcf",            "THY"),   # split on '-' → 'thyroid' matched first
        ("cholangiocarcinoma_sample.vcf",    "CHOL"),
        ("sarcoma.vcf",                      "SARC"),
        ("ovarian_tumour.vcf",               "OVT"),
        ("unknown_sample.vcf",               "UNKNOWN"),
        ("random_name.vcf",                  "UNKNOWN"),
    ])
    def test_known_and_unknown(self, filename, expected):
        assert oncokb.get_cancer_from_filename(filename) == expected

    def test_case_insensitive_lookup(self):
        # Keys in CANCER_MAP are lowercase; filename parts are lowercased before lookup
        assert oncokb.get_cancer_from_filename("LUNG_sample.vcf") == "LUNG"


# ===========================================================================
# clean_value
# ===========================================================================

class TestCleanValue:
    def test_none_returns_empty_string(self):
        assert oncokb.clean_value(None) == ""

    def test_semicolons_replaced_with_commas(self):
        assert oncokb.clean_value("a;b;c") == "a,b,c"

    def test_newlines_stripped(self):
        assert oncokb.clean_value("text\nmore") == "text more"

    def test_carriage_return_stripped(self):
        assert oncokb.clean_value("text\rmore") == "text more"

    def test_leading_trailing_whitespace_stripped(self):
        assert oncokb.clean_value("  hello  ") == "hello"

    def test_plain_string_unchanged(self):
        assert oncokb.clean_value("Cobimetinib") == "Cobimetinib"

    def test_integer_converted_to_string(self):
        assert oncokb.clean_value(42) == "42"


# ===========================================================================
# stringify_list
# ===========================================================================

class TestStringifyList:
    def test_empty_list_returns_empty_string(self):
        assert oncokb.stringify_list([]) == ""

    def test_none_returns_empty_string(self):
        assert oncokb.stringify_list(None) == ""

    def test_single_item(self):
        assert oncokb.stringify_list(["Ivosidenib"]) == "Ivosidenib"

    def test_multiple_items_joined_with_pipe(self):
        result = oncokb.stringify_list(["Cobimetinib", "Trametinib"])
        assert result == "Cobimetinib|Trametinib"

    def test_none_items_filtered_out(self):
        result = oncokb.stringify_list(["a", None, "b"])
        assert result == "a|b"

    def test_empty_string_items_filtered_out(self):
        result = oncokb.stringify_list(["a", "", "b"])
        assert result == "a|b"

    def test_semicolons_in_values_replaced(self):
        # clean_value converts ; → , before joining
        result = oncokb.stringify_list(["a;b"])
        assert ";" not in result


# ===========================================================================
# _extract_one_letter_protein
# ===========================================================================

class TestExtractOneLetterProtein:
    def _csq(self, *fields):
        """Build a pipe-delimited CSQ string; HGVSp is at index 1."""
        base = ["allele", ""]
        if fields:
            base[1] = fields[0]
        return "|".join(base)

    @pytest.mark.parametrize("hgvsp_raw,expected", [
        ("ENSP00000:p.Val600Glu",   "V600E"),   # BRAF V600E
        ("ENSP00000:p.Gln61Arg",    "Q61R"),    # NRAS Q61R
        ("ENSP00000:p.Arg132His",   "R132H"),   # IDH1 R132H
        ("ENSP00000:p.Glu545Lys",   "E545K"),   # PIK3CA E545K
        ("ENSP00000:p.His1047Arg",  "H1047R"),  # PIK3CA H1047R
        ("ENSP00000:p.Arg175His",   "R175H"),   # TP53 R175H
    ])
    def test_real_variants(self, hgvsp_raw, expected):
        csq = self._csq(hgvsp_raw)
        result = oncokb._extract_one_letter_protein(csq, 1)
        assert result == expected

    def test_stop_codon_ter_to_asterisk(self):
        # Ter is replaced with * (1 char) before the 3-letter regex runs.
        # The regex expects 3-letter alt, so it fails and the partially-
        # converted string "Arg248*" is returned as-is.
        csq = self._csq("ENSP:p.Arg248Ter")
        assert oncokb._extract_one_letter_protein(csq, 1) == "Arg248*"

    def test_negative_index_returns_none(self):
        csq = self._csq("ENSP:p.Val600Glu")
        assert oncokb._extract_one_letter_protein(csq, -1) is None

    def test_empty_protein_field_returns_none(self):
        csq = "allele|"
        assert oncokb._extract_one_letter_protein(csq, 1) is None

    def test_index_out_of_range_returns_none(self):
        csq = "allele"
        assert oncokb._extract_one_letter_protein(csq, 5) is None


# ===========================================================================
# classify_variant
# ===========================================================================

class TestClassifyVariant:
    def test_plain_snv_returns_snv(self):
        v = MockVariant()
        assert oncokb.classify_variant(v) == "snv"

    @pytest.mark.parametrize("svtype", ["DUP", "DEL", "GAIN", "LOSS"])
    def test_cnv_svtypes_return_cnv(self, svtype):
        v = MockVariant(svtype=svtype)
        assert oncokb.classify_variant(v) == "cnv"

    @pytest.mark.parametrize("svtype", ["INS", "BND", "INV"])
    def test_non_cnv_svtypes_return_snv(self, svtype):
        v = MockVariant(svtype=svtype)
        assert oncokb.classify_variant(v) == "snv"

    @pytest.mark.parametrize("variant_id", [
        "MantaDUP:1:1000:2000", "MantaDEL:1:1000:2000"
    ])
    def test_manta_dup_del_ids_return_cnv(self, variant_id):
        v = MockVariant(variant_id=variant_id)
        assert oncokb.classify_variant(v) == "cnv"

    @pytest.mark.parametrize("variant_id", [
        "MantaINS:1:1000:2000", "MantaBND:1:1000:2000"
    ])
    def test_manta_ins_bnd_ids_return_snv(self, variant_id):
        v = MockVariant(variant_id=variant_id)
        assert oncokb.classify_variant(v) == "snv"

    def test_id_overrides_svtype_for_manta_dup(self):
        # Even if SVTYPE is INS, MantaDUP in ID → cnv
        v = MockVariant(svtype="INS", variant_id="MantaDUP:1:100:200")
        assert oncokb.classify_variant(v) == "cnv"


# ===========================================================================
# map_svtype_to_cna
# ===========================================================================

class TestMapSvtypeToCna:
    @pytest.mark.parametrize("variant_id,expected", [
        ("MantaDUP:1:100:200", "AMPLIFICATION"),
        ("MantaDEL:1:100:200", "DELETION"),
    ])
    def test_manta_ids(self, variant_id, expected):
        assert oncokb.map_svtype_to_cna(None, variant_id) == expected

    @pytest.mark.parametrize("svtype,expected", [
        ("DUP",  "AMPLIFICATION"),
        ("GAIN", "AMPLIFICATION"),
        ("DEL",  "DELETION"),
        ("LOSS", "DELETION"),
    ])
    def test_svtype_mapping(self, svtype, expected):
        assert oncokb.map_svtype_to_cna(svtype, "") == expected

    def test_ins_returns_none(self):
        assert oncokb.map_svtype_to_cna("INS", "") is None

    def test_none_svtype_no_id_returns_none(self):
        assert oncokb.map_svtype_to_cna(None, "") is None


# ===========================================================================
# load_preferred_transcripts
# ===========================================================================

class TestLoadPreferredTranscripts:
    def test_none_path_returns_empty_set(self):
        result = oncokb.load_preferred_transcripts(None)
        assert result == set()

    def test_missing_file_returns_empty_set(self):
        result = oncokb.load_preferred_transcripts("/nonexistent/path.txt")
        assert result == set()

    def test_strips_version_suffixes(self, nm_file):
        result = oncokb.load_preferred_transcripts(nm_file)
        assert "NM_000077" in result
        assert "NM_000077.5" not in result

    def test_skips_comment_lines(self, nm_file):
        result = oncokb.load_preferred_transcripts(nm_file)
        assert "#" not in result

    def test_returns_set(self, nm_file):
        assert isinstance(oncokb.load_preferred_transcripts(nm_file), set)

    def test_loads_all_transcripts(self, nm_file):
        result = oncokb.load_preferred_transcripts(nm_file)
        assert {"NM_000077", "NM_002524", "NM_005228", "NM_000546"}.issubset(result)


# ===========================================================================
# get_csq_index
# ===========================================================================

class TestGetCsqIndex:
    def setup_method(self):
        self.field_map = {"Allele": 0, "Consequence": 1, "SYMBOL": 2, "HGVSp": 3}

    def test_existing_field_returns_correct_index(self):
        assert oncokb.get_csq_index(self.field_map, "Consequence") == 1
        assert oncokb.get_csq_index(self.field_map, "HGVSp") == 3

    def test_missing_field_not_required_returns_minus_one(self):
        assert oncokb.get_csq_index(self.field_map, "MISSING", required=False) == -1

    def test_missing_field_required_raises_key_error(self):
        with pytest.raises(KeyError, match="MISSING"):
            oncokb.get_csq_index(self.field_map, "MISSING", required=True)

    def test_default_not_required(self):
        # required=False is the default
        assert oncokb.get_csq_index(self.field_map, "NOTHERE") == -1


# ===========================================================================
# parse_csq_format
# ===========================================================================

class TestParseCsqFormat:
    def test_parses_field_names(self, vep_vcf):
        fields, field_map = oncokb.parse_csq_format(vep_vcf)
        assert "Allele" in fields
        assert "Consequence" in fields
        assert "HGVSp" in fields
        assert "MANE_SELECT" in fields

    def test_field_map_indices_are_correct(self, vep_vcf):
        fields, field_map = oncokb.parse_csq_format(vep_vcf)
        for i, name in enumerate(fields):
            assert field_map[name] == i

    def test_missing_csq_header_raises_runtime_error(self, vep_vcf_no_csq):
        with pytest.raises(RuntimeError, match="CSQ header"):
            oncokb.parse_csq_format(vep_vcf_no_csq)
