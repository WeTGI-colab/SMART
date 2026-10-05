"""
Unit tests for utils/build_dbnsfp.py

Covers:
  - the predictor list: no tool already annotated by SMART, no repeated tool
  - select_columns(): header lookup and clear errors for missing columns
  - chromosome_files(): chromosome ordering for a tabix-ready concatenation
"""
import importlib.util
import os
import re

import pytest

UTILS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "utils")
_spec = importlib.util.spec_from_file_location("build_dbnsfp", os.path.join(UTILS_DIR, "build_dbnsfp.py"))
build = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(build)


def tool_of(column):
    """Tool name of a dbNSFP column: 'AlphaMissense_score' -> 'AlphaMissense'."""
    return re.sub(r"_(score|pred|Confidence)$", "", column).replace("-phred_coding", "")


class TestPredictorList:
    def test_no_predictor_already_annotated_by_smart(self):
        # SMART annotates SIFT and PolyPhen (VEP), REVEL and CADD (plugins)
        for column in build.PREDICTOR_COLUMNS:
            assert not re.match(r"^(SIFT|SIFT4G|Polyphen2|REVEL|CADD)_", column), column

    def test_no_repeated_columns(self):
        assert len(build.PREDICTOR_COLUMNS) == len(set(build.PREDICTOR_COLUMNS))

    def test_one_version_per_tool(self):
        tools = {tool_of(c) for c in build.PREDICTOR_COLUMNS}
        for variant, base in [("VARITY_ER", "VARITY_R"), ("BayesDel_addAF", "BayesDel_noAF"),
                              ("Eigen-PC", "Eigen"), ("GERP_91_mammals", "GERP++_RS")]:
            assert variant not in tools
        assert not any("LOO" in c or "rankscore" in c for c in build.PREDICTOR_COLUMNS)

    def test_one_alignment_per_conservation_method(self):
        assert [c for c in build.PREDICTOR_COLUMNS if c.startswith("phyloP")] == ["phyloP100way_vertebrate"]
        assert [c for c in build.PREDICTOR_COLUMNS if c.startswith("phastCons")] == ["phastCons100way_vertebrate"]


class TestSelectColumns:
    def header(self, extra=()):
        names = build.PLUGIN_COLUMNS + ["SIFT_score", "REVEL_score"] + build.PREDICTOR_COLUMNS + list(extra)
        return "\t".join(names) + "\n", names

    def test_indices_follow_the_header(self):
        header, names = self.header()
        idx = build.select_columns(header)
        assert [names[i - 1] for i in idx] == build.PLUGIN_COLUMNS + build.PREDICTOR_COLUMNS

    def test_missing_column_is_reported(self):
        header, _ = self.header()
        with pytest.raises(ValueError, match="AlphaMissense_score"):
            build.select_columns(header.replace("AlphaMissense_score", "AM_score"))


def test_chromosome_files_are_ordered(tmp_path):
    for c in ["10", "2", "X", "1", "M"]:
        (tmp_path / f"dbNSFP5.0a_variant.chr{c}.gz").write_bytes(b"")
    names = [os.path.basename(p) for p in build.chromosome_files(str(tmp_path))]
    assert names == [f"dbNSFP5.0a_variant.chr{c}.gz" for c in ["1", "2", "10", "X", "M"]]
