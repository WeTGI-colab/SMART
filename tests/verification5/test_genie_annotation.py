"""
Unit tests for scripts/genie_annotation.py — GENIE case counts as annotation columns.
"""
import importlib.util
import os
import sys

SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "scripts")
sys.path.insert(0, SCRIPTS_DIR)
_spec = importlib.util.spec_from_file_location("genie_annotation", os.path.join(SCRIPTS_DIR, "genie_annotation.py"))
ga = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ga)
svig = ga.svig


def entry(n, lineage, code, vclass="missense", pos=100, alt="A"):
    tissue = "|".join({"Myeloid": "Myeloid", "Lymphoid": "Lymphoid"}.get(l, "Lung") for l in lineage.split("|"))
    return svig.GenieEntry(vclass, "chr1", pos, pos, "G", alt, "R", "H", lineage, n, code, tissue)


def test_split_counts_counts_patients_once_per_group():
    entries = [entry(5, "Myeloid", "AML"), entry(3, "Solid", "LUAD"), entry(2, "Myeloid|Solid", "AML|LUAD")]
    # a patient with haematological and solid samples counts in both groups but once in all
    assert ga.split_counts(entries) == (10, 7, 5)


def test_top_tumour_types():
    entries = [entry(5, "Myeloid", "AML"), entry(3, "Solid", "LUAD"), entry(2, "Myeloid|Solid", "AML|LUAD")]
    assert ga.top_tumour_types(entries) == "AML 7; LUAD 5"


def test_cohort_label():
    sizes = {("all", "all"): 242866, ("lineage", "Haematological"): 24926, ("lineage", "Solid"): 215540}
    assert ga.cohort_label(sizes) == "all 242,866; haematological 24,926; solid 215,540 patients"


def test_annotate_row_same_nucleotide_and_no_o4_class():
    res = svig.Resources()
    res.genie = type("G", (), {"by_gene": {"GENE1": [entry(4, "Myeloid", "AML"), entry(9, "Solid", "LUAD", alt="T")]}})()
    row = {"Hugo_Symbol": "GENE1", "Chromosome": "chr1", "POS": "100", "REF": "G", "ALT": "A",
           "Consequence": "synonymous_variant", "HGVSp": ""}
    out = ga.annotate_row(row, res, "cohort")
    assert (out["GENIE_SAME_NUCLEOTIDE_ALL"], out["GENIE_SAME_NUCLEOTIDE_HAEM"], out["GENIE_SAME_NUCLEOTIDE_SOLID"]) == (4, 4, 0)
    assert out["GENIE_O4_MATCH_TYPE"] == "" and out["GENIE_TOP_TUMOUR_TYPES"] == "AML 4"


def test_missing_table_is_skipped(tmp_path, capsys):
    maf = tmp_path / "s.maf"
    maf.write_text("Hugo_Symbol\nGENE1\n")
    assert ga.main(["--maf", str(maf), "--refs", str(tmp_path)]) == 0
    assert "skipped" in capsys.readouterr().out
    assert maf.read_text() == "Hugo_Symbol\nGENE1\n"
