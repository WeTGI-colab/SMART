"""
Unit tests for scripts/oncokb_maf_query.py — splitting MafAnnotator queries between
protein change and genomic change, and merging the outputs back without losing rows.
"""
import importlib.util
import os

import pytest

SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "scripts")
_spec = importlib.util.spec_from_file_location("oncokb_maf_query", os.path.join(SCRIPTS_DIR, "oncokb_maf_query.py"))
mq = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mq)

COLUMNS = ["Hugo_Symbol", "NM_Transcript", "Chromosome", "Start_Position", "End_Position",
           "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "HGVSp_Short"]
ROWS = [
    ["BRAF", "NM_004333", "chr7", "140753336", "140753336", "A", "A", "T", "p.V600E"],
    ["NF1", "NM_001042492", "chr17", "31159042", "31159043", "A", "A", "AT", "p.Y80Lfs*27"],
    ["NF1", "NM_000267", "chr17", "31327535", "31327535", "C", "C", "T", "p.R1748*"],
]


def genomic(gene, nm, transcripts):
    return bool(transcripts.get(gene) and nm.split(".")[0] != transcripts[gene])


@pytest.fixture
def table(tmp_path):
    path = tmp_path / "table.tsv"
    path.write_text("\t".join(COLUMNS) + "\n" + "".join("\t".join(r) + "\n" for r in ROWS))
    return path


def test_split_by_transcript(table):
    _, protein, gen = mq.split(str(table), {"NF1": "NM_000267"}, genomic)
    assert [r["HGVSp_Short"] for r in protein] == ["p.V600E", "p.R1748*"]
    assert [r["HGVSp_Short"] for r in gen] == ["p.Y80Lfs*27"]
    # VCF-style insertion: End_Position made consistent with the reference allele for the query
    assert gen[0]["End_Position"] == "31159042"


def test_merge_restores_order_positions_and_records_query(table, tmp_path):
    fields, protein, gen = mq.split(str(table), {"NF1": "NM_000267"}, genomic)
    for name, rows in (("p.maf", protein), ("g.maf", gen)):
        annotated = [dict(r, ONCOGENIC="Likely Oncogenic") for r in rows]
        mq.write_table(str(tmp_path / name), fields + ["ONCOGENIC"], annotated)
    names, rows = mq.merge(str(tmp_path / "p.maf"), str(tmp_path / "g.maf"))
    assert [r["HGVSp_Short"] for r in rows] == [r[8] for r in ROWS]
    assert rows[1]["End_Position"] == "31159043"
    assert rows[1][mq.QUERY_COLUMN].startswith("Genomic_Change")
    assert rows[0][mq.QUERY_COLUMN] == "HGVSp_Short"
    assert mq.ROW_ID not in names and mq.ORIGINAL_END not in names and names[-1] == mq.QUERY_COLUMN


def test_merge_labels_protein_fallback(table, tmp_path):
    fields, protein, gen = mq.split(str(table), {"NF1": "NM_000267"}, genomic)
    mq.write_table(str(tmp_path / "p.maf"), fields, protein)
    mq.write_table(str(tmp_path / "g.maf"), fields, gen)
    _, rows = mq.merge(str(tmp_path / "p.maf"), str(tmp_path / "g.maf"), genomic_failed=True)
    assert "genomic-change query failed" in rows[1][mq.QUERY_COLUMN]


def test_merge_fails_when_rows_are_lost(table, tmp_path):
    fields, protein, _ = mq.split(str(table), {"NF1": "NM_000267"}, genomic)
    mq.write_table(str(tmp_path / "p.maf"), fields, protein)
    with pytest.raises(ValueError, match="missing"):
        mq.merge(str(tmp_path / "p.maf"), str(tmp_path / "absent.maf"))


def test_merge_reassembles_rows_split_by_oncokb_line_breaks(table, tmp_path):
    """OncoKB free-text fields can contain line breaks; MafAnnotator writes them unquoted."""
    fields, protein, gen = mq.split(str(table), {"NF1": "NM_000267"}, genomic)
    names = fields + ["ONCOKB_VARIANT_SUMMARY"]
    body = []
    for i, r in enumerate(protein):
        summary = "First line.\nSecond line." if i == 0 else "ok"
        body.append("\t".join([r.get(c, "") for c in fields] + [summary]))
    (tmp_path / "p.maf").write_text("\t".join(names) + "\n" + "\n".join(body) + "\n")
    mq.write_table(str(tmp_path / "g.maf"), names, [dict(r, ONCOKB_VARIANT_SUMMARY="ok") for r in gen])
    _, rows = mq.merge(str(tmp_path / "p.maf"), str(tmp_path / "g.maf"))
    assert len(rows) == len(ROWS)
    assert rows[0]["ONCOKB_VARIANT_SUMMARY"] == "First line. Second line."


def test_written_cells_have_no_line_breaks(tmp_path):
    mq.write_table(str(tmp_path / "o.tsv"), ["a", "b"], [{"a": "x\ny", "b": None}])
    assert (tmp_path / "o.tsv").read_text() == "a\tb\nx y\t\n"


def test_split_reassembles_input_rows_split_by_line_breaks(tmp_path):
    """A line break inside a cell of the input table must not create or lose a variant."""
    path = tmp_path / "table.tsv"
    broken = list(ROWS[0])
    broken[8] = "p.V600E"
    body = ["\t".join(broken[:8]) + "\tFirst line.\nSecond line."] + ["\t".join(r) for r in ROWS[1:]]
    path.write_text("\t".join(COLUMNS) + "\n" + "\n".join(body) + "\n")
    _, protein, gen = mq.split(str(path), {"NF1": "NM_000267"}, genomic)
    assert len(protein) + len(gen) == len(ROWS)
    assert protein[0]["HGVSp_Short"] == "First line. Second line."


def test_short_line_that_starts_a_variant_is_not_a_continuation(tmp_path):
    """A variant split across lines still starts with the genome build, as SMART's tables do;
    a continuation line does not, so the two cases can be told apart."""
    path = tmp_path / "table.tsv"
    cols = ["NCBI_Build"] + COLUMNS
    rows = [["GRCh38"] + r for r in ROWS]
    # the line break falls inside a cell (as OncoKB free text does), not at a tab
    split_variant = ("\t".join(rows[1][:3]) + "\tchr\n17\t" + "\t".join(rows[1][4:]))
    body = "\t".join(rows[0]) + "\n" + split_variant + "\n" + "\t".join(rows[2])
    path.write_text("\t".join(cols) + "\n" + body + "\n")
    names, out = mq.read_annotated(str(path))
    assert len(out) == 3
    assert [r["Hugo_Symbol"] for r in out] == [r[0] for r in ROWS]
    assert out[1]["Chromosome"] == "chr 17" and out[1]["HGVSp_Short"] == ROWS[1][8]
