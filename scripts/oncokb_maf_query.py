#!/usr/bin/env python3

# =============================================================================
# oncokb_maf_query.py
# =============================================================================
#
# DESCRIPTION:
#   Chooses how OncoKB's MafAnnotator queries each variant. OncoKB curates protein
#   changes on one RefSeq transcript per gene; when the transcript selected by SMART
#   is a different one, the protein numbering can differ (NF1 NM_001042492 p.R1769*
#   is NM_000267 p.R1748*) and a protein-change query does not recognise the variant.
#
#     split  writes rows whose selected transcript is OncoKB's (or whose gene OncoKB
#            does not curate) to a table queried by HGVSp_Short, and the other rows
#            to a table queried by genomic change (GRCh38).
#     merge  joins the two MafAnnotator outputs back in the original row order and adds
#            ONCOKB_MAF_QUERY recording the query used for each row.
#
# USAGE:
#   oncokb_maf_query.py split --table T.tsv --token TOKEN --protein P.tsv --genomic G.tsv
#   oncokb_maf_query.py merge --protein P.maf --genomic G.maf --out OUT.maf
# =============================================================================

import argparse
import csv
import importlib.util
import os
import sys

ROW_ID = "SMART_ROW_ID"
ORIGINAL_END = "SMART_ORIGINAL_END_POSITION"
QUERY_COLUMN = "ONCOKB_MAF_QUERY"


def _oncokb_module():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "oncokb2.0.py")
    spec = importlib.util.spec_from_file_location("oncokb2_0", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_table(path):
    with open(path, newline="") as fh:
        lines = [l for l in fh if not l.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    return reader.fieldnames or [], list(reader)


def read_annotated(path):
    """Read a tab-separated table whose cells may contain unquoted line breaks (OncoKB free text
    reaches both the input table and the MafAnnotator output), reassembling split variants."""
    with open(path, newline="") as fh:
        lines = [l.rstrip("\r\n") for l in fh if not l.startswith("#")]
    fieldnames = lines[0].split("\t")
    rows, pending, starts = [], [], set()
    for line in lines[1:]:
        fields = line.split("\t")
        if pending:                       # completing a variant split across lines
            pending[-1] += " " + fields[0]
            pending += fields[1:]
        elif len(fields) < len(fieldnames) and rows and fields[0] not in starts:
            # line break inside the last column of the previous variant
            rows[-1][fieldnames[-1]] += " " + line
            continue
        else:
            pending = fields
        if len(pending) >= len(fieldnames):
            rows.append(dict(zip(fieldnames, pending)))
            starts.add(pending[0])        # values that begin a variant (e.g. the genome build)
            pending = []
    if pending:                           # trailing continuation or a truncated last line
        if rows and len(pending) < len(fieldnames) and pending[0] not in starts:
            rows[-1][fieldnames[-1]] += " " + "\t".join(pending)
        else:
            rows.append(dict(zip(fieldnames, pending + [""] * (len(fieldnames) - len(pending)))))
    return fieldnames, rows


def clean(value):
    """Line breaks and tabs inside a cell would break the MAF for downstream steps."""
    return str(value).replace("\r", " ").replace("\n", " ").replace("\t", " ") if value is not None else ""


def write_table(path, fieldnames, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows({k: clean(v) for k, v in r.items()} for r in rows)


def split(table, oncokb_transcripts, use_genomic_query):
    """Return (fieldnames, protein_rows, genomic_rows); each row carries its original index.
    The table is read tolerantly: OncoKB free text can carry line breaks into a cell."""
    fieldnames, rows = read_annotated(table)
    protein, genomic = [], []
    for i, row in enumerate(rows):
        row[ROW_ID] = str(i)
        row[ORIGINAL_END] = row.get("End_Position", "")
        if use_genomic_query(row.get("Hugo_Symbol", ""), row.get("NM_Transcript", ""), oncokb_transcripts):
            # SMART keeps VCF-style alleles (with the anchor base); OncoKB's genomic-change query
            # needs End_Position consistent with them. MAF-style insertion ends (start + 1) make
            # the batch request fail.
            try:
                row["End_Position"] = str(int(row["Start_Position"]) + len(row["Reference_Allele"]) - 1)
            except (KeyError, ValueError):
                pass
            genomic.append(row)
        else:
            protein.append(row)
    return fieldnames + [ROW_ID, ORIGINAL_END], protein, genomic


def merge(protein_maf, genomic_maf, genomic_failed=False):
    """Return (fieldnames, rows) in original order with ONCOKB_MAF_QUERY added. Raises
    ValueError if any input row is missing from the MafAnnotator outputs."""
    fieldnames, rows = [], []
    genomic_label = ("HGVSp_Short (selected transcript differs from OncoKB's; genomic-change query failed)"
                     if genomic_failed else "Genomic_Change (selected transcript differs from OncoKB's)")
    for path, label in ((protein_maf, "HGVSp_Short"), (genomic_maf, genomic_label)):
        if not path or not os.path.exists(path):
            continue
        names, part = read_annotated(path)
        fieldnames += [n for n in names if n not in fieldnames]
        for row in part:
            row[QUERY_COLUMN] = label
            if ORIGINAL_END in row:
                row["End_Position"] = row.pop(ORIGINAL_END)
        rows += part
    rows.sort(key=lambda r: int(r[ROW_ID]))
    if [int(r[ROW_ID]) for r in rows] != list(range(len(rows))):
        raise ValueError("rows missing or duplicated after OncoKB annotation")
    fieldnames = [n for n in fieldnames if n not in (ROW_ID, ORIGINAL_END)] + [QUERY_COLUMN]
    return fieldnames, rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)
    s = sub.add_parser("split")
    s.add_argument("--table", required=True)
    s.add_argument("--token", required=True)
    s.add_argument("--protein", required=True)
    s.add_argument("--genomic", required=True)
    m = sub.add_parser("merge")
    m.add_argument("--protein", required=True)
    m.add_argument("--genomic")
    m.add_argument("--genomic-failed", action="store_true",
                   help="the genomic table was annotated by protein change after the genomic query failed")
    m.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    if args.command == "split":
        oncokb = _oncokb_module()
        transcripts = oncokb.load_oncokb_transcripts(args.token)
        fieldnames, protein, genomic = split(args.table, transcripts, oncokb.use_genomic_query)
        write_table(args.protein, fieldnames, protein)
        write_table(args.genomic, fieldnames, genomic)
        genes = sorted({r.get("Hugo_Symbol", "") for r in genomic})
        print(f"OncoKB MafAnnotator queries: {len(protein)} by protein change, {len(genomic)} by genomic change"
              + (f" ({', '.join(genes)})" if genes else ""))
    else:
        fieldnames, rows = merge(args.protein, args.genomic, args.genomic_failed)
        write_table(args.out, fieldnames, rows)


if __name__ == "__main__":
    sys.exit(main())
