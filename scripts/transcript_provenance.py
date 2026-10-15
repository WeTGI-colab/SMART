#!/usr/bin/env python3
"""
Answer one question per column: did this value come from the transcript we reported?

Why
---
A user supplies `--transcripts-file` and reasonably expects every annotation in the output to
describe that transcript. For many columns it does. For some it cannot, because the source is
keyed on genomic position and knows nothing about transcripts. And for some it does not,
although the source *is* transcript-specific -- which is the case nobody can see from the
output and the one this report exists to surface.

Three categories, and the third is the point:

  resolved      the value demonstrably describes the reported transcript. Either the value
                carries the transcript itself (HGVSc, HGVSp, SpliceAI's GENE--REFSEQ symbol)
                or the annotation step matched on it (LOEUF runs with match_by=transcript).

  positional    the source has no transcript to offer: ClinVar, CIViC, COSMIC, Cancer
                Hotspots, CADD and gnomAD are keyed on position and allele. The value is not
                wrong and is not a gap; it simply is not transcript-specific. Said explicitly
                so that its absence from the "resolved" list is not read as a defect.

  unresolved    the source holds one value per transcript, the output keeps them all, and
                nothing marks which element is ours. dbNSFP is the case: it returns a
                '&'-separated list parallel to `Ensembl_transcriptid`, and the plugin is run
                without `transcript_match=1`, so for TP53 the reported ENST00000269305 sits at
                index 23 of 28 while the first value in the cell belongs to ENST00000413465 --
                a different score. Anyone reading the cell as a single value, or taking the
                first element, reads another transcript's prediction.

                Also unresolved: REVEL, run with `no_match=1` so a score is returned even when
                the transcript differs, and the output records no transcript to check against;
                and SpliceAI in `--spliceai lookup` mode, whose tables carry no transcript at
                all (GENCODE v24 canonical, issue #2).

The per-transcript columns are found by structure rather than by a list: a column whose element
count equals `Ensembl_transcriptid`'s, on every row that has such a list, is parallel to it.
That way a new dbNSFP column is covered the day it appears instead of the day someone
remembers to add it.
"""
import argparse
import csv
import json
import os
import sys

SENTINELS = {"", ".", "-", "nan", "NaN", "None", "NA", "[]"}
TX_LIST_COL = "Ensembl_transcriptid"

# Sources keyed on genomic position and allele: there is no transcript to honour. Prefix-based
# so a new column from the same source is covered.
POSITIONAL_PREFIXES = ("ClinVar_", "CIViC_", "COSMIC_", "CancerHotspots_", "CADD_",
                       "gnomAD", "GENIE_")
POSITIONAL_EXACT = {"CADD_PHRED", "CADD_RAW", "Existing_variation", "DOMAINS", "PUBMED",
                    "CLIN_SIG", "SOMATIC", "PHENO"}

# Values that carry the transcript inside them, so the claim is self-evidencing.
CARRIES_TRANSCRIPT = {"HGVSc": "ENST", "HGVSp": "ENSP", "SpliceAI_pred_SYMBOL": "--"}

# Transcript-specific in the source but not resolved to ours by the way the step is run.
UNRESOLVED_REASONS = {
    "REVEL": "the REVEL plugin runs with no_match=1, so a score is returned even when the "
             "transcript differs, and the output records no transcript to check it against",
}


def read_table(path):
    csv.field_size_limit(10 ** 7)
    with open(path, newline="") as fh:
        lines = fh.readlines()
    body = [l for l in lines if not l.startswith("#") and l.strip()]
    rows = [r for r in csv.reader(body, delimiter="\t") if r and any(f.strip() for f in r)]
    if not rows:
        raise SystemExit(f"{path}: no data rows")
    return rows[0], rows[1:]


def cell(row, header_index, col):
    i = header_index.get(col)
    if i is None or i >= len(row):
        return ""
    return (row[i] or "").strip()


def per_transcript_columns(header, rows, idx):
    """Columns that hold one value per transcript, found by being parallel to the ID list.

    Required on *every* row that carries a transcript list, not merely one: several
    unrelated columns match the element count by coincidence on a single row.
    """
    counts, total = {}, 0
    for row in rows:
        tx = cell(row, idx, TX_LIST_COL)
        if not tx or tx in SENTINELS or "&" not in tx:
            continue
        n = len(tx.split("&"))
        total += 1
        for col in header:
            v = cell(row, idx, col)
            if v and v not in SENTINELS and len(v.split("&")) == n:
                counts[col] = counts.get(col, 0) + 1
    return sorted(c for c, k in counts.items() if k == total and total), total


def whitelist_adherence(header, rows, idx, transcripts):
    """Rows whose selected transcript is not one the user supplied.

    The root guarantee, and the one with a silent failure path. `vcf2table.py` walks the CSQ
    blocks looking for a transcript in the whitelist and, finding none, falls back to the
    first block:

        if not nm_selected and first_vep is not None:
            # No match found -- use the first CSQ record

    So a variant in a gene the list does not cover, or one where VEP did not annotate the
    listed transcript, is reported on an arbitrary transcript with nothing saying so. Every
    other guarantee in this report is downstream of this one: a predictor faithfully resolved
    to the selected transcript is still the wrong answer if the selected transcript was not
    the one asked for.

    Returns (n_ok, [(gene, selected, feature), ...]).
    """
    want = {t.split(".")[0] for t in transcripts if t}
    ok, bad = 0, []
    for row in rows:
        sel = cell(row, idx, "NM_Transcript")
        if sel and sel.split(".")[0] in want:
            ok += 1
        else:
            bad.append((cell(row, idx, "Hugo_Symbol"), sel or "(none)",
                        cell(row, idx, "Feature")))
    return ok, bad


def unresolved_remnants(header, rows, idx):
    """Columns still holding '&' lists while the transcript list itself is resolved.

    Once dbNSFP runs with `transcript_match=1`, `Ensembl_transcriptid` holds the single
    matched transcript and the predictors hold single values -- so the parallelism that
    `per_transcript_columns` relies on no longer exists, and that function correctly finds
    nothing. A column that is *still* a list in that state has not been resolved by the
    flag, which is the remaining case worth reporting: dbNSFP keys MutationTaster on a
    different transcript set, so the plugin's documentation states the flag has no effect
    on those fields.

    Returns {column: number of rows still holding a list}.
    """
    out = {}
    for row in rows:
        tx = cell(row, idx, TX_LIST_COL)
        if not tx or tx in SENTINELS or "&" in tx:
            continue                       # unresolved or absent: not this case
        for col in header:
            # Only columns that could have been per-transcript. PHENO, PUBMED, SOMATIC,
            # DOMAINS and the OncoKB entries are lists by nature -- several phenotypes, several
            # citations -- not one value per transcript, so a list there is not a remnant of
            # anything and reporting it would bury the two columns that are.
            if col == TX_LIST_COL or col in POSITIONAL_EXACT:
                continue
            if col.startswith(POSITIONAL_PREFIXES) or col.startswith(("ONCOKB_", "TX_")):
                continue
            v = cell(row, idx, col)
            if v and v not in SENTINELS and "&" in v:
                out[col] = out.get(col, 0) + 1
    return out


def transcript_list_agrees_with_feature(header, rows, idx):
    """Rows where the resolved dbNSFP transcript is the one the row reports.

    With `transcript_match=1` this is self-evidencing: `Ensembl_transcriptid` holds the
    transcript dbNSFP matched, `Feature` holds the one VEP selected from the user's list, and
    if they agree then the predictors in that row describe the transcript the user asked for.
    Returns (agreeing, disagreeing) counts.
    """
    ok = bad = 0
    for row in rows:
        tx = cell(row, idx, TX_LIST_COL)
        feat = cell(row, idx, "Feature")
        if not tx or tx in SENTINELS or "&" in tx or not feat:
            continue
        if tx.split(".")[0] == feat.split(".")[0]:
            ok += 1
        else:
            bad += 1
    return ok, bad


def classify(col, per_tx):
    if col in CARRIES_TRANSCRIPT:
        return "resolved"
    if col == "LOEUF":
        return "resolved"
    if col in UNRESOLVED_REASONS:
        return "unresolved"
    if col in per_tx:
        return "unresolved"
    if col in POSITIONAL_EXACT or col.startswith(POSITIONAL_PREFIXES):
        return "positional"
    return "other"


def divergences(header, rows, idx, per_tx):
    """For each row, where our transcript's value differs from the cell's first element.

    The first element is what a reader takes when the cell looks like a single value, so the
    difference between it and ours is the size of the mistake available to be made.
    """
    out = []
    for row in rows:
        tx = cell(row, idx, TX_LIST_COL)
        ours = cell(row, idx, "Feature")
        if not tx or "&" not in tx or not ours:
            continue
        txs = tx.split("&")
        if ours not in txs:
            out.append(dict(gene=cell(row, idx, "Hugo_Symbol"), transcript=ours,
                            column=TX_LIST_COL, index=None, n=len(txs),
                            ours=None, first=txs[0], note="reported transcript absent"))
            continue
        at = txs.index(ours)
        for col in per_tx:
            if col in (TX_LIST_COL, "Ensembl_geneid", "Ensembl_proteinid", "genename"):
                continue
            v = cell(row, idx, col)
            if not v or v in SENTINELS:
                continue
            parts = v.split("&")
            if len(parts) != len(txs) or at >= len(parts):
                continue
            mine, first = parts[at].strip(), parts[0].strip()
            if mine != first:
                out.append(dict(gene=cell(row, idx, "Hugo_Symbol"), transcript=ours,
                                column=col, index=at, n=len(txs),
                                ours=mine, first=first, note=""))
    return out


def report(header, rows, out=sys.stdout, limit=12, transcripts=()):
    idx = {c: i for i, c in enumerate(header)}
    if TX_LIST_COL not in idx:
        print(f"{TX_LIST_COL} absent: dbNSFP was not run, so there are no per-transcript "
              f"lists to resolve.", file=out)
        per_tx, n_tx_rows = [], 0
    else:
        per_tx, n_tx_rows = per_transcript_columns(header, rows, idx)

    populated = [c for c in header
                 if any(cell(r, idx, c) not in SENTINELS for r in rows)]
    groups = {}
    for col in populated:
        groups.setdefault(classify(col, per_tx), []).append(col)

    print(f"Transcript provenance: {len(rows)} rows, {len(populated)} populated columns",
          file=out)
    if transcripts:
        ok, bad = whitelist_adherence(header, rows, idx, transcripts)
        print(f"  selected transcript is one the user supplied: {ok}/{len(rows)}", file=out)
        if bad:
            print(f"  *** {len(bad)} row(s) were annotated on a transcript that was NOT "
                  f"supplied -- vcf2table fell back to the first CSQ block:", file=out)
            for gene, sel, feat in bad[:10]:
                print(f"        {gene:10} selected={sel:18} Feature={feat}", file=out)
            print(f"      Every other guarantee below is downstream of this one.", file=out)
    for kind in ("resolved", "positional", "unresolved", "other"):
        cols = groups.get(kind, [])
        print(f"  {kind:12} {len(cols)}", file=out)

    if groups.get("unresolved"):
        print(f"\n  NOT resolved to the reported transcript ({len(groups['unresolved'])} "
              f"columns):", file=out)
        for col in groups["unresolved"]:
            why = UNRESOLVED_REASONS.get(
                col, "dbNSFP returns one value per transcript and the plugin is run without "
                     "transcript_match=1, so the cell keeps them all and marks none")
            print(f"    - {col}: {why}" if col in UNRESOLVED_REASONS else f"    - {col}",
                  file=out)
        if per_tx:
            # Only true while dbNSFP is run without transcript_match=1. With the flag on,
            # these columns are resolved and the only unresolved one left is REVEL, whose
            # reason is printed per column above -- so this paragraph would be false.
            print(f"\n    Reason for the dbNSFP columns above: one value per transcript, "
                  f"parallel to {TX_LIST_COL}; the plugin is run without "
                  f"transcript_match=1, so the cell keeps every transcript's value and marks "
                  f"none as the reported one. Pass --dbnsfp-transcript-match 1 to resolve "
                  f"them (issue #29).", file=out)

    remnants = unresolved_remnants(header, rows, idx) if TX_LIST_COL in idx else {}
    if remnants:
        print(f"\n  Still a list although the transcript is resolved "
              f"({len(remnants)} columns):", file=out)
        for col, n in sorted(remnants.items()):
            print(f"    - {col} (on {n} rows)", file=out)
        if any(c.startswith("MutationTaster") for c in remnants):
            print(f"    dbNSFP keys MutationTaster on a different transcript set, so "
                  f"transcript_match=1 has no effect on those fields -- the plugin documents "
                  f"this. See issue #29.", file=out)
        others = sorted(c for c in remnants if not c.startswith("MutationTaster"))
        if others:
            print(f"    No documented reason for: {', '.join(others)} -- worth checking "
                  f"whether these are per-transcript or multi-valued by nature.", file=out)

    if TX_LIST_COL in idx:
        ok, bad = transcript_list_agrees_with_feature(header, rows, idx)
        if ok or bad:
            print(f"\n  dbNSFP matched the reported transcript on {ok} row(s)"
                  + (f", and a DIFFERENT one on {bad}" if bad else "")
                  + ".", file=out)
            if bad:
                print(f"    A disagreement means the predictors on those rows do not "
                      f"describe the transcript the user supplied.", file=out)

    div = divergences(header, rows, idx, per_tx) if per_tx else []
    if div:
        print(f"\n  Where the reported transcript's value differs from the cell's first "
              f"element ({len(div)} instances; the first element is what a reader takes):",
              file=out)
        for d in div[:limit]:
            if d["note"]:
                print(f"    {d['gene']:8} {d['transcript']:18} {d['note']}", file=out)
            else:
                print(f"    {d['gene']:8} {d['transcript']:18} {d['column']:24} "
                      f"ours[{d['index']}/{d['n']}]={d['ours']:<10} first={d['first']}",
                      file=out)
        if len(div) > limit:
            print(f"    ... and {len(div) - limit} more", file=out)
    return groups, div


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--maf", required=True)
    ap.add_argument("--transcripts", help="the --transcripts-file the run used; enables "
                                         "the check that every row was annotated on a "
                                         "transcript the user actually supplied")
    ap.add_argument("--json", help="write the full classification here")
    ap.add_argument("--fail-on-unresolved", action="store_true",
                    help="exit 1 if any column is unresolved (off by default: the current "
                         "pipeline has unresolved columns by design and they are reported)")
    a = ap.parse_args()
    header, rows = read_table(a.maf)
    txs = ()
    if a.transcripts:
        with open(a.transcripts) as fh:
            txs = tuple(l.strip() for l in fh if l.strip())
    groups, div = report(header, rows, transcripts=txs)
    if a.json:
        with open(a.json, "w") as fh:
            json.dump({"groups": groups, "divergences": div}, fh, indent=1, sort_keys=True)
        print(f"\nwritten: {a.json}", file=sys.stderr)
    if a.fail_on_unresolved and groups.get("unresolved"):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
