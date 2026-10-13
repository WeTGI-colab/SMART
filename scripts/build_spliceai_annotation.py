#!/usr/bin/env python3
"""
Build a SpliceAI gene-annotation file restricted to a chosen set of transcripts.

Why
---
SpliceAI's published precomputed tables carry no transcript: one score per gene, computed on
GENCODE v24 canonical (December 2015), within a +/-50 bp window. For 16 of our 519 panel
transcripts that is a different transcript from the one we report -- MET, HRAS, AKT1, FGFR1,
CDKN2A, CUX1, TP63, MITF, MUTYH and TCF3 among them -- and nothing in the output says which case
applies. See issue #2.

Running SpliceAI ourselves fixes that, because the model takes a gene annotation as input: give
it the transcripts we report and the scores describe those transcripts. This script builds that
annotation.

Input
-----
A GENCODE annotation already in SpliceAI's own format, as shipped with a SpliceAI-lookup
deployment (`gencode.v43.annotation.txt.gz`). Its NAME column packs six fields:

    GENE---ENSG---ENST---<canonical yes|no>---<biotype>---<RefSeq ids, comma separated>

so the RefSeq accessions we select on are already there and no GTF parsing, no Ensembl database
and no extra dependency is needed. The exon coordinates are taken across unchanged, which is the
point: they are the structures the published tool uses, so a score computed against them is
comparable with one from the hosted service.

Output
------
SpliceAI's annotation format, which is what `spliceai -A` reads:

    #NAME   CHROM   STRAND  TX_START    TX_END  EXON_START  EXON_END
    BRAF--NM_004333 chr7    -   140719327   140924929   140719327,...   140719465,...

NAME is written as `GENE--REFSEQ` to match the convention of the laboratory's existing
GRCh37 preferred-transcript file, so a score can be traced to the transcript it was computed on.

Matching
--------
By RefSeq accession, ignoring the version, because a panel list and an annotation release drift
in version while meaning the same transcript. Every transcript asked for and not found is
reported: the annotation cannot invent a structure, and a silent omission here would reproduce
the failure this work exists to remove.
"""
import argparse
import gzip
import os
import sys
from collections import defaultdict


def read_requested(path):
    """The transcript accessions to keep, as {accession without version: as written}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            t = line.strip()
            if not t or t.startswith("#"):
                continue
            out[t.split(".")[0]] = t
    return out


def parse_name(name):
    """GENE---ENSG---ENST---canonical---biotype---refseqs -> dict, or None."""
    p = name.split("---")
    if len(p) < 6:
        return None
    return {"gene": p[0], "ensg": p[1], "enst": p[2],
            "canonical": p[3] == "yes", "biotype": p[4],
            "refseqs": [r.strip() for r in p[5].split(",") if r.strip() and r.strip() != "."]}


def build(annotation, requested, out_fh, report=sys.stderr):
    """Write the annotation rows for `requested`. Returns (written, missing, ambiguous)."""
    found = defaultdict(list)
    header = None
    opener = gzip.open if annotation.endswith(".gz") else open
    with opener(annotation, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                header = line.rstrip("\n").split("\t")
                continue
            fields = line.rstrip("\n").split("\t")
            meta = parse_name(fields[0])
            if not meta:
                continue
            for rs in meta["refseqs"]:
                base = rs.split(".")[0]
                if base in requested:
                    found[base].append((meta, rs, fields))

    out_fh.write("#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n")
    written, ambiguous = 0, []
    for base in sorted(found):
        rows = found[base]
        if len(rows) > 1:
            # one accession mapping to several structures: prefer the canonical, and say so
            canonical = [r for r in rows if r[0]["canonical"]]
            chosen = canonical[0] if canonical else rows[0]
            ambiguous.append((base, [r[0]["enst"] for r in rows],
                              chosen[0]["enst"], bool(canonical)))
        else:
            chosen = rows[0]
        meta, rs, fields = chosen
        out_fh.write("\t".join([f"{meta['gene']}--{rs}", fields[1], fields[2],
                                fields[3], fields[4], fields[5], fields[6]]) + "\n")
        written += 1

    missing = sorted(set(requested) - set(found))
    print(f"transcripts requested : {len(requested)}", file=report)
    print(f"written               : {written}", file=report)
    print(f"not in the annotation : {len(missing)}", file=report)
    for m in missing:
        print(f"    MISSING {requested[m]}", file=report)
    if ambiguous:
        print(f"\nresolved to one structure : {len(ambiguous)}", file=report)
        for base, all_enst, chosen, was_canonical in ambiguous:
            why = "canonical" if was_canonical else "first, none canonical"
            print(f"    {requested[base]:18} {len(all_enst)} structures -> {chosen} ({why})",
                  file=report)
    return written, missing, ambiguous


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--annotation", required=True,
                   help="GENCODE annotation in SpliceAI format, e.g. gencode.v43.annotation.txt.gz")
    p.add_argument("--transcripts", required=True,
                   help="one RefSeq accession per line (e.g. reference_lists/TSO500_transcripts_list.txt)")
    p.add_argument("--out", required=True, help="output annotation file for 'spliceai -A'")
    a = p.parse_args()

    requested = read_requested(a.transcripts)
    if not requested:
        print("ERROR: no transcripts read", file=sys.stderr)
        return 1
    with open(a.out, "w") as fh:
        written, missing, _ = build(a.annotation, requested, fh)
    print(f"\nwritten to {a.out}", file=sys.stderr)
    return 0 if written else 1


if __name__ == "__main__":
    sys.exit(main())
