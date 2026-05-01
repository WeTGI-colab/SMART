#!/usr/bin/env python3
"""
get_oncokb_transcripts.py — Fetch OncoKB's canonical transcript list.

Calls /utils/allCuratedGenes and writes a transcript file (one NM_ accession
per line) containing the grch38RefSeq transcript that OncoKB uses internally
for each curated gene.  The output can be passed directly to the SMART pipeline
via --transcripts-file.

Usage:
    python3 scripts/get_oncokb_transcripts.py --token TOKEN [options]

Options:
    --token   TOKEN     OncoKB API token (required)
    --output  FILE      Output transcript file (default: oncokb_transcripts.txt)
    --tsv     FILE      Also write a TSV summary with gene + transcript columns
    --no-version        Strip version suffix from NM accessions (e.g. NM_005228
                        instead of NM_005228.5).  Default: keep version when
                        OncoKB provides one.
"""

import argparse
import json
import sys
import requests

ONCOKB_BASE = "https://www.oncokb.org/api/v1"
ENDPOINT     = f"{ONCOKB_BASE}/utils/allCuratedGenes"


def fetch_curated_genes(token: str) -> list[dict]:
    headers = {
        "Authorization": f"Bearer {token}",
        "accept": "application/json",
    }
    resp = requests.get(ENDPOINT, headers=headers, timeout=30)
    if resp.status_code == 401:
        sys.exit("ERROR: OncoKB returned 401 Unauthorized — check your token.")
    if not resp.ok:
        sys.exit(f"ERROR: OncoKB returned {resp.status_code}: {resp.text[:200]}")
    return resp.json()


def main():
    parser = argparse.ArgumentParser(
        description="Fetch OncoKB canonical transcripts (grch38RefSeq) for all curated genes."
    )
    parser.add_argument("--token",      required=True, help="OncoKB API token")
    parser.add_argument("--output",     default="oncokb_transcripts.txt",
                        help="Output transcript file (default: oncokb_transcripts.txt)")
    parser.add_argument("--tsv",        default=None,
                        help="Optional TSV summary file (gene + transcript columns)")
    parser.add_argument("--no-version", action="store_true",
                        help="Strip version suffix from NM accessions")
    args = parser.parse_args()

    print(f"Fetching curated gene list from OncoKB...")
    genes = fetch_curated_genes(args.token)
    print(f"  {len(genes)} curated genes returned.")

    # Collect results
    rows = []           # (hugo, nm, enst) tuples for TSV
    transcripts = []    # NM_ accessions for the transcript file
    no_refseq = []      # genes with no grch38RefSeq

    for g in sorted(genes, key=lambda x: x.get("hugoSymbol", "")):
        hugo  = g.get("hugoSymbol", "")
        nm    = g.get("grch38RefSeq", "") or ""
        enst  = g.get("grch38Isoform", "") or ""

        if not nm:
            no_refseq.append(hugo)
            rows.append((hugo, "", enst))
            continue

        if args.no_version:
            nm = nm.split(".")[0]

        rows.append((hugo, nm, enst))
        transcripts.append(nm)

    # Write transcript file
    with open(args.output, "w") as fh:
        fh.write("\n".join(transcripts) + "\n")
    print(f"  Wrote {len(transcripts)} transcripts → {args.output}")

    # Write optional TSV summary
    if args.tsv:
        with open(args.tsv, "w") as fh:
            fh.write("HugoSymbol\tgrch38RefSeq\tgrch38Isoform\n")
            for hugo, nm, enst in rows:
                fh.write(f"{hugo}\t{nm}\t{enst}\n")
        print(f"  Wrote TSV summary → {args.tsv}")

    # Report genes without a RefSeq transcript
    if no_refseq:
        print(f"\n  {len(no_refseq)} genes have no grch38RefSeq in OncoKB:")
        for h in no_refseq:
            print(f"    {h}")

    print("\nDone.")


if __name__ == "__main__":
    main()
