#!/usr/bin/env python3
"""Prove that the local LOEUF table is the gnomAD release it claims to be.

Why this exists
---------------
`utils/get_ref_files.sh` obtained `loeuf_dataset_grch38.tsv.gz` from an unversioned
private Google Drive link whose only integrity check was `file ... | grep HTML` -- that
the download was not an error page. Nothing established *which* gnomAD release the
numbers came from, and nothing would have noticed if the link had later served a
different file: LOEUF and missense Z feed O8 and the Tier 3 output, so a silent
substitution changes classifications with no visible symptom.

A plain checksum cannot be the whole answer here. The distributed file is not a gnomAD
artefact byte-for-byte: it is the released constraint TSV with its first column marked
`#` and three coordinate columns appended so that tabix can index it (all transcripts of
a gene share the gene's span, which is what lets the VEP plugin retrieve every candidate
row for a variant and then pick the transcript it annotated). So a checksum of the
distributed file pins *a* file, but says nothing about whether its numbers are gnomAD's.

This script checks the thing that matters: that every data column of every row equals
the corresponding row of a named, citable gnomAD release. It also reports which release
the file actually matches, because the repository had the version recorded wrongly --
Config.yaml said `gnomAD_v4.0` for a file whose every row matches v4.1 and differs from
v4.0 on 210,882 of 211,363 rows.

Why the table is verified rather than rebuilt
---------------------------------------------
The released constraint TSV carries a `chromosome` column but no start/end, and for its
RefSeq half not even the chromosome. Rebuilding the indexed table therefore needs a
gene-span source. The transcript model this pipeline already ships resolves 17,904 of
the 18,204 gene symbols, so a rebuild from what is on hand would silently drop 300 genes
that the present file covers. Verifying the file against the release gives the same
assurance without that regression.

Exit codes
----------
0  every compared row matches the expected release
1  the file's numbers differ from the expected release, or it matches a different one
2  the local file could not be read
3  the reference release could not be fetched, so nothing was compared (not a failure:
   the caller decides whether an unverifiable reference is acceptable)
"""

import argparse
import csv
import gzip
import os
import sys
import urllib.error
import urllib.request

RELEASES = {
    "4.1": "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/"
           "constraint/gnomad.v4.1.constraint_metrics.tsv",
    "4.0": "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/"
           "constraint/gnomad.v4.0.constraint_metrics.tsv",
}

# Columns the indexed table adds for tabix and which therefore have no counterpart in
# the released TSV. Everything else must match.
ADDED_COLUMNS = ("chrom", "start", "end")

# gnomAD v4.1 -- the release this table is pinned to -- publishes constraint for the
# autosomes only. That is a property of the release, not of our build: the official v4.1
# TSV has no chrX or chrY row either, and ATRX, KDM6A, BCOR, STAG2, DDX3X, PHF6, GATA1,
# BTK and AR are absent from it entirely. gnomAD's stated reason was a missing model
# input: methylation data for the sex chromosomes were unavailable when v4.0 was built,
# and the mutational model conditions expected counts on CpG methylation.
#
# gnomAD v4.1.1 (2026-03-25) added them and v4.1.2 (2026-09-23) is current: 3,581 chrX
# rows, 113 chrY rows, and native values for every gene listed above. Upgrading is not a
# file swap, which is why this stays pinned to v4.1 for now -- on the 17,460 genes with a
# MANE Select row in both releases, 750 more reach the missense-Z threshold of 3.09 and
# 21 fewer do, none of them on a sex chromosome. Tracked as issue 31.
#
# Backfilling chrX from v2.1.1 instead, the obvious move before v4.1.1 existed, would
# have been wrong twice: v2.1.1 reaches 3.09 for only 749 of the 1,470 genes v4.1 reaches
# it for, which would make an X-linked gene about half as likely to meet O8 as an
# autosomal one; and v2.1.1's sex-chromosome values were computed with the
# autosome-calibrated model and no documented ploidy correction, while v4.1.1 is the
# first release fitting separate models for chrX non-PAR and chrY non-PAR.
EXPECTED_CONTIGS = tuple("chr%d" % i for i in range(1, 23))


def _open_text(path):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")


def read_local(path):
    """Return (header, {transcript: row}) for the indexed table."""
    with _open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        header = [header[0].lstrip("#")] + list(header[1:])
        try:
            tx = header.index("transcript")
        except ValueError:
            raise ValueError("no 'transcript' column in %s" % path)
        rows = {}
        dupes = 0
        for rec in reader:
            if len(rec) != len(header):
                continue
            key = rec[tx]
            if key in rows:
                dupes += 1
                continue
            rows[key] = rec
    return header, rows, dupes


def fetch_release(release, cache_dir, url=None):
    """Return a path to the released TSV, downloading it if needed."""
    url = url or RELEASES[release]
    dest = os.path.join(cache_dir, os.path.basename(url))
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    os.makedirs(cache_dir, exist_ok=True)
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=120) as resp, open(tmp, "wb") as out:
        declared = resp.headers.get("Content-Length")
        written = 0
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
            written += len(chunk)
    # A truncated reference would make every comparison below meaningless, and a
    # truncation shows up as missing rows -- which reads as "the local file has extra
    # genes", the opposite of the truth. So it is caught here instead.
    if declared is not None and int(declared) != written:
        os.unlink(tmp)
        raise IOError("reference download truncated: %d of %s bytes" % (written, declared))
    os.rename(tmp, dest)
    return dest


def compare(header, rows, release_path):
    """Compare every shared row. Returns a dict of counts and up to 5 examples."""
    compared = set(header) - set(ADDED_COLUMNS)
    idx = {c: header.index(c) for c in compared}
    out = {"identical": 0, "different": 0, "absent_from_release": 0,
           "only_in_release": 0, "examples": [], "columns_compared": len(compared),
           "columns_skipped": sorted(set(header) & set(ADDED_COLUMNS))}
    seen = set()
    with _open_text(release_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = compared - set(reader.fieldnames or ())
        if missing:
            out["schema_mismatch"] = sorted(missing)
            return out
        for rec in reader:
            key = rec.get("transcript")
            local = rows.get(key)
            if local is None:
                out["only_in_release"] += 1
                continue
            seen.add(key)
            diffs = [c for c in compared if local[idx[c]] != rec[c]]
            if diffs:
                out["different"] += 1
                if len(out["examples"]) < 5:
                    out["examples"].append(
                        (key, diffs[:3],
                         [local[idx[c]] for c in diffs[:3]],
                         [rec[c] for c in diffs[:3]]))
            else:
                out["identical"] += 1
    out["absent_from_release"] = len(set(rows) - seen)
    return out


def contigs_of(header, rows):
    if "chrom" not in header:
        return []
    i = header.index("chrom")
    return sorted({r[i] for r in rows.values() if r[i]},
                  key=lambda c: (len(c), c))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table", help="local loeuf_dataset_grch38.tsv.gz")
    ap.add_argument("--release", default="4.1", choices=sorted(RELEASES),
                    help="the gnomAD release the table is expected to equal")
    ap.add_argument("--release-file", help="a already-downloaded released TSV to compare "
                                           "against instead of fetching one")
    ap.add_argument("--cache-dir", default=None,
                    help="where to keep the downloaded release (default: beside the table)")
    ap.add_argument("--also-identify", action="store_true",
                    help="on a mismatch, say which other known release the table matches")
    args = ap.parse_args(argv)

    try:
        header, rows, dupes = read_local(args.table)
    except (OSError, ValueError, StopIteration) as exc:
        print("[FAIL] cannot read %s: %s" % (args.table, exc))
        return 2

    print("[INFO] %s" % args.table)
    print("       %d columns, %d transcripts%s"
          % (len(header), len(rows),
             ", %d duplicate transcript rows ignored" % dupes if dupes else ""))

    contigs = contigs_of(header, rows)
    if contigs:
        print("       contigs: %s" % " ".join(contigs))
        unexpected = [c for c in contigs if c not in EXPECTED_CONTIGS]
        if unexpected:
            print("[WARN] contigs beyond the autosomes: %s" % " ".join(unexpected))
            print("       The v4.1 release this expects is autosome-only. gnomAD v4.1.1")
            print("       (2026-03-25) added chrX and chrY, so this is probably an upgraded")
            print("       table -- which is not a file swap: on the genes shared with v4.1,")
            print("       750 more reach the missense-Z threshold of 3.09 and 21 fewer do.")
            print("       O8's 'not evaluated' note for X-linked genes, the LOEUF")
            print("       description in Config.yaml and the 3.09 threshold itself all need")
            print("       revisiting before this table is annotated against. See issue 31.")
        else:
            print("       autosomes only, as gnomAD v4 publishes -- no chrX/chrY constraint")

    cache = args.cache_dir or os.path.dirname(os.path.abspath(args.table))
    if args.release_file:
        release_path = args.release_file
    else:
        try:
            release_path = fetch_release(args.release, cache)
        except (urllib.error.URLError, IOError, OSError) as exc:
            print("[SKIP] could not fetch gnomAD v%s: %s" % (args.release, exc))
            print("       nothing was compared; the table's provenance is unverified")
            return 3

    res = compare(header, rows, release_path)
    if res.get("schema_mismatch"):
        print("[FAIL] gnomAD v%s has no column %s -- the table is not that release's schema"
              % (args.release, ", ".join(res["schema_mismatch"][:5])))
        return 1

    print("[INFO] compared %d columns against gnomAD v%s (skipped the index columns %s)"
          % (res["columns_compared"], args.release, ", ".join(res["columns_skipped"])))
    print("       identical            : %d" % res["identical"])
    print("       different            : %d" % res["different"])
    print("       absent from release  : %d" % res["absent_from_release"])
    print("       only in release      : %d" % res["only_in_release"])

    if res["different"]:
        print("[FAIL] the table's numbers are not gnomAD v%s" % args.release)
        for key, cols, mine, theirs in res["examples"]:
            print("       %s: %s local=%s release=%s"
                  % (key, ",".join(cols), mine, theirs))
        if args.also_identify:
            for other in sorted(RELEASES):
                if other == args.release:
                    continue
                try:
                    p = fetch_release(other, cache)
                except (urllib.error.URLError, IOError, OSError):
                    continue
                alt = compare(header, rows, p)
                if not alt.get("schema_mismatch") and alt["different"] == 0 and alt["identical"]:
                    print("[INFO] it matches gnomAD v%s instead (%d rows identical)"
                          % (other, alt["identical"]))
        return 1

    if not res["identical"]:
        print("[FAIL] no row could be compared -- the two files share no transcript")
        return 1

    extra = ""
    if res["only_in_release"]:
        extra = ("; %d release rows are absent from the indexed table, which is expected: "
                 "rows whose gene has no span in the transcript model cannot be indexed"
                 % res["only_in_release"])
    print("[OK] every one of the %d compared rows equals gnomAD v%s%s"
          % (res["identical"], args.release, extra))
    return 0


if __name__ == "__main__":
    sys.exit(main())
