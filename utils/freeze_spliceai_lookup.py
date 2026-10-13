#!/usr/bin/env python3
"""
Freeze the Broad's SpliceAI-lookup answers for a small set of public variants.

Why frozen rather than queried by the tests
-------------------------------------------
The hosted service is the one thing that can say whether our splice annotation is right rather
than merely self-consistent, because it is an independent implementation of the same model. But
it must not be called on every test run: it is rate limited, it asks explicitly that batch work
be done locally, and it returned HTTP 500 on four of 93 rows during this work -- a live test
would fail for reasons that are not ours.

So the answers are recorded once, with the date and the endpoint, and the tests compare against
the record. That also makes a later difference meaningful: if the frozen values stop matching,
either our annotation changed or their Gencode release did, and both are worth knowing.

Only public variants are listed -- the service's own front-page example, well-known hotspots and
positions chosen to exercise one property each. Nothing patient-derived is sent to a third party.

Note that the service runs a modified SpliceAI (github.com/bw2/SpliceAI) which supports
multi-nucleotide variants and exposes raw REF/ALT site scores that the released package does not,
so agreement with a local run of the release is expected for substitutions and not for MNVs.

Usage
-----
    python3 utils/freeze_spliceai_lookup.py tests/verification5/data/spliceai_lookup_frozen.json

Run it deliberately, not in CI. Queries are paced at 8 s.
"""
import json
import sys
import time
import urllib.request
from datetime import date

BASE = "https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/"
PACE_SECONDS = 8

# Each entry states the property it is there to exercise, so a later reader can tell whether a
# changed value matters or is incidental.
WANT = [
    ("16-2317763-T-A",  "ABCA3",
     "strongest signal at 106 bp, outside the precomputed +/-50 bp window"),
    ("7-140753336-A-T", "BRAF",   "V600E; coding hotspot, no splice signal expected"),
    ("9-21974826-A-G",  "CDKN2A", "two reported transcripts, p16INK4a and p14ARF"),
    ("7-116771900-T-C", "MET",    "our reported transcript is not the gene-level one"),
    ("6-31740453-G-T",  "MSH5",   "one position, three genes, nine transcripts"),
    ("2-208248388-C-T", "IDH1",   "control: agrees with the precomputed table"),
]

# DS_*_REF and DS_*_ALT are the raw site scores, which a delta cannot replace: they say whether
# a splice site was already there. That is the evidence issue #21 needs.
KEEP = ("t_id", "t_priority", "t_refseq_ids", "g_name", "t_strand",
        "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL",
        "DS_AG_REF", "DS_AG_ALT", "DS_AL_REF", "DS_AL_ALT")


def fetch(variant):
    url = f"{BASE}?hg=38&variant=chr{variant}&content-type=application/json"
    with urllib.request.urlopen(url, timeout=90) as r:
        return json.loads(r.read().decode())


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    dest = sys.argv[1]
    out = {"_provenance": {
        "source": "https://spliceailookup.broadinstitute.org/ (hosted SpliceAI-lookup API)",
        "endpoint": BASE,
        "queried": str(date.today()),
        "note": ("Frozen so the tests do not call the service on every run: it is rate limited "
                 "and asks that batch work be done locally. Refresh with "
                 "utils/freeze_spliceai_lookup.py. The service runs a modified SpliceAI "
                 "(github.com/bw2/SpliceAI) on its own Gencode release, so a difference here may "
                 "be a change at their end rather than ours."),
    }, "variants": {}}

    for variant, gene, why in WANT:
        try:
            payload = fetch(variant)
        except Exception as e:
            print(f"  {variant:18} FAILED {type(e).__name__}", file=sys.stderr)
            continue
        scores = [{k: s.get(k) for k in KEEP} for s in payload.get("scores", [])]
        out["variants"][variant] = {"gene": gene, "why": why,
                                    "distance": payload.get("distance"),
                                    "mask": payload.get("mask"), "scores": scores}
        print(f"  {variant:18} {gene:8} {len(scores)} transcript record(s)")
        time.sleep(PACE_SECONDS)

    if not out["variants"]:
        print("nothing fetched; the fixture was not written", file=sys.stderr)
        return 1
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(f"\nwritten: {dest}  ({len(out['variants'])} variants)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
