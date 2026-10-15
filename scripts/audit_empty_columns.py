#!/usr/bin/env python3
"""
Report every output column that never carries a value, and require each to be accounted for.

Why this exists
---------------
`ONCOKB_VUS` was empty on every row of every run for years. The cause was one character of
case: the code read `data.get("VUS")` where OncoKB's field is `vus`. Nothing caught it and
nothing could have, because the three states are indistinguishable in the output --

  * the annotation source has no value for this variant,
  * the column does not apply to this kind of variant,
  * the column is broken and will never carry a value,

all look like an empty cell. `dict.get` returns None for an absent key exactly as it does for
a null value, and the pipeline writes "." for both.

The only way to tell them apart is to ask a different question: *is there any variant, in any
run, for which this column carried a value?* A column that never does is either genuinely
inapplicable or silently broken, and someone has to say which. This script makes that
judgement explicit and committed, rather than implicit and forgotten.

How it is used
--------------
`empty_columns_allowed.tsv` lists one row per column (or glob), with a `kind` and a reason:

  structural   no variant could fill it as the pipeline runs today. Changing this means
               changing the pipeline, not the test data -- e.g. everything that needs a
               tumour-type query, which SMART does not make (issue #26).
  cohort       a variant would fill it; the test cohort simply has none. This is a to-do
               list: each row names the kind of variant that would cover it.

The audit then fails in both directions, and the second is the one that matters here:

  * an always-empty column with no allowlist row -- possibly the next ONCOKB_VUS;
  * an allowlist row that no longer matches anything -- the column now carries data, so the
    entry is stale and must go. That is what happens when a mutation is added to the cohort,
    which is the point: adding variants is what turns an unexamined blank into a covered one.
"""
import argparse
import csv
import fnmatch
import os
import sys

# Values that mean "nothing here". "." is MafAnnotator's, "-" is VEP's, the rest are what a
# str() of a missing value looks like by the time it reaches a TSV.
SENTINELS = {"", ".", "-", "nan", "NA", "NaN", "None", "null"}


def read_table(path):
    """(header, rows) from a MAF/TSV, skipping the leading '#' metadata lines.

    Deliberately not pandas' comment="#": it honours the marker mid-line too, so a row is
    truncated at the '#' inside a value such as ClinVar's UniProtKB:P01111#VAR_006847 and
    every later column comes back empty -- which this audit would then report as dozens of
    broken columns (issue #4).
    """
    csv.field_size_limit(10 ** 7)
    with open(path, newline="") as fh:
        lines = fh.readlines()
    # Drop the leading metadata lines by looking at the start of the line only. A previous
    # version relied on them having a single field, which worked by accident and discarded a
    # legitimate single-column table along with them. The marker is only a comment at the
    # start of a line: '#' appears inside values too, such as ClinVar's
    # UniProtKB:P01111#VAR_006847, and treating those as comments truncates the row and loses
    # every later column (issue #4).
    body = [l for l in lines if not l.startswith("#") and l.strip()]
    rows = [r for r in csv.reader(body, delimiter="\t") if r and any(f.strip() for f in r)]
    if not rows:
        raise SystemExit(f"{path}: no data rows")
    return rows[0], rows[1:]


def always_empty(header, rows):
    """The columns for which no row carries a non-sentinel value."""
    out = []
    for i, col in enumerate(header):
        seen = {(r[i] if i < len(r) else "").strip() for r in rows}
        if not (seen - SENTINELS):
            out.append(col)
    return out


def load_allowlist(path):
    """[(pattern, kind, reason)] from the TSV, ignoring comments and the header."""
    entries = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3 or parts[0] == "pattern":
                continue
            entries.append((parts[0].strip(), parts[1].strip(), parts[2].strip()))
    return entries


def audit(header, rows, allow):
    """(unexplained, stale, covered) -- see the module docstring for what each means."""
    empty = always_empty(header, rows)
    unexplained = [c for c in empty
                   if not any(fnmatch.fnmatchcase(c, p) for p, _, _ in allow)]
    stale = [(p, k, r) for p, k, r in allow
             if not any(fnmatch.fnmatchcase(c, p) for c in empty)
             and any(fnmatch.fnmatchcase(c, p) for c in header)]
    return unexplained, stale, empty


def report(header, rows, allow, out=sys.stdout):
    unexplained, stale, empty = audit(header, rows, allow)
    kinds = {}
    for col in empty:
        for p, k, r in allow:
            if fnmatch.fnmatchcase(col, p):
                kinds.setdefault(k, []).append((col, r))
                break
    print(f"Empty-column audit: {len(header)} columns, {len(rows)} rows", file=out)
    print(f"  never carry a value : {len(empty)}", file=out)
    for k in sorted(kinds):
        print(f"    {k:10} {len(kinds[k])}", file=out)
    if kinds.get("cohort"):
        print("\n  Coverable by adding a variant -- these are gaps in the test data, "
              "not in the code:", file=out)
        seen = set()
        for col, reason in kinds["cohort"]:
            if reason in seen:
                continue
            seen.add(reason)
            print(f"    - {reason}", file=out)
    if unexplained:
        print(f"\n  UNACCOUNTED ({len(unexplained)}) -- each is either inapplicable or the "
              f"next ONCOKB_VUS:", file=out)
        for c in unexplained:
            print(f"    ? {c}", file=out)
    if stale:
        print(f"\n  STALE ALLOWLIST ROWS ({len(stale)}) -- the column now carries data, so "
              f"remove the entry:", file=out)
        for p, k, r in stale:
            print(f"    x {p}  ({k}: {r})", file=out)
    ok = not unexplained and not stale
    print(f"\n{'OK' if ok else 'FAILED'}", file=out)
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--maf", required=True)
    # Repeatable, because the two MAFs audited do not have the same columns and are not
    # populated the same way. post_analysis pads the OncoKB indexed columns out to a fixed
    # count, so ONCOKB_TX_*_abstracts is wholly populated in FINAL_Table and partly empty
    # in the Tier 1 MAF: one list describing both would be stale for one table while being
    # needed by the other. The shared list carries what is true of any output; a second
    # list carries what is true only of the later table.
    ap.add_argument("--allowlist", required=True, action="append",
                    help="allowlist TSV; repeat to combine several")
    a = ap.parse_args()
    header, rows = read_table(a.maf)
    allow = [row for path in a.allowlist for row in load_allowlist(path)]
    return 0 if report(header, rows, allow) else 1


if __name__ == "__main__":
    sys.exit(main())
