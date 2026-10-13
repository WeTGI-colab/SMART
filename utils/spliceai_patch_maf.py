#!/usr/bin/env python3
"""
Replace a MAF's SpliceAI columns with locally computed, transcript-aware scores.

Purpose
-------
To find out whether computing SpliceAI changes anything that matters, the classifier has to be
run twice on the same variants: once with the precomputed scores the pipeline uses today, and
once with scores computed against the transcripts we actually report. This produces the second
MAF, so the two Tier 4 outputs differ in exactly one input.

Columns replaced: `SpliceAI_pred_DS_{AG,AL,DG,DL}`, `SpliceAI_pred_DP_{AG,AL,DG,DL}` and
`SpliceAI_pred_SYMBOL`. The symbol is set to the annotation's name for the row, `GENE--REFSEQ`,
so the patched MAF records which transcript each score was computed on -- the thing the
precomputed tables never carried.

`SpliceAI_cutoff` is deliberately left alone: it is the plugin's own PASS/FAIL verdict against
its threshold, not a value from the source, and recomputing it here would invent an annotation
the pipeline did not make.

Rows the local run did not score
--------------------------------
SpliceAI scores single-base substitutions whose REF matches the reference genome, so an indel or
a malformed record gets nothing. `--on-missing` decides what the patched MAF then says:

  blank (default)  empty, which is what a local-only pipeline would produce. Honest, and it
                   means a difference in the result can be attributed to the change.
  keep             the precomputed value, which mixes the two sources in one column and makes
                   any downstream difference ambiguous. Offered only for a sensitivity check.

Either way the count is reported, because a silent substitution is the failure this work exists
to remove.
"""
import argparse
import csv
import os
import re
import sys

DS = ["SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
      "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"]
DP = ["SpliceAI_pred_DP_AG", "SpliceAI_pred_DP_AL",
      "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL"]
SYMBOL = "SpliceAI_pred_SYMBOL"
ORDER = ["ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL",
         "DP_AG", "DP_AL", "DP_DG", "DP_DL"]


def read_local(path):
    """{(chrom,pos,ref,alt): record} from a SpliceAI-annotated VCF.

    Where several annotation rows overlap a variant the highest delta wins, which is the
    convention the pipeline's own columns follow. With a one-transcript-per-gene annotation
    this almost never arises; it is handled so that it cannot pass silently.
    """
    out, multi = {}, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            m = re.search(r"SpliceAI=([^;\t]+)", f[7])
            if not m:
                continue
            recs = []
            for rec in m.group(1).split(","):
                p = rec.split("|")
                if len(p) == len(ORDER):
                    recs.append(dict(zip(ORDER, p)))
            if not recs:
                continue
            if len(recs) > 1:
                multi += 1

            def peak(r):
                vals = []
                for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL"):
                    try:
                        vals.append(float(r[k]))
                    except (TypeError, ValueError):
                        pass
                return max(vals) if vals else -1.0

            out[(f[0], f[1], f[3], f[4])] = max(recs, key=peak)
    return out, multi


def patch(maf_in, local, maf_out, on_missing="blank", report=sys.stderr):
    csv.field_size_limit(10 ** 7)
    with open(maf_in, newline="") as fh:
        lines = fh.readlines()
    # keep the leading metadata lines (#SMART_VERSION ...) exactly as they are
    lead = 0
    while lead < len(lines) and lines[lead].startswith("#"):
        lead += 1
    hdr = lines[lead].rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(hdr)}
    need = ("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2")
    if any(c not in idx for c in need):
        raise SystemExit(f"MAF lacks {[c for c in need if c not in idx]}")
    targets = [c for c in DS + DP + [SYMBOL] if c in idx]
    if not targets:
        raise SystemExit("MAF has no SpliceAI columns to replace")

    patched = missing = 0
    out = lines[:lead + 1]
    for line in lines[lead + 1:]:
        f = line.rstrip("\n").split("\t")
        if len(f) < len(hdr):
            f += [""] * (len(hdr) - len(f))
        key = (f[idx["Chromosome"]], f[idx["Start_Position"]],
               f[idx["Reference_Allele"]], f[idx["Tumor_Seq_Allele2"]])
        rec = local.get(key)
        if rec:
            for c in DS + DP:
                if c in idx:
                    f[idx[c]] = rec[c.replace("SpliceAI_pred_", "")]
            if SYMBOL in idx:
                f[idx[SYMBOL]] = rec["SYMBOL"]
            patched += 1
        else:
            missing += 1
            if on_missing == "blank":
                for c in targets:
                    f[idx[c]] = ""
        out.append("\t".join(f) + "\n")

    with open(maf_out, "w") as fh:
        fh.writelines(out)
    total = patched + missing
    print(f"rows                        : {total}", file=report)
    print(f"  patched with local scores : {patched}", file=report)
    print(f"  not scored locally        : {missing}  ({on_missing})", file=report)
    return patched, missing


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--maf", required=True)
    p.add_argument("--spliceai-vcf", required=True, help="output of the local spliceai run")
    p.add_argument("--out", required=True)
    p.add_argument("--on-missing", choices=("blank", "keep"), default="blank")
    a = p.parse_args()
    local, multi = read_local(a.spliceai_vcf)
    print(f"local scores read           : {len(local)}", file=sys.stderr)
    if multi:
        print(f"  variants with >1 annotation row: {multi} (highest delta kept)",
              file=sys.stderr)
    patch(a.maf, local, a.out, a.on_missing)
    print(f"\nwritten: {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
