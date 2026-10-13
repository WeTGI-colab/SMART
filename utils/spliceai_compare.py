#!/usr/bin/env python3
"""
Score variants with a locally run SpliceAI and compare with the precomputed tables.

What this measures
------------------
SMART annotates splicing by looking variants up in Illumina's precomputed score files. Those
carry no transcript -- one score per gene, computed on GENCODE v24 canonical (December 2015) --
and were precomputed within a +/-50 bp window, verified by scanning the files (max |delta
position| = 50 over 199,971 records). Issue #2 has the detail.

Running the model locally removes both limits: the annotation is ours, so the score describes
the transcript we report, and the distance is a parameter. The open question is not whether that
is more correct but whether it changes anything that matters, and at what cost. This script
produces those numbers:

  - how many variants change score at all
  - how many cross 0.2, the delta from which the classifier awards O6 and from which
    ANNOTATION_FLAGS keeps a LOW/MODIFIER row
  - how many cross 0.5
  - where the strongest signal lies beyond 50 bp, which the precomputed tables cannot see
  - wall-clock per variant, which is the one thing that could make this impractical

How it works
------------
No model code is reimplemented. The `spliceai` CLI reads a VCF and writes a VCF, and `-A` takes
a gene annotation; `scripts/build_spliceai_annotation.py` builds ours from the panel transcript
list. Because SpliceAI reports the annotation's NAME column as its SYMBOL field, and our
annotation names rows `GENE--REFSEQ`, **the output states which transcript each score belongs
to** -- the thing the precomputed tables never recorded.

Usage
-----
    python3 utils/spliceai_compare.py \\
        --maf  path/to/Final_result_tier1.maf \\
        --annotation  .../tso500_grch38.spliceai_annotation.txt \\
        --reference   .../hg38.fa \\
        --precomputed .../SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz \\
        --distance 500 --out comparison.tsv

`--spliceai` points at the executable if it is not on PATH (a venv, typically).
"""
import argparse
import csv
import os
import re
import subprocess
import sys
import time

FIELDS = ["ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL",
          "DP_AG", "DP_AL", "DP_DG", "DP_DL"]
DS = ["DS_AG", "DS_AL", "DS_DG", "DS_DL"]
EMPTY = ("", ".", "-", "nan", "NaN", "None", "NA")

# the thresholds this is measured against: the classifier awards O6 from 0.2, and
# ANNOTATION_FLAGS keeps a LOW/MODIFIER row from the same value
THRESHOLDS = (0.2, 0.5)
PRECOMPUTED_WINDOW = 50


def read_maf(path):
    """Substitutions and indels from a MAF, without losing values containing '#'.

    `comment="#"` is wrong on a MAF: a comment character is honoured mid-line, so a row is
    truncated at the '#' inside a value such as ClinVar's UniProtKB:P01111#VAR_006847. See
    issue #4.
    """
    csv.field_size_limit(10 ** 7)
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if len(r) > 1]
    hdr = rows[0]
    need = ("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2")
    missing = [c for c in need if c not in hdr]
    if missing:
        raise SystemExit(f"MAF lacks {missing}")
    idx = {c: hdr.index(c) for c in hdr}
    out = []
    for r in rows[1:]:
        alt = r[idx["Tumor_Seq_Allele2"]]
        if alt.startswith("<") or not alt or alt in EMPTY:
            continue                      # symbolic alleles are not sequence
        out.append({"chrom": r[idx["Chromosome"]], "pos": r[idx["Start_Position"]],
                    "ref": r[idx["Reference_Allele"]], "alt": alt,
                    "gene": r[idx["Hugo_Symbol"]] if "Hugo_Symbol" in idx else "",
                    "feature": r[idx["Feature"]] if "Feature" in idx else ""})
    return out


def write_vcf(variants, path, reference):
    """A minimal VCF for the model. Contigs are declared from the reference's .fai."""
    contigs = []
    fai = reference + ".fai"
    if os.path.exists(fai):
        with open(fai) as fh:
            for line in fh:
                f = line.split("\t")
                contigs.append(f"##contig=<ID={f[0]},length={f[1]}>")
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("\n".join(contigs) + ("\n" if contigs else ""))
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in variants:
            fh.write(f"{v['chrom']}\t{v['pos']}\t.\t{v['ref']}\t{v['alt']}\t.\t.\t.\n")


def parse_spliceai_vcf(path):
    """{(chrom,pos,ref,alt): [record dicts]} from a SpliceAI-annotated VCF."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            key = (f[0], f[1], f[3], f[4])
            recs = []
            m = re.search(r"SpliceAI=([^;\t]+)", f[7])
            if m:
                for rec in m.group(1).split(","):
                    parts = rec.split("|")
                    if len(parts) == len(FIELDS):
                        recs.append(dict(zip(FIELDS, parts)))
            out[key] = recs
    return out


def tabix_precomputed(path, chrom, pos, ref, alt):
    """The precomputed record(s) for this variant, or []."""
    c = chrom[3:] if chrom.startswith("chr") else chrom
    r = subprocess.run(["tabix", path, f"{c}:{pos}-{pos}"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return []
    out = []
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) > 7 and f[3] == ref and alt in f[4].split(","):
            m = re.search(r"SpliceAI=([^;\t]+)", f[7])
            if m:
                for rec in m.group(1).split(","):
                    parts = rec.split("|")
                    if len(parts) == len(FIELDS):
                        out.append(dict(zip(FIELDS, parts)))
    return out


def fnum(v):
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return None


def max_delta(records):
    """(highest delta, its type, its position, the record it came from) or (None, ...)."""
    best = (None, None, None, None)
    for rec in records:
        for k in DS:
            x = fnum(rec.get(k))
            if x is not None and (best[0] is None or x > best[0]):
                best = (x, k, fnum(rec.get(k.replace("DS", "DP"))), rec)
    return best


def crossings(a, b):
    """Which thresholds the two values fall on opposite sides of."""
    out = []
    for t in THRESHOLDS:
        if (a is not None and a >= t) != (b is not None and b >= t):
            out.append(t)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--maf", required=True)
    p.add_argument("--annotation", required=True, help="from build_spliceai_annotation.py")
    p.add_argument("--reference", required=True, help="hg38.fa, with its .fai")
    p.add_argument("--precomputed", help="the precomputed SNV table, to compare against")
    p.add_argument("--distance", type=int, default=500)
    p.add_argument("--mask", type=int, default=0, choices=(0, 1))
    p.add_argument("--spliceai", default="spliceai", help="the executable")
    p.add_argument("--out", default="spliceai_comparison.tsv")
    p.add_argument("--workdir", default=None, help="where to put the intermediate VCFs")
    a = p.parse_args()

    variants = read_maf(a.maf)
    if not variants:
        raise SystemExit("no usable variants in the MAF")
    work = a.workdir or os.path.dirname(os.path.abspath(a.out)) or "."
    os.makedirs(work, exist_ok=True)
    in_vcf = os.path.join(work, "spliceai_input.vcf")
    out_vcf = os.path.join(work, "spliceai_output.vcf")
    write_vcf(variants, in_vcf, a.reference)
    print(f"{len(variants)} variants -> {in_vcf}")

    cmd = [a.spliceai, "-I", in_vcf, "-O", out_vcf, "-R", a.reference,
           "-A", a.annotation, "-D", str(a.distance), "-M", str(a.mask)]
    print("running: " + " ".join(cmd))
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0
    if r.returncode != 0:
        sys.stderr.write(r.stdout[-4000:] + "\n" + r.stderr[-4000:] + "\n")
        raise SystemExit(f"spliceai failed ({r.returncode})")
    print(f"scored in {elapsed:.1f}s — {elapsed / len(variants):.2f}s per variant")

    local = parse_spliceai_vcf(out_vcf)
    rows, changed, crossed, beyond_window, no_pre = [], 0, [], 0, 0
    for v in variants:
        key = (v["chrom"], v["pos"], v["ref"], v["alt"])
        lrecs = local.get(key, [])
        lmax, ltype, lpos, lrec = max_delta(lrecs)
        prec = tabix_precomputed(a.precomputed, **v_key(v)) if a.precomputed else []
        pmax, ptype, ppos, _ = max_delta(prec)
        if a.precomputed and not prec:
            no_pre += 1
        x = crossings(pmax, lmax)
        if x:
            crossed.append((v, pmax, lmax, x))
        if lmax is not None and pmax is not None and abs(lmax - pmax) > 1e-9:
            changed += 1
        if lpos is not None and abs(lpos) > PRECOMPUTED_WINDOW:
            beyond_window += 1
        rows.append({
            "Chromosome": v["chrom"], "Start_Position": v["pos"],
            "Reference_Allele": v["ref"], "Tumor_Seq_Allele2": v["alt"],
            "Hugo_Symbol": v["gene"], "Reported_Feature": v["feature"],
            "Local_max_DS": "" if lmax is None else f"{lmax:g}",
            "Local_type": ltype or "", "Local_DP": "" if lpos is None else lpos,
            "Local_transcript": (lrec or {}).get("SYMBOL", ""),
            "Local_n_transcripts": len(lrecs),
            "Precomputed_max_DS": "" if pmax is None else f"{pmax:g}",
            "Precomputed_type": ptype or "", "Precomputed_DP": "" if ppos is None else ppos,
            "Crosses": ";".join(str(t) for t in x),
        })

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"\n{'':4}variants scored locally        : "
          f"{sum(1 for r in rows if r['Local_max_DS'] != '')}/{len(rows)}")
    if a.precomputed:
        print(f"{'':4}absent from the precomputed table: {no_pre}")
        print(f"{'':4}score differs                  : {changed}")
        print(f"{'':4}crosses 0.2 or 0.5             : {len(crossed)}")
        for v, pm, lm, x in crossed[:20]:
            print(f"      {v['gene']:10} {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}  "
                  f"precomputed {pm} -> local {lm}   crosses {x}")
    print(f"{'':4}strongest signal beyond +/-{PRECOMPUTED_WINDOW} bp : {beyond_window}"
          f"   (invisible to the precomputed tables)")
    print(f"\nwritten: {a.out}")
    return 0


def v_key(v):
    return {"chrom": v["chrom"], "pos": v["pos"], "ref": v["ref"], "alt": v["alt"]}


if __name__ == "__main__":
    sys.exit(main())
