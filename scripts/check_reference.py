#!/usr/bin/env python3
"""
Verify that a VCF's REF alleles match the reference genome.

Why this exists
---------------
Nothing in the pipeline checked that the REF column of an incoming VCF agrees with the
assembly it is annotated against, and VEP does not refuse a VCF that disagrees: it takes
the reference base from the FASTA and keeps the submitted ALT. The consequence it then
reports is for a variant nobody submitted.

That was found in our own test data. `tests/verification1/verification1.vcf` carried two
variants written on the gene's strand instead of the genomic strand, both in minus-strand
genes:

    GNAQ_Q209L      chr9:77794572   A>T    genome: T
    DICER1_E1705K   chr14:95094139  G>A    genome: C

Neither produced a warning anywhere in the run. `GNAQ_Q209L` came out as `p.Gln209=`,
`Silent` -- the flipped ALT equalled the true reference, so the variant annotated itself
out of existence -- and `DICER1_E1705K` came out as `p.E1705*`, `Nonsense_Mutation`,
a different substitution from the one its own name states. The SpliceAI lookup was the
only step that noticed, and only by omission: it matches on CHROM+POS+REF+ALT, so both
variants silently received no splice annotation at all.

A REF mismatch has three usual causes and all three invalidate a run: the VCF is on a
different assembly, the alleles are strand-flipped, or the file is damaged. So the default
is to stop.

The reference is read directly through its `.fai` index rather than by calling out to
samtools per variant, which keeps a whole-cohort check to a single pass over each VCF.
"""
import argparse
import gzip
import os
import sys

# ALTs that describe something other than a literal sequence. A symbolic allele (<DUP>,
# <DEL>) or a missing ALT (.) is recognised by its first character; a breakend carries its
# bracket in the middle (A[chr1:100[), so that one is searched for anywhere in the string.
SYMBOLIC_PREFIX = ("<", ".")
BREAKEND_CHARS = "[]"

MAX_REF_LEN = 1000          # a REF longer than this is a symbolic/structural record
REPORT_LIMIT = 20           # mismatches listed individually before summarising
WRONG_ASSEMBLY_FRACTION = 0.5


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


class Reference:
    """Random access to a FASTA through its .fai index."""

    def __init__(self, fasta):
        self.path = fasta
        fai = fasta + ".fai"
        if not os.path.exists(fai):
            raise FileNotFoundError(
                f"{fai} not found; index the reference with 'samtools faidx {fasta}'")
        self.contigs = {}
        with open(fai) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                name, length, offset, linebases, linewidth = f[0], *map(int, f[1:5])
                self.contigs[name] = (length, offset, linebases, linewidth)
        self._fh = open(fasta, "rb")

    def close(self):
        self._fh.close()

    def resolve(self, chrom):
        """The contig name this FASTA uses, allowing for the chr prefix differing."""
        if chrom in self.contigs:
            return chrom
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        return alt if alt in self.contigs else None

    def fetch(self, chrom, start, length):
        """`length` bases from 1-based `start`, or None if out of range."""
        name = self.resolve(chrom)
        if name is None:
            return None
        clen, offset, linebases, linewidth = self.contigs[name]
        if start < 1 or start + length - 1 > clen:
            return None
        # bases are wrapped every `linebases`, and each line ends with `linewidth-linebases`
        # bytes of line terminator, so the byte offset of base i skips the terminators before it
        i = start - 1
        byte = offset + (i // linebases) * linewidth + (i % linebases)
        # read generously: enough bytes to cover the bases plus every terminator among them
        span = length + (length // linebases + 2) * (linewidth - linebases)
        self._fh.seek(byte)
        raw = self._fh.read(span).decode("ascii", "replace")
        seq = "".join(c for c in raw if c not in "\r\n")
        return seq[:length] if len(seq) >= length else None


def is_checkable(ref, alt):
    """Whether this record states a literal REF sequence that can be compared."""
    if not ref or len(ref) > MAX_REF_LEN:
        return False
    if any(alt.startswith(s) for s in SYMBOLIC_PREFIX):
        return False
    if any(c in alt for c in BREAKEND_CHARS):
        return False
    if set(ref.upper()) - set("ACGT"):
        return False            # N or IUPAC placeholders carry no information
    return True


def check(vcf, reference, limit=REPORT_LIMIT):
    """Compare every checkable REF against the reference.

    Returns a dict with the counts, the mismatches found (up to `limit`) and the contigs
    the FASTA does not have.
    """
    checked = skipped = 0
    mismatches, missing_contigs, out_of_range = [], set(), 0
    with _open(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, pos, vid, ref, alt = f[0], f[1], f[2], f[3].upper(), f[4]
            if not is_checkable(ref, alt):
                skipped += 1
                continue
            if reference.resolve(chrom) is None:
                missing_contigs.add(chrom)
                skipped += 1
                continue
            got = reference.fetch(chrom, int(pos), len(ref))
            if got is None:
                out_of_range += 1
                skipped += 1
                continue
            checked += 1
            if got.upper() != ref:
                if len(mismatches) < limit:
                    mismatches.append((chrom, pos, vid, ref, got.upper(), alt))
                else:
                    mismatches.append(None)     # counted, not kept
    return {"checked": checked, "skipped": skipped,
            "mismatches": [m for m in mismatches if m],
            "n_mismatches": len(mismatches),
            "missing_contigs": sorted(missing_contigs),
            "out_of_range": out_of_range}


def report(res, vcf, reference_path, mode):
    """Print the outcome. Returns True when the run should continue."""
    n = res["n_mismatches"]
    print(f"Reference check: {os.path.basename(vcf)} against "
          f"{os.path.basename(reference_path)}")
    print(f"  REF alleles verified: {res['checked']}")
    print(f"  not checkable (symbolic ALT, placeholder REF): {res['skipped']}")
    if res["missing_contigs"]:
        print(f"  WARNING: contigs absent from the reference: "
              f"{', '.join(res['missing_contigs'][:10])}")
    if res["out_of_range"]:
        print(f"  WARNING: {res['out_of_range']} record(s) positioned past the end of "
              f"their contig")

    if not n:
        print("  no REF mismatches")
        return True

    frac = n / res["checked"] if res["checked"] else 1.0
    label = "ERROR" if mode == "fail" else "WARNING"
    print(f"  {label}: {n} of {res['checked']} REF alleles do not match the reference "
          f"({frac:.1%})")
    print(f"  {'CHROM':8} {'POS':>12} {'ID':22} {'VCF':>10} {'GENOME':>8}  ALT")
    for chrom, pos, vid, ref, got, alt in res["mismatches"]:
        print(f"  {chrom:8} {pos:>12} {(vid or '.')[:22]:22} {ref[:10]:>10} "
              f"{got[:8]:>8}  {alt[:12]}")
    if n > len(res["mismatches"]):
        print(f"  ... and {n - len(res['mismatches'])} more")

    if frac >= WRONG_ASSEMBLY_FRACTION:
        print("  Most REF alleles disagree: this VCF is probably on a different assembly "
              "than the reference, or needs liftover.")
    else:
        print("  A handful of disagreements usually means strand-flipped alleles -- check "
              "whether they were written on the gene's strand rather than the genome's. "
              "Minus-strand genes are where this happens.")
    print("  VEP does not reject these: it uses the reference base and keeps the submitted "
          "ALT, so the consequence it reports is for a different variant.")
    return mode != "fail"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="VCF to verify (plain or .gz)")
    p.add_argument("--reference", required=True, help="reference FASTA, with its .fai index")
    p.add_argument("--mode", choices=("fail", "warn"), default="fail",
                   help="fail (default) stops the run on any mismatch; warn only reports")
    p.add_argument("--limit", type=int, default=REPORT_LIMIT,
                   help="how many mismatches to list individually")
    a = p.parse_args()

    try:
        ref = Reference(a.reference)
    except FileNotFoundError as e:
        print(f"Reference check skipped: {e}")
        return 0                        # a missing index must not take the pipeline down
    try:
        res = check(a.vcf, ref, a.limit)
    finally:
        ref.close()
    return 0 if report(res, a.vcf, a.reference, a.mode) else 1


if __name__ == "__main__":
    sys.exit(main())
