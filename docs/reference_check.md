# Checking a VCF's REF alleles against the reference genome

Before VEP runs, SMART verifies that the REF column of the VCF it is about to annotate
agrees with the reference genome. A mismatch stops the run.

## Why

Nothing checked this before, and no other step catches it. **VEP does not reject a VCF whose
REF disagrees with the FASTA:** it takes the reference base from the genome and keeps the
submitted ALT. The consequence it then reports belongs to a variant nobody submitted, and
nothing in the output says so.

This was found in SMART's own test data. `tests/verification1/verification1.vcf` contains
three records whose REF does not match GRCh38:

| ID | Position | VCF REF | Genome | What SMART reported |
|---|---|---|---|---|
| `GNAQ_Q209L` | chr9:77794572 | `A` | `T` | `p.Gln209=`, `Silent` |
| `DICER1_E1705K` | chr14:95094139 | `G` | `C` | `p.E1705*`, `Nonsense_Mutation` |
| `EGFR_del19` | chr7:55174772 | `GAATTAAGAGAAGCAT` | `GGAATTAAGAGAAGCA` | no SpliceAI annotation |

The first two are minus-strand genes whose alleles were written on the gene's strand rather
than the genome's. For `GNAQ_Q209L` the flipped ALT happened to equal the true reference, so
the variant annotated itself out of existence and a missense hotspot was reported as
synonymous. `DICER1_E1705K` became a different substitution that creates a stop codon, so a
missense was reported as a nonsense. The third is an off-by-one REF on a deletion.

None of them produced a warning anywhere in the run. The only step that noticed was the
SpliceAI lookup, and only by omission — it matches on CHROM+POS+REF+ALT, so these records
silently received no splice annotation at all. That is what prompted this check.

A REF mismatch has three usual causes, and all three invalidate a run:

- the VCF is on a different assembly than the reference (needs liftover);
- the alleles are strand-flipped, which is how the two above arose;
- the file is truncated or damaged.

## What it does

`scripts/check_reference.py` reads the reference through its `.fai` index — one seek per
record rather than a subprocess per variant — and compares every REF it can. It reports:

- how many REF alleles were verified;
- how many records were not checkable;
- every mismatch, with the VCF's REF and the genome's, up to a limit;
- contigs the reference does not contain, and positions past the end of a contig.

When most REF alleles disagree it says the assembly is probably wrong. When a handful
disagree it says to look for strand-flipped alleles, because that is the likelier cause.

Records that state no literal sequence are skipped rather than failed: a symbolic ALT
(`<DUP>`, `<DEL>`), a breakend (`A[chr1:100[`), a missing ALT (`.`), or a REF that is a
placeholder (`N`, or any IUPAC ambiguity code). Their REF carries no information to compare.

## Where it runs

Inside the per-sample loop, on `$VEP_INPUT` — the file VEP is about to read — and after
liftover rather than before it. That matters: with liftover enabled the input VCF is on
GRCh37 and would mismatch GRCh38 by design, so checking the original would fail every such
run. Checking what VEP receives verifies the thing that is actually annotated, in the
assembly the reference describes.

If the reference has no `.fai` index the check prints why and the run continues. A missing
index is a deployment problem and should not stop a pipeline that would otherwise work;
index it with `samtools faidx`.

## Controlling it

```
--ref-check fail    stop the run on any mismatch (default)
--ref-check warn    report and carry on
```

The default stops, because a single wrong allele produces a wrong consequence and therefore
a wrong classification, and in a diagnostic context that should not pass quietly. `warn`
exists for the case where the mismatches are known and understood and the run is wanted
anyway.

## Tests

`tests/verification5/test_check_reference.py` (42 tests): FASTA reading through a
hand-built `.fai` including line-wrapped spans, both chromosome naming conventions,
soft-masked lower-case reference, which records are checkable, detection of a strand flip,
an off-by-one multi-base REF and a wrong assembly, gzipped input, counting beyond the report
limit, and the exit status for each mode. Two tests assert the three bad records are still
present in `verification1.vcf` and record the correct alleles for whoever fixes it.
