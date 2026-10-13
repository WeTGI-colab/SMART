# Computing SpliceAI locally

Exploratory work on branch `spliceai-local-compute`, for [issue #2](https://github.com/WeTGI-colab/SMART/issues/2).

## Why

SMART annotates splicing by looking variants up in Illumina's precomputed score files. Those
have two limits that no configuration can change, because both are baked into the files:

- **No transcript.** The INFO is `ALLELE|SYMBOL|DS_AG|...` — one score per gene, computed on
  GENCODE v24 canonical (December 2015). The tool's own help confirms it: `-A grch38` is
  *"GENCODE V24 canonical annotation file in package"*. For 16 of our 519 panel transcripts that
  is a different transcript from the one we report, MET, HRAS, AKT1, FGFR1 and CDKN2A among them.
- **A ±50 bp window**, measured on the files themselves: the largest absolute delta position over
  199,971 records is exactly 50. The current default elsewhere is 500, so any splice effect
  further out is reported as `0.00` — indistinguishable from "no predicted effect".

Running the model takes an annotation as input, so both limits go away together: the scores
describe the transcripts we report, and the distance becomes a parameter.

## How

Two scripts, and no model code of our own.

**`utils/build_spliceai_annotation.py`** builds the gene annotation `spliceai -A` reads. It does
not parse a GTF: a SpliceAI-lookup deployment already ships `gencode.v43.annotation.txt.gz`, which
is GENCODE v43 in SpliceAI's own format with the RefSeq accessions packed into its NAME column, so
the build is a filter. 516 of our 519 panel transcripts are present; the three that are not are
named in the output rather than dropped silently.

```bash
python3 utils/build_spliceai_annotation.py \
    --annotation .../annotations/gencode.v43.annotation.txt.gz \
    --transcripts reference_lists/TSO500_transcripts_list.txt \
    --out .../tso500_grch38.spliceai_annotation.txt
```

Rows are named `GENE--REFSEQ`, following the convention of the laboratory's existing GRCh37
preferred-transcript file. That matters more than it looks: **SpliceAI reports the annotation's
NAME as its SYMBOL field**, so the output states which transcript each score was computed on —
the thing the precomputed tables never recorded.

**`utils/spliceai_compare.py`** writes a VCF from a MAF, runs the model, and compares the result
with the precomputed table variant by variant, reporting how many scores change, how many cross
0.2 (the delta from which the classifier awards O6) or 0.5, and how many have their strongest
signal beyond ±50 bp.

## Which SpliceAI, and a caveat on the validation below

The hosted service at spliceailookup.broadinstitute.org does not run the released package. Its
own page says so:

> Under the hood, this service runs **modified versions** of SpliceAI and Pangolin which have
> added support for multi-nucleotide variants (MNVs), provide access to raw scores in addition to
> delta scores, and have other customizations: https://github.com/bw2/SpliceAI

The work here uses the released `spliceai 1.3.1`, so the agreement demonstrated below holds for a
single-nucleotide substitution, where the fork and the release behave identically. It does **not**
transfer to multi-nucleotide variants, and that is not hypothetical: of the 1,054 variants in the
development cohort, **7 are MNVs and the released package scored none of them**, silently.

```
MPL    chr1:43346917    GG>AA      no score
MPL    chr1:43349337    TG>AA      no score
RAD21  chr8:116851977   CA>TG      no score
RUNX1  chr21:34792613   GC>CT      no score
BCOR   chrX:40064601    TG>CT      no score
BTK    chrX:101353928   GGA>TCC    no score
BTK    chrX:101353930   GGA>TCC    no score
```

Indels are unaffected — 233 of 237 scored — so this is specific to MNVs, and it is 0.7% of the
cohort in genes that matter. **Adoption would therefore require the fork, not the released
package**, or MNVs would go unannotated without saying so, which is the failure mode this whole
line of work exists to remove. It also keeps the two from being interchangeable for verification:
an MNV scored by the fork cannot be checked against the release.

The same page states, in their words, why the precomputed tables are not the thing to build on:

> We no longer use the Illumina precomputed SpliceAI score tables since they are based on an old
> version of Gencode.

## Validation

The decisive test is whether a local run agrees with an independent implementation of the same
model. `16:2317763 T>A` in ABCA3, the Broad's SpliceAI-lookup front-page example:

| | DS_AG | DS_AL | DS_DG | DS_DL |
|---|---|---|---|---|
| Broad SpliceAI-lookup, 500 bp | 0.01 @ 1 | **0.21 @ 106** | 0.04 @ 128 | 0.00 @ 497 |
| this pipeline, 500 bp | 0.01 @ 1 | **0.21 @ 106** | 0.04 @ 128 | 0.00 @ 497 |
| our precomputed table, ±50 bp | 0.01 @ 1 | 0.08 @ −2 | 0.00 @ 2 | 0.00 @ −2 |

Identical to the decimal on all four deltas and all four positions, and the output names the
transcript: `ABCA3--NM_001089.3`. It also reproduces the consequence — the strongest signal sits
at 106 bp, outside the precomputed window, so the table reports 0.08 and the computation 0.21,
crossing the threshold at which we flag a row and award O6.

## The validation is now a test

The comparison above is asserted rather than described. `tests/verification5/test_spliceai_lookup_agreement.py`
(10 tests) checks our annotation against the Broad's service for six public variants whose answers
are frozen in `tests/verification5/data/spliceai_lookup_frozen.json`, with the endpoint and the
date recorded in the fixture.

Frozen rather than queried on each run, for three reasons: the service is rate limited and asks
that batch work be done locally; it returned HTTP 500 on four of 93 rows during this work, so a
live test would fail for reasons that are not ours; and a dated record is what makes a later
difference meaningful — if the values stop matching, either our annotation changed or their
Gencode release did. `utils/freeze_spliceai_lookup.py` refreshes it deliberately.

The six variants each exercise one property: ABCA3 the window, MET and CDKN2A the transcript,
MSH5 a position with three genes and nine transcripts, BRAF a coding hotspot with no splice
signal, IDH1 a control that agrees with the precomputed table. Only public variants are listed,
so nothing patient-derived reaches a third party.

Three groups of test. The first asserts what the frozen answers say about the published tool —
that ABCA3's strongest signal is beyond 50 bp and reaches 0.2, that MET and CDKN2A come back with
both transcripts labelled, and that the raw REF/ALT site scores are present, since those are the
evidence issue #21 needs. The second locks in how the precomputed tables diverge, so replacing
them would be noticed. The third runs the model locally and requires agreement to two decimals,
and it is the one the branch rests on; it skips where the model or the reference FASTA is absent.

## Cost

On an M-series Mac, CPU only: **0.71 s per variant** for 18 variants including model load, with
inference itself about 44 ms per model and five models per variant. A 300-variant panel sample is
therefore a few minutes rather than the microseconds a `tabix` lookup costs.

Against that, the reference data shrinks from **91 GB to about 8 GB**. SpliceAI's precomputed
tables are the single largest component of our 316 GB reference set, ahead of CADD; the model plus
a GENCODE annotation is 8 GB and scores any variant, on any transcript in the annotation, at any
distance.

## What this does not yet show

On the 18 variants of `verification1`, nothing material changed: 4 scores differed, **0 crossed
0.2 or 0.5**, and although 8 had their strongest signal beyond ±50 bp, those scores were 0.01–0.04
— noise at distance rather than signal. That cohort is coding hotspots, which is not where this
matters. Whether it matters for us depends on a cohort with splice-region and intronic variants,
and that measurement has not been made.

Three variants received no local score: the EGFR in-frame deletion, and the GNAQ and DICER1
records whose REF allele does not match GRCh38 ([#19](https://github.com/WeTGI-colab/SMART/issues/19)).
The second of those is the fourth independent detector of that fixture defect, after the SpliceAI
lookup's omission, `verification1/verify.py`'s SKIP and the VEP API's 400.

## Not decided

Nothing here changes the pipeline. Before it could:

- the per-variant latency needs measuring on the hardware that would run it, not on a laptop;
- a cohort with splice-region variants needs scoring both ways, to turn "more correct" into a number;
- the model and annotation versions would have to be pinned and recorded per run, replacing the
  reproducibility a static file gave for free;
- the benchmark figures in the manuscript were produced with the precomputed scores, so changing
  the source means re-running them and saying so;
- the three unaddressable transcripts need a stated fallback rather than a silent gap;
- and the fork would have to replace the released package, for the MNV reason above, which
  means depending on one person's repository rather than Illumina's release.
