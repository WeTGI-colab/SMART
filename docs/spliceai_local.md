# Computing SpliceAI locally

**This is what SMART does by default.** `entrypoint.sh` sets `SPLICEAI_MODE="compute"`, runs the
model per sample against an annotation built from your `--transcripts-file`, and the image ships
`tensorflow-cpu==2.15.*` and `spliceai==1.3.1` to do it. The precomputed tables are still reachable
with `--spliceai lookup`.

This page is the detail behind that default, for [issue
#2](https://github.com/WeTGI-colab/SMART/issues/2): why computing beats looking up, what was
measured, what the measurements do **not** cover, and the one known regression that ships with it.
It began as exploratory work on a branch called `spliceai-local-compute`; that branch has been
merged and deleted, and the findings below are kept because they are the evidence for the default,
not a proposal for one. The parts that are still open are collected under
[Open questions](#open-questions) at the end.

## Why compute rather than look up

The alternative — `--spliceai lookup`, and SMART's behaviour before this work — reads Illumina's
precomputed score files. Those have two limits that no configuration can change, because both are
baked into the files:

- **No transcript.** The INFO is `ALLELE|SYMBOL|DS_AG|...` — one score per gene, computed on
  GENCODE v24 canonical (December 2015). The tool's own help confirms it: `-A grch38` is
  *"GENCODE V24 canonical annotation file in package"*. For 16 of our 519 panel transcripts that
  is a different transcript from the one we report, MET, HRAS, AKT1, FGFR1 and CDKN2A among them.
- **A ±50 bp window**, measured on the files themselves: the largest absolute delta position over
  199,971 records is exactly 50. The current default elsewhere is 500, so any splice effect
  further out is reported as `0.00` — indistinguishable from "no predicted effect".

Running the model takes an annotation as input, so both limits go away together: the scores
describe the transcripts we report, and the distance becomes a parameter. That is the whole
argument for the default.

## How

Two scripts, and no model code of our own.

**`scripts/build_spliceai_annotation.py`** builds the gene annotation `spliceai -A` reads. It does
not parse a GTF: a SpliceAI-lookup deployment already ships `gencode.v43.annotation.txt.gz`, which
is GENCODE v43 in SpliceAI's own format with the RefSeq accessions packed into its NAME column, so
the build is a filter. 516 of our 519 panel transcripts are present; the three that are not are
named in the output rather than dropped silently.

```bash
python3 scripts/build_spliceai_annotation.py \
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

## How it runs in the pipeline

`entrypoint.sh` takes two flags:

```
--spliceai compute|lookup    compute against the reported transcripts (default) or look up
                             the precomputed tables
--spliceai-mask 0|1          raw (0) or masked (1). Default 1
```

**The defaults are `compute` and `mask=1`.** Masked because that is what the model's developers
recommend for variant interpretation and what this laboratory uses; computed because the tables
carry no transcript and stop at ±50 bp.

How it fits together, with no change to any column name:

1. Once per run, the annotation is built from the same `--transcripts-file` already used for
   prioritisation, so a score describes the transcript the report names.
2. Per sample, the model runs on the file VEP is about to annotate — after liftover, so in the
   assembly the reference describes.
3. VEP runs **without** the SpliceAI plugin in compute mode.
4. After MafAnnotator, the scores are merged into the `SpliceAI_pred_*` columns. The names are the
   plugin's, so the tier split, the `ANNOTATION_FLAGS` threshold, O6 and the cryptic-site check
   behave the same whichever mode produced them.

Rows the model did not score are left **empty** rather than carrying a precomputed value, so the
table never mixes the two sources silently in one column.

`compute` requires `gencode.v43.annotation.txt.gz` under
`$REF_DIR/SpliceAI-lookup/references_files/annotations/`, and fails with a message naming
`--spliceai lookup` as the alternative rather than falling back silently. The image ships
`tensorflow-cpu==2.15.*` and `spliceai==1.3.1`; the TensorFlow version is pinned because the
released package uses `tf.keras`, which Keras 3 removes.

**MNVs are reported, not hidden.** The released model silently returns nothing for a
multi-nucleotide variant, so the run counts them and warns per sample. In the development cohort
that is 7 of 1,054. Using the Broad's fork would fix it; see [Open questions](#open-questions).

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

**Correction to a figure quoted earlier in this work.** I repeatedly gave the comparison as
91 GB against about 8 GB. That was wrong: 8 GB is the whole SpliceAI-lookup reference bundle,
which includes its own copies of hg19 and hg38. Compute mode uses the pipeline's existing
`$REF_DIR/liftover/hg38.fa`, so the only additional reference file is
`gencode.v43.annotation.txt.gz` at **10 MB**, with the model riding in the Docker image. The real
comparison is **91 GB against 10 MB**.


On an M-series Mac, CPU only: **0.71 s per variant** for 18 variants including model load, with
inference itself about 44 ms per model and five models per variant. A 300-variant panel sample is
therefore a few minutes rather than the microseconds a `tabix` lookup costs.

Against that, the reference data compute mode *reads* is **10 MB rather than 91 GB**. SpliceAI's
precomputed tables are the single largest component of our 316 GB reference set, ahead of CADD; the
model ships in the image and the GENCODE annotation is 10 MB, and together they score any variant,
on any transcript in the annotation, at any distance.

**That saving is not yet realised on disk.** `$SPLICEAI_SNV` and `$SPLICEAI_INDEL` are in
`entrypoint.sh`'s `REQUIRED_PATHS`, which is checked before anything runs and with no reference to
`SPLICEAI_MODE`, so a default run still refuses to start without the 91 GB present even though it
will not read it. Making the requirement mode-dependent is the obvious follow-up and is listed
under [Open questions](#open-questions). Until then, treat the figure as "what the annotation
costs to compute", not "what you can delete".

## What this does not yet show

On `verification1`, nothing material changed: 4 scores differed, **0 crossed
0.2 or 0.5**, and although 8 had their strongest signal beyond ±50 bp, those scores were 0.01–0.04
— noise at distance rather than signal. That cohort is coding hotspots, which is not where this
matters. Whether it matters for us depends on a cohort with splice-region and intronic variants,
and that measurement has not been made.

Three variants received no local score when this was measured: the EGFR in-frame deletion, and
the GNAQ and DICER1 records whose REF allele did not match GRCh38
([#19](https://github.com/WeTGI-colab/SMART/issues/19)) — the fourth independent detector of that
fixture defect, after the SpliceAI lookup's omission, `verification1/verify.py`'s SKIP and the VEP
API's 400. **The fixture has since been corrected**, so the two substitutions now score; an indel
still does not, which is the model's own limit and not a defect. The measurement above was taken
before the correction and has not been repeated.

## The regression this default ships with

**This ships a measured regression until [#21](https://github.com/WeTGI-colab/SMART/issues/21) is
fixed.** On the development cohort, computing rather than looking up moved three TET2 variants from
Likely Oncogenic to VUS — 391 true positives to 388 — not because the scores are wrong but because
supplying scores that were previously absent trips the O2 cryptic-site check on an acceptor gain at
+1 bp. Masking does not prevent it. So the gain from computing is real and currently masked by a
rule that penalises it, and the two should be settled together.

## Open questions

Compute is the shipped default; these are the things that were *not* settled by shipping it, kept
as a list because each one is a reason to read a computed score with a particular caution.

- **The MNV gap needs the fork.** The released package silently scores no multi-nucleotide variant.
  The run warns per sample, so the gap is visible rather than silent, but it is still a gap.
  Replacing the release with [the Broad's fork](https://github.com/bw2/SpliceAI) would close it and
  means depending on one person's repository rather than Illumina's release — which is why it has
  not been done yet.
- **[#21](https://github.com/WeTGI-colab/SMART/issues/21) is open**, so the measured regression in
  the section above is live.
- **The 91 GB is still required on disk** although compute mode does not read it; see the note
  under [Cost](#cost).
- **Per-variant latency has only been measured on a laptop.** 0.71 s per variant on an M-series Mac,
  CPU only. The hardware that runs this in service has not been measured.
- **No cohort with splice-region and intronic variants has been scored both ways.** `verification1`
  is coding hotspots and showed nothing material (see above), which is the wrong cohort for the
  question. Until that measurement exists, "more correct" is an argument from the model's inputs,
  not a number.
- **The model and annotation versions are pinned in the image but not recorded per run.**
  `tensorflow-cpu==2.15.*` and `spliceai==1.3.1` are in the Dockerfile and `gencode.v43` is in the
  annotation's filename, but `logs/run_provenance.txt` records the mode and the mask, not those
  versions — so a computed score is reproducible from the image tag, which is weaker than the
  reproducibility a static table gave for free.
- **The manuscript's benchmark figures were produced with the precomputed scores.** Re-running them
  against computed scores, and saying which was used, is outstanding.
- **Three of the 519 panel transcripts are absent from the GENCODE v43 annotation.**
  `build_spliceai_annotation.py` names them rather than dropping them silently, but there is no
  stated fallback for a variant on one of them.
