# Where each column's description, source and version come from

Every column of the SMART outputs carries three pieces of metadata in the second header row of
Tiers 2 and 3: a description, the source it came from, and the version of that source. This note
says where each of the three is taken from, and why.

Prompted by [#16](https://github.com/WeTGI-colab/SMART/issues/16), which reported that the variant
caller's `hotspot` flag was attributed to the wrong resource.

## The problem this solves

Before this, all three came from `Config.yaml` alone — a hand-maintained list of field definitions.
A hand-maintained list drifts from what the pipeline actually does. **The four figures below are the
measurements that prompted the change, taken at that time; they are history, not the current
state.** The current state is in the next section but one.

- `hotspot` was described as a Cancer Hotspots flag, version `changv2`. It is written by the
  variant caller into the input VCF. A clinical scientist reading Tier 3 would have attributed
  it to the wrong database.
- `FractionInformativeReads` was attributed to Manta. DRAGEN writes it.
- Six further fields DRAGEN writes had no entry at all.
- Of the 807 columns in the Tier 1 MAF of the run measured then, 475 had no entry.

None of this was anyone's carelessness. It is what happens to a list that has to be updated by a
person every time a resource, a caller or a plugin changes.

## Where the metadata actually lives

A VCF documents its own fields, and does so for the run that produced them:

| Written by | What it declares |
|---|---|
| the variant caller | one `##INFO`/`##FORMAT` header per field, plus `##source` naming itself |
| VEP | `##VEP="v114.0" … cache="…"`, and the CSQ field order |
| a VEP custom annotation | one `##INFO` per field, naming **the resource file used in this run** |
| a VEP plugin | a line of its own, e.g. `##SpliceAI_pred_DS_AG=SpliceAI (v1.3) …` |
| OncoKB | one `##INFO` per field it adds, and `ONCOKB_DATA_VERSION` per record |

The caller's headers survive VEP annotation untouched, so the annotated VCF carries all of it.

SMART used to read exactly one of those headers — CSQ, for the field order — and discard the rest.
`scripts/vcf_metadata.py` now harvests them.

## The rule

| | Taken from |
|---|---|
| **description** | the VCF where it has one; `Config.yaml` fills the gaps |
| **source** | the VCF |
| **version** | the VCF |
| **tier** | `Config.yaml` always — it is our editorial decision, not a property of the data |

Two qualifications, both deliberate:

**A generic "VEP" does not override a specific source.** A core CSQ field such as `gnomADe_AF`
has no header of its own, so all the file says is that VEP delivered it. `Config.yaml` naming
gnomAD is more informative and is kept; the version still comes from the run. Without this rule
every population-frequency and predictor column would have been relabelled "VEP".

**Columns no VCF describes keep `Config.yaml`.** That is the MAF standard columns, everything
SMART computes, and the OncoKB MafAnnotator columns, which come from a table and have no
headers — 848 of the 1,165 Tier 1 columns in the committed `tests/verification1` run resolve to
source `OncoKB`, so this is the majority of the sheet. Their *version*, however, is the OncoKB data
version observed in the run: an analysis cannot span two OncoKB releases, so the release is a
property of the run, not of the tool.

## What it currently resolves

Measured on the committed run under `tests/verification1/output/`, by passing each Tier 1 column
through `vcf_metadata.resolve()` with that run's own `logs/verification1.vcf_metadata.json`:

| | Count |
|---|---|
| Tier 1 columns | 1,165 |
| resolving to source `OncoKB` (the MafAnnotator table) | 848 |
| resolving to a named annotation source (VEP, ClinVar, CIViC, gnomAD, CancerHotspots, SpliceAI, COSMIC, CADD, dbNSFP, GENIE, SMART, the input VCF) | 298 |
| resolving to no source at all | 19 |
| `Config.yaml` entries behind that: explicit field definitions | 427 |
| `Config.yaml` entries behind that: wildcard patterns | 36 |

The last two come straight out of the file:

```bash
python3 -c "import yaml; c=yaml.safe_load(open('Config.yaml')); \
print(len(c['fields']), 'fields,', len(c['field_patterns']), 'patterns')"
```

The comparison with "475 of 807 had no entry" is the point of the whole exercise: 19 of 1,165
columns now resolve to nothing, where 475 of 807 did. The column count itself grows with the
cohort, because the OncoKB evidence columns are padded to the widest row in a run, which is why
these figures are stated against a named run rather than in the abstract.
`tests/verification5/test_docs_claims.py` recomputes them and fails if this table drifts.

## What is checked

| Check | Level | Why |
|---|---|---|
| An input `INFO` name that collides with a VEP or plugin field | **error** | `vcf2table` applies INFO after the CSQ fields, so the caller's value silently replaces the annotation's. A wrong value, not a missing one |
| Fields declared in the input VCF that never reach the table | warning | annotation lost somewhere between the caller and the table |
| More than one `ONCOKB_DATA_VERSION` in one run | **error** | the run straddles two OncoKB releases |
| Samples in one run disagreeing on VEP version or a resource file | **error** | the merged tables would describe some of their own rows wrongly |
| Liftover enabled | warning | GATK rewrites records; the caller's INFO fields may not survive |
| No `##VEP` header | warning | VEP version unknown for this run |
| `Config.yaml` claiming a source the VCF contradicts | warning | the VCF is used; the config entry is stale |

`FORMAT` entries are deliberately excluded from the collision and loss checks: `post_analysis`
expands and renames them (the caller's `AF` becomes `VAF`), so they never claim a column of their
own name. This matters — letting the caller's `AF` claim the `AF` column would label a population
frequency as a tumour allele fraction on a clinical sheet.

## How it runs

The annotated and OncoKB VCFs are deleted per sample once the tables are written, so the harvest
happens inside the sample loop, before the cleanup:

```
vcf_metadata.py --input <sample>.vcf.gz \
                --annotated <sample>_annotated.vcf \
                --oncokb <sample>.oncoKB.vcf \
                --out logs/<sample>.vcf_metadata.json
```

`post_analysis.py --vcf-metadata logs` then merges the per-sample sidecars, checks the samples
agree, and builds the two-row header from the result. Without `--vcf-metadata` the pipeline
behaves exactly as before and describes every field from `Config.yaml`.

## Standalone use

```bash
# what a VCF says about itself
python3 scripts/vcf_metadata.py --input sample.vcf.gz

# the full harvest for one sample
python3 scripts/vcf_metadata.py --input s.vcf.gz --annotated s_annotated.vcf \
    --oncokb s.oncoKB.vcf --oncokb-version v7.6 --out s.vcf_metadata.json

# combine a run and report where its samples disagree
python3 scripts/vcf_metadata.py --merge logs/*.vcf_metadata.json --summary
```

## Test fixtures

`tests/verification5/data/dragen_real_subset.vcf` is a de-identified extract of a real run: five
somatic calls with the headers intact. De-identification is part of the contract and is asserted
by a test — the run timestamp, the inferred sex karyotype and a centre name in a file path were
removed, the contigs reduced to those used, and only variants with an allele fraction below 0.4
kept, because a run of heterozygous germline genotypes is a genetic fingerprint.

`dragen_like.vcf` is a hand-written file kept for the cases the real one cannot exercise. It is
also a standing reminder of why this module exists: when its field descriptions were written from
knowledge rather than read from a file, **all eight were subtly wrong** — `RU` was "Repeat unit
for the short tandem repeat" where DRAGEN says "Tandem repeat unit (bases).", and `SoftClipRatio`
was described as a ratio of reads where DRAGEN says it is a ratio of bases.
