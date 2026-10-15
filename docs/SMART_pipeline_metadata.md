# SMART Pipeline Metadata
## For Data Submission — University Hospital Southampton

---

## 1. Pipeline Overview

**Pipeline name:** SMART — Somatic Mutation Annotation and Reporting Tool  
**Version:** 2.1.0 (the single source is the repository's `VERSION` file; the pipeline stamps it
into the first line of every Tier 1 MAF)  
**Repository:** https://github.com/WeTGI-colab/SMART  
**DOI:** https://doi.org/10.5281/zenodo.20206503  
**Containerisation:** Docker (`monkiky/smart:2.1.0`; images are published as
`monkiky/smart:<version>`)  
**Contact:** open an issue at https://github.com/WeTGI-colab/SMART/issues

SMART is a Dockerised pipeline for somatic variant annotation and clinical reporting. Starting from
somatic VCF files, it performs PASS filtering, optional coordinate liftover (hg19 → GRCh38), a REF
check against the reference genome, SpliceAI computation against the reported transcripts,
comprehensive functional annotation via Ensembl VEP, clinical annotation via OncoKB, and produces
analysis-ready output tables in four audience-targeted tiers.

**Every version and file name in this document is recorded per run**, in
`output/logs/run_provenance.txt` and `output/logs/<sample>.vcf_metadata.json`. Where a figure below
is stated as a fact about a particular run, it is taken from the committed run output under
`tests/verification1/output/` and is asserted by
`tests/verification5/test_docs_claims.py`. **Read the run's own provenance file in preference to
this document**: it cannot go stale.

---

## 2. Reference Genome

| Parameter | Value |
|-----------|-------|
| Assembly | GRCh38 / hg38 |
| Source | UCSC Genome Browser |
| File | hg38.fa (`$REF_DIR/liftover/hg38.fa`) |
| Version identifier | UCSC_hg38 |
| Coordinates | All variants annotated and reported in GRCh38 |
| Input coordinates | GRCh38 (liftover from hg19 available if required) |

---

## 3. Pipeline Steps and Parameters

The steps below are in execution order. Steps 1, 3–9 and 11–12 are on by default; steps 2 and 10
are conditional as stated.

### Step 1 — PASS Filter

| Parameter | Value |
|-----------|-------|
| Tool | `zcat` \| `awk` \| `bgzip`, then `tabix` to index (no bcftools — the filter is one `awk` expression on column 7) |
| Action | Retain only variants with FILTER = PASS |
| Flag | `--pass` (default ON) / `--no-pass` to disable |

### Step 2 — Coordinate Liftover (optional)

| Parameter | Value |
|-----------|-------|
| Tool | GATK LiftoverVcf v4.6.0.0 |
| Chain file | hg19ToHg38.over.chain (UCSC) |
| Action | Convert hg19/GRCh37 → GRCh38 coordinates |
| Flag | `--liftover` (default ON) / `--no-liftover` to skip |

### Step 3 — SpliceAI (computed)

| Parameter | Value |
|-----------|-------|
| Tool | `spliceai` 1.3.1 on `tensorflow-cpu` 2.15, both shipped in the image |
| Action | Run the model on the file VEP is about to annotate, against a gene annotation built once per run from `--transcripts-file`, so each score names the transcript it was computed on |
| Distance / masking | `-D 500`, `-M 1` (masked) |
| Flags | `--spliceai compute` (default) / `--spliceai lookup` for Illumina's precomputed tables; `--spliceai-mask 0\|1` (default 1) |
| Caveat | The released model returns nothing for a multi-nucleotide variant; the run counts and warns about them per sample rather than leaving the gap silent. See [`spliceai_local.md`](spliceai_local.md) |

### Step 4 — Reference (REF) Check

| Parameter | Value |
|-----------|-------|
| Script | `scripts/check_reference.py` |
| Action | Compare every checkable REF allele of the VCF VEP is about to read against hg38, through its `.fai` |
| Flag | `--ref-check fail` (default, stops the run) / `--ref-check warn` |
| Why | VEP does not reject a REF that disagrees with the FASTA; it silently annotates a variant nobody submitted. See [`reference_check.md`](reference_check.md) |

### Step 5 — Functional Annotation

| Parameter | Value |
|-----------|-------|
| Tool | Ensembl VEP v114.0 (the Dockerfile pins `VEP_VERSION=114.0`; the annotated VCF's `##VEP=` header records what actually ran) |
| Cache | `homo_sapiens_merged/114_GRCh38`, offline. **Merged, not `homo_sapiens/`** — the pipeline runs `--merged`, which reads the merged Ensembl+RefSeq cache and will not find a plain one |
| Transcript selection | MANE Select / MANE Plus Clinical (preferred); fallback to first transcript |
| Custom transcript whitelist | User-supplied NM_ accession list |
| Output format | VCF (uncompressed) |
| Key flags | `--vcf`, `--everything`, `--offline`, `--merged`, `--assembly GRCh38`, `--no_escape`, `--fasta`, `--dir`, `--dir_plugins` |

Three flags an earlier version of this document listed are **not** passed, and listing them was
misleading about which cache a reader needs: `--everything` already implies canonical and MANE
selection, and `--offline` already implies cache use. What matters and is easy to miss is
`--merged`, which is passed and decides that the cache directory must be `homo_sapiens_merged/`.
The `##VEP=` header of each annotated VCF records the version, cache path and assembly that
actually ran.

**VEP Plugins and custom annotations:**

| Plugin / Annotation | Version | Description |
|---------------------|---------|-------------|
| SpliceAI | plugin not used in the default `compute` mode (see step 3); v1.3 tables in `lookup` mode | Splice-site impact prediction (delta scores for AG/AL/DG/DL) |
| REVEL | v1.3 | Missense pathogenicity meta-score (0–1); run with `no_match=1` |
| LOEUF | gnomAD v4.1 | Loss-of-function constraint, **per transcript, not per gene** (93% of multi-transcript genes differ; looked up with `match_by=transcript`). Autosomes only -- v4.1 publishes no chrX/chrY constraint |
| dbNSFP | 5.0a (subset built by `utils/build_dbnsfp.py`) | Missense and LoF predictor ensemble, run with `transcript_match=1` so the values are the reported transcript's. Optional (academic licence) |
| CADD | whole-genome SNVs + gnomAD r4.0 indels | Genome-wide deleteriousness score. Optional (~82 GB) |
| ClinVar | `clinvar.vcf.gz`, VEP `--custom` | Clinical significance + somatic oncogenicity (ClinVar_ONC, ClinVar_SCI) |
| CIViC | `civic_grch38.vcf.gz`, VEP `--custom` (built by `utils/civic_formating.py`) | Clinical variant evidence |
| COSMIC CMC | `cosmic_cmc_grch38.vcf.gz`, VEP `--custom` | Cancer Mutation Census counts and tiers. Optional (COSMIC licence) |
| CancerHotspots | changv2_gao_nc (GRCh38) | Statistically significant somatic mutational hotspots |
| gnomAD | v4.1 (via VEP cache) | Population allele frequencies (exomes + genomes) |

### Step 6 — Clinical Annotation

| Parameter | Value |
|-----------|-------|
| Tool | `scripts/oncokb2.0.py` against the OncoKB public API |
| OncoKB data version | Queried from `/api/v1/info` at the start of every run and recorded in `logs/run_provenance.txt`; also written per record as `ONCOKB_DATA_VERSION`. The committed run records **v7.6 (09/30/2026)** |
| API version | Queried the same way; the committed run records **v1.7.0** |
| Tumour type mode | `generic` (pan-cancer; omits tumour type from API query) |
| Endpoints used | `/annotate/mutations/byProteinChange` (SNV/indel), `/annotate/mutations/byGenomicChange` (fallback when the reported transcript is not OncoKB's), `/annotate/copyNumberAlterations` (CNA) |

### Step 7 — VCF to Table

| Parameter | Value |
|-----------|-------|
| Script | `scripts/vcf2table.py` |
| Action | Parse VEP-annotated VCF to structured CSV |
| Transcript selection | 3-tier logic: (1) preferred NM whitelist, (2) MANE Select/Plus Clinical, (3) first VEP transcript |
| Multi-transcript output | One row per preferred transcript when variant overlaps multiple whitelist isoforms |
| Provenance columns | `TRANSCRIPT_REQUESTED` and `TRANSCRIPT_MATCH` state per row whether the annotated transcript is one you supplied. See [`transcript_specificity.md`](transcript_specificity.md) |

### Step 8 — MAF Standardisation

| Parameter | Value |
|-----------|-------|
| Tool | OncoKB MafAnnotator |
| Action | Convert CSV to MAF format; standardise column names |

### Step 9 — VCF Metadata Harvest

| Parameter | Value |
|-----------|-------|
| Script | `scripts/vcf_metadata.py` |
| Action | Read the `##INFO`/`##FORMAT`/`##VEP`/plugin headers of the input, annotated and OncoKB VCFs into `logs/<sample>.vcf_metadata.json`, before the VCFs are cleaned up |
| Why | The description, source and version in the Tier 2/3 header rows then come from the run that produced the data rather than from a hand-maintained list. See [`field_provenance.md`](field_provenance.md) |

### Step 10 — GENIE Case Counts (conditional)

| Parameter | Value |
|-----------|-------|
| Script | `scripts/genie_annotation.py` |
| Action | Annotate each variant with AACR Project GENIE case counts, used by the classifier for O4 |
| Condition | Runs when `$REF_DIR/GENIE/genie_o4_counts.tsv.gz` exists. The table is user-built (`utils/build_genie_o4.py`) because GENIE is under a data-use agreement |
| Flag | `--o4-on-target code\|tissue\|lineage\|all` (default `lineage`) |

### Step 11 — Post-analysis and Tiering

| Parameter | Value |
|-----------|-------|
| Script | `scripts/post_analysis.py` |
| Action | Merge per-sample MAFs; expand JSON fields; add `SPLICEAI_MAX_DS` and `ANNOTATION_FLAGS`; apply tier filtering; build the two-row header from the harvested VCF metadata |
| Output files | `Final_result_tier1.maf` (all fields), `Final_result_tier2.tsv` (bioinformatics), `Final_result_tier3.tsv` (clinical) |

### Step 12 — Tier 4: SVIG-UK Oncogenicity Classification

| Parameter | Value |
|-----------|-------|
| Script | `scripts/svig_uk_classifier.py` |
| Action | Read the Tier 1 MAF and score each variant against the ACGS/SVIG-UK somatic oncogenicity codes, writing `Final_result_tier4.tsv` and `Final_result_tier4.json` |
| Tumour type | `--tumour-type` takes an OncoTree code (AML, MDS, LUAD…) or a broad group (`myeloid`, `lymphoid`, `haem`, `solid`, `unknown` — default `unknown`). A per-variant `TUMOUR_TYPE` column in the input takes precedence |
| Missing resources | Classifier resources under the reference directory are optional; a code whose resource is absent is reported as not evaluated in `SVIG_UK_NOTES` rather than scored as absent evidence |
| Flag | On by default; `--no-classification` to skip |

---

## 4. Software Versions

Recorded per run in `output/logs/run_provenance.txt`. The values below are those of the committed
run under `tests/verification1/output/`.

| Software | Version | Purpose |
|----------|---------|---------|
| SMART | 2.1.0 | Somatic variant annotation pipeline |
| Docker | 29.3.1 | Container runtime |
| Ensembl VEP | 114.0 | Functional annotation |
| GATK | 4.6.0.0 | Coordinate liftover |
| bcftools | 1.13 | VCF manipulation |
| SpliceAI | 1.3.1 (on tensorflow-cpu 2.15) | Splice prediction, computed per transcript |
| OncoKB API | v1.7.0 (queried per run) | Clinical annotation |
| OncoKB Annotator | pinned at image build | MAF standardisation |
| samtools | system (Docker image) | FASTA/BAM utilities |
| tabix / bgzip | system (Docker image) | VCF indexing |
| Python | 3.10.12 | Pipeline scripts |
| pandas | ≥2.0 | Data processing |
| cyvcf2 | ≥0.31 | VCF parsing |

---

## 5. Reference Databases

The file name of each is recorded per run in `output/logs/run_provenance.txt`, and the resource
version that annotated each column is recorded per run in
`output/logs/<sample>.vcf_metadata.json` and reproduced in the Tier 2/3 header rows.

| Database | Version / File | Source | Used for | Required? |
|----------|---------------|--------|----------|-----------|
| GRCh38 reference genome | UCSC hg38 (`liftover/hg38.fa`) | UCSC Genome Browser | Annotation reference, REF check, SpliceAI, liftover target | yes |
| hg19→hg38 liftover chain | `liftover/hg19ToHg38.over.chain` | UCSC Genome Browser | Coordinate conversion | yes |
| Ensembl VEP cache | `homo_sapiens_merged/114_GRCh38` | Ensembl FTP | Functional annotation | yes |
| SpliceAI GENCODE annotation | `gencode.v43.annotation.txt.gz` (10 MB) | SpliceAI-lookup reference bundle | `--spliceai compute`: names the transcript each score describes | yes, in compute mode |
| SpliceAI precomputed scores | `spliceai_scores.raw.{snv,indel}.hg38.vcf.gz` (91 GB) | Illumina BaseSpace | `--spliceai lookup` | yes (checked in both modes) |
| REVEL scores | `new_tabbed_revel_grch38.tsv.gz` | REVEL website | Missense pathogenicity | yes |
| ClinVar | `clinvar.vcf.gz` | NCBI FTP | Germline + somatic clinical significance | yes |
| CIViC | `civic_grch38.vcf.gz` | civicdb.org | Clinical variant evidence | yes |
| gnomAD (via VEP cache) | v4.1 | Broad Institute | Population allele frequencies | yes |
| gnomAD LOEUF constraints | v4.1 (`loeuf_dataset_grch38.tsv.gz`), **per transcript** | Broad Institute | Transcript constraint scores | yes |
| CancerHotspots | changv2_gao_nc, GRCh38 (Chang et al. 2018) | cancerhotspots.org | Somatic mutational hotspots | yes |
| COSMIC CMC | `cosmic_cmc_grch38.vcf.gz` | COSMIC (licence) | Cancer Mutation Census counts and tiers | optional |
| CADD | `whole_genome_SNVs.tsv.gz`, `gnomad.genomes.r4.0.indel.tsv.gz` | CADD website | Genome-wide deleteriousness | optional |
| dbNSFP | `dbNSFP5.0a_grch38.smart.gz` | dbNSFP (academic licence) | Predictor ensemble, per transcript | optional |
| AACR Project GENIE | `genie_o4_counts.tsv.gz` (user-built) | AACR (data-use agreement) | Case counts for SVIG-UK code O4 | optional |
| OncoTree | OncoTree release JSON | oncotree.mskcc.org | Resolving `--tumour-type` to a code, tissue and lineage for Tier 4 | optional |
| SVIG-UK canonical variants | `svig_uk_canonical_variants.tsv` | ACGS/SVIG-UK | Code O1, and the `svig_canonical` annotation flag | optional |
| OncoKB | queried live; the committed run records v7.6 (09/30/2026) | oncokb.org | Therapeutic/diagnostic/prognostic levels | yes (token) |

---

## 6. Output Files

Each pipeline run produces the following per cohort, under the directory passed to
`--output-dir`. The column counts are those of the committed run under
`tests/verification1/output/` and vary with the cohort, because the OncoKB evidence columns are
padded to the widest row in the run.

| File | Format | Description |
|------|--------|-------------|
| `output/Final_result_tier1.maf` | MAF | All non-redundant annotation fields (1,165 columns in the committed run). First line: `#SMART_VERSION 2.1.0` |
| `output/Final_result_tier2.tsv` | TSV | Selected fields for bioinformaticians (700 columns). Two-row header (field names + source/version metadata) |
| `output/Final_result_tier3.tsv` | TSV | Clinically relevant fields (78 columns). Two-row header |
| `output/Final_result_tier4.tsv` | TSV | SVIG-UK oncogenicity classification: class, score, and the codes awarded, per variant |
| `output/Final_result_tier4.json` | JSON | The same classification with the per-code evidence that produced it |
| `variant_counts.txt` | TXT | Per-sample variant counts at each pipeline stage |
| `logs/run_provenance.txt` | TXT | Every tool version, reference file and configuration flag for this run |
| `logs/<sample>.log` | TXT | Full per-sample log: the VEP command, OncoKB call details, counts, warnings |
| `logs/<sample>.vcf_metadata.json` | JSON | The field descriptions, sources and versions harvested from this sample's VCF headers |
| `AnnotatedVcf/` | VCF (uncompressed) and `_summary.html` | VEP-annotated VCFs, kept with `--keep-tmp` |
| `SpliceAI/`, `FilteredVcf/`, `OncoKB_VCF/`, `Table/`, `FINAL_Table/` | — | Intermediates, kept with `--keep-tmp` / `--keep-tables` |

---

## 7. Exact Command Lines Used

`--output-dir` is **required**: the run aborts with `ERROR: --output-dir is required.` without it,
because without a mounted volume the results are lost when the container exits. The command below
is the one `tests/run_all_verifications.sh` issues, which produced the committed output this
document cites.

```bash
docker run --rm \
  -v /path/to/SMART/tests/verification1:/data \
  -v /path/to/SMART/tests/verification1/output:/results \
  -v /path/to/refs:/refs:ro \
  monkiky/smart:2.1.0 \
  <ONCOKB_TOKEN> \
  --transcripts-file /data/verification1_transcripts.txt \
  --config /data/Config.yaml \
  --ref-dir /refs \
  --output-dir /results \
  --input-dir /data \
  --no-liftover \
  --keep-tmp \
  --keep-tables
```

**Parameter explanation:**

| Flag | Value | Meaning |
|------|-------|---------|
| `--output-dir` | /results | **Required.** Where every output directory is created. Mount a host volume here |
| `--ref-dir` | /refs | **Required.** The reference set |
| `--transcripts-file` | verification1_transcripts.txt | NM_ accession whitelist for transcript prioritisation. Required unless `--vep-only` |
| `--config` | Config.yaml | Field tier definitions (controls which columns appear in each output tier). Required unless `--vep-only` |
| `--no-liftover` | ON | Input VCFs are already in GRCh38 — liftover step skipped |
| `--pass` | ON (default) | Only PASS-filtered variants retained |
| `--ref-check` | fail (default) | A REF allele that disagrees with hg38 stops the run |
| `--spliceai` | compute (default) | SpliceAI computed against the reported transcripts |
| `--spliceai-mask` | 1 (default) | Masked scores |
| `--dbnsfp-transcript-match` | 1 (default) | dbNSFP predictors for the reported transcript only |
| `--tumour-type` | unknown (default) | Tumour type for the Tier 4 classification |
| `--o4-on-target` | lineage (default) | Which GENIE entries count towards O4 |
| `--keep-tmp` | ON | Intermediate VCFs kept after the run |
| `--keep-tables` | ON | Per-sample intermediate MAF tables kept after merging |
| `--jobs` | 1 (default) | Samples processed sequentially |

**Pipeline configuration recorded at runtime:**

Every SMART run writes this to `output/logs/run_provenance.txt` as well as to the log, so a
run's versions and settings are recoverable from its own output afterwards. This is the committed
run's file, verbatim:

```
============================================================
SMART — Somatic Mutation Annotation and Reporting Tool  v2.1.0
============================================================
Tool versions:
  VEP:               114.0
  GATK:              4.6.0.0
  bcftools:          bcftools 1.13
  Python:            Python 3.10.12
  OncoKB API:        v1.7.0

Reference versions:
  ClinVar:           clinvar.vcf.gz
  CIViC:             civic_grch38.vcf.gz
  CancerHotSpots:    hg38.hotspots_changv2_gao_nc.vcf.gz
  COSMIC CMC:        cosmic_cmc_grch38.vcf.gz
  SpliceAI SNV:      spliceai_scores.raw.snv.hg38.vcf.gz
  SpliceAI INDEL:    spliceai_scores.raw.indel.hg38.vcf.gz
  REVEL:             new_tabbed_revel_grch38.tsv.gz
  CADD:              whole_genome_SNVs.tsv.gz
  dbNSFP:            dbNSFP5.0a_grch38.smart.gz
  GENIE counts:      genie_o4_counts.tsv.gz
  LOEUF:             loeuf_dataset_grch38.tsv.gz
  OncoKB data:       v7.6 (09/30/2026)


Config:
  PASS filter:      ENABLED
  On REF mismatch:  stop the run (--ref-check fail)
  dbNSFP values:    the reported transcript's only (transcript_match=1)
  SpliceAI:         computed against the reported transcripts, masked scores
  Liftover:         DISABLED
  VEP only:         NO
  Clean tmp:        DISABLED
  Clean tables:     DISABLED
  Parallel jobs:    1
  Transcript file:  /data/verification1_transcripts.txt
  Config file:      /data/Config.yaml
  Reference dir:    /refs
  Input dir:        /data
============================================================
```

---

## 8. Reproducibility

The entire pipeline runs inside a single Docker container (`monkiky/smart:2.1.0`):

```bash
docker pull monkiky/smart:2.1.0
```

To read the version out of the image, and out of any SMART output file:

```bash
docker run --rm --entrypoint cat monkiky/smart:2.1.0 /opt/smart/VERSION
# → 2.1.0

head -1 Final_result_tier1.maf
# → #SMART_VERSION 2.1.0
```

Three further records make a run reproducible without reference to this document:

- `output/logs/run_provenance.txt` — every tool version, reference file and flag, as shown above;
- the `##VEP=` header of the annotated VCF — VEP's own version, API, cache path, assembly and the
  release of every source VEP itself supplies (gnomAD, dbSNP, GENCODE, SIFT, PolyPhen);
- `output/logs/<sample>.vcf_metadata.json` — the description, source and version of every field,
  read from the headers of the VCFs that carried it.

All pipeline parameters are controlled via command-line flags. A minimal run command:

```bash
docker run --rm \
  -v /path/to/data:/data \
  -v /path/to/results:/results \
  -v /path/to/refs:/refs:ro \
  monkiky/smart:2.1.0 \
  <ONCOKB_TOKEN> \
  --transcripts-file /data/transcripts.txt \
  --config /data/Config.yaml \
  --ref-dir /refs \
  --output-dir /results \
  --input-dir /data \
  --no-liftover
```

---

*Document generated: May 2026. Last checked against the code and the committed run output:
October 2026.*  
*SMART v2.1.0 · University Hospital Southampton*
