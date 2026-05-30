# SMART Pipeline Metadata
## For Data Submission — University Hospital Southampton

---

## 1. Pipeline Overview

**Pipeline name:** SMART — Somatic Mutation Annotation and Reporting Tool  
**Version:** 1.0.1  
**Repository:** https://github.com/Manuel-DominguezCBG/SMART  
**DOI:** https://doi.org/10.5281/zenodo.20206503  
**Containerisation:** Docker (`monkiky/smart:1.0.1`)  
**Contact:** Manuel Dominguez — manolo.biomero@gmail.com

SMART is a Dockerised pipeline for somatic variant annotation and clinical reporting. Starting from somatic VCF files, it performs PASS filtering, optional coordinate liftover (hg19 → GRCh38), comprehensive functional annotation via Ensembl VEP, clinical annotation via OncoKB, and produces analysis-ready output tables in three audience-targeted tiers.

---

## 2. Reference Genome

| Parameter | Value |
|-----------|-------|
| Assembly | GRCh38 / hg38 |
| Source | UCSC Genome Browser |
| File | hg38.fa |
| Version identifier | UCSC_hg38 |
| Coordinates | All variants annotated and reported in GRCh38 |
| Input coordinates | GRCh38 (liftover from hg19 available if required) |

---

## 3. Pipeline Steps and Parameters

### Step 1 — PASS Filter

| Parameter | Value |
|-----------|-------|
| Tool | bcftools / awk |
| Action | Retain only variants with FILTER = PASS |
| Flag | `--pass` (default ON) / `--no-pass` to disable |

### Step 2 — Coordinate Liftover (optional)

| Parameter | Value |
|-----------|-------|
| Tool | GATK LiftoverVcf v4.6.0.0 |
| Chain file | hg19ToHg38.over.chain (UCSC) |
| Action | Convert hg19/GRCh37 → GRCh38 coordinates |
| Flag | `--liftover` (default ON) / `--no-liftover` to skip |

### Step 3 — Functional Annotation

| Parameter | Value |
|-----------|-------|
| Tool | Ensembl VEP v114.2 |
| Cache | homo_sapiens 114 GRCh38 (offline) |
| Transcript selection | MANE Select / MANE Plus Clinical (preferred); fallback to first transcript |
| Custom transcript whitelist | User-supplied NM_ accession list |
| Output format | VCF |
| Key flags | `--everything`, `--canonical`, `--mane`, `--offline`, `--cache` |

**VEP Plugins and custom annotations:**

| Plugin / Annotation | Version | Description |
|---------------------|---------|-------------|
| SpliceAI | v1.3 | Splice-site impact prediction (delta scores for AG/AL/DG/DL) |
| REVEL | v1.3 | Missense pathogenicity meta-score (0–1) |
| LOEUF | gnomAD v4.0 | Gene-level loss-of-function constraint score |
| ClinVar | 2026-03-21 | Clinical significance + somatic oncogenicity (ClinVar_ONC, ClinVar_SCI) |
| CIViC | nightly 2026-03-29 | Clinical variant evidence (converted to VCF) |
| CancerHotspots | changv2_gao_nc (GRCh38) | Statistically significant somatic mutational hotspots |
| gnomAD | v4.1 (via VEP cache) | Population allele frequencies (exomes + genomes, 807,162 individuals) |

### Step 4 — Clinical Annotation

| Parameter | Value |
|-----------|-------|
| Tool | OncoKB Annotator (MafAnnotator) |
| OncoKB data version | v7.0 (as of run date — version recorded per run in output MAF header) |
| API version | v5.4 |
| Tumour type mode | `generic` (pan-cancer; omits tumour type from API query) |
| Endpoints used | `/annotate/mutations/byProteinChange` (SNV/indel), `/annotate/copyNumberAlterations` (CNA) |

### Step 5 — VCF to Table

| Parameter | Value |
|-----------|-------|
| Script | `scripts/vcf2table.py` (SMART v1.0.1) |
| Action | Parse VEP-annotated VCF to structured CSV |
| Transcript selection | 3-tier logic: (1) preferred NM whitelist, (2) MANE Select/Plus Clinical, (3) first VEP transcript |
| Multi-transcript output | One row per preferred transcript when variant overlaps multiple whitelist isoforms |

### Step 6 — MAF Standardisation

| Parameter | Value |
|-----------|-------|
| Tool | OncoKB MafAnnotator |
| Action | Convert CSV to MAF format; standardise column names |

### Step 7 — Post-analysis and Tiering

| Parameter | Value |
|-----------|-------|
| Script | `scripts/post_analysis.py` (SMART v1.0.1) |
| Action | Merge per-sample MAFs; expand JSON fields; apply tier filtering |
| Output files | `Final_result_tier1.maf` (all fields), `Final_result_tier2.tsv` (bioinformatics), `Final_result_tier3.tsv` (clinical) |

---

## 4. Software Versions

| Software | Version | Purpose |
|----------|---------|---------|
| SMART | 1.0.1 | Somatic variant annotation pipeline |
| Docker | 29.3.1 | Container runtime |
| Ensembl VEP | 114.2 | Functional annotation |
| GATK | 4.6.0.0 | Coordinate liftover |
| OncoKB Annotator | latest (pinned at build) | Clinical annotation |
| bcftools | system (Docker image) | VCF manipulation |
| samtools | system (Docker image) | BAM/VCF utilities |
| tabix / bgzip | system (Docker image) | VCF indexing |
| Python | 3.10 | Pipeline scripts |
| pandas | ≥2.0 | Data processing |
| cyvcf2 | ≥0.31 | VCF parsing |

---

## 5. Reference Databases

| Database | Version / Date | Source | Used for |
|----------|---------------|--------|---------|
| GRCh38 reference genome | UCSC hg38 | UCSC Genome Browser | Alignment reference, liftover target |
| hg19→hg38 liftover chain | UCSC | UCSC Genome Browser | Coordinate conversion |
| Ensembl VEP cache | Release 114 (GRCh38) | Ensembl FTP | Functional annotation |
| SpliceAI scores | v1.3 | Illumina BaseSpace | Splice-site impact (SNV + indel) |
| REVEL scores | v1.3 | REVEL website | Missense pathogenicity |
| ClinVar | 2026-03-21 | NCBI FTP | Germline + somatic clinical significance |
| CIViC | nightly 2026-03-29 | civicdb.org | Clinical variant evidence |
| gnomAD (via VEP) | v4.1 | Broad Institute | Population allele frequencies |
| gnomAD LOEUF constraints | v4.0 | Broad Institute | Gene constraint scores |
| CancerHotspots | v2 (Chang et al. 2018) | cancerhotspots.org | Somatic mutational hotspots |
| OncoKB | v7.0 | oncokb.org | Therapeutic/diagnostic/prognostic levels |

---

## 6. Output Files

Each pipeline run produces the following output files per cohort:

| File | Format | Description |
|------|--------|-------------|
| `Final_result_tier1.maf` | MAF | All non-redundant annotation fields (~1,028 columns depending on OncoKB evidence). First line: `#SMART_VERSION 1.0.1` |
| `Final_result_tier2.tsv` | TSV | Selected fields for bioinformaticians (~670 columns). Two-row header (field names + source/version metadata) |
| `Final_result_tier3.tsv` | TSV | Clinically relevant fields (~77 columns). Two-row header |
| `variant_counts.txt` | TXT | Per-sample variant counts at each pipeline stage |
| `AnnotatedVcf/` | VCF.gz | VEP-annotated VCFs (if `--keep-tmp`) |

---

## 7. Exact Command Lines Used

The following is the exact Docker command used to process the samples in the most recent run:

```bash
docker run --rm \
  -v /path/to/data:/data \
  -v /Volumes/ExternalSSD/refs:/refs:ro \
  monkiky/smart:1.0.1 \
  <ONCOKB_TOKEN> \
  --transcripts-file /data/transcripts.txt \
  --config /data/Config.yaml \
  --ref-dir /refs \
  --input-dir /data \
  --no-liftover \
  --keep-tables
```

**Parameter explanation:**

| Flag | Value | Meaning |
|------|-------|---------|
| `--no-liftover` | ON | Input VCFs are already in GRCh38 — liftover step skipped |
| `--pass` | ON (default) | Only PASS-filtered variants retained |
| `--keep-tables` | ON | Per-sample intermediate MAF tables kept after merging |
| `--clean-tmp` | ON (default) | Temporary VCF files removed after run |
| `--jobs` | 1 (default) | Samples processed sequentially |
| `--transcripts-file` | transcripts.txt | NM_ accession whitelist for transcript prioritisation |
| `--config` | Config.yaml | Field tier definitions (controls which columns appear in each output tier) |

**Pipeline configuration recorded at runtime:**

Every SMART run prints the following configuration summary to the log before processing begins:

```
============================================================
SMART — Somatic Mutation Annotation and Reporting Tool  v1.0.1
============================================================
  PASS filter:      ENABLED
  Liftover:         DISABLED
  VEP only:         NO
  Clean tmp:        ENABLED
  Clean tables:     ENABLED
  Parallel jobs:    1
  Transcript file:  /data/transcripts.txt
  Config file:      /data/Config.yaml
  Reference dir:    /refs
  Input dir:        /data
============================================================
```

The per-sample log files (written to `output/logs/<sample>.log`) record the complete VEP command, OncoKB API call details, variant counts at each stage, and any warnings. These logs are retained with the output data.

---

## 8. Reproducibility

The entire pipeline runs inside a single Docker container (`monkiky/smart:1.0.1`). The exact image digest is:

```
docker pull monkiky/smart:1.0.1
```

To verify the version of any SMART output file:
```bash
head -1 Final_result_tier1.maf
# → #SMART_VERSION 1.0.1
```

All pipeline parameters are controlled via command-line flags. A representative run command:

```bash
docker run --rm \
  -v /path/to/data:/data \
  -v /path/to/refs:/refs:ro \
  monkiky/smart:1.0.1 \
  <ONCOKB_TOKEN> \
  --transcripts-file /data/transcripts.txt \
  --config /data/Config.yaml \
  --ref-dir /refs \
  --input-dir /data \
  --no-liftover
```

---

*Document generated: May 2026*  
*SMART v1.0.1 · University Hospital Southampton*
