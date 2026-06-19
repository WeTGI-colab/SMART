# SMART — Utility Scripts

This directory contains standalone helper scripts that run **outside** the Docker
container, either to set up the pipeline environment or to generate inputs for it.
None of these scripts require the Docker image to be built first.

---

## `get_ref_files.sh` — Download reference files

Downloads and prepares every reference dataset the SMART pipeline needs.
Run this once before building or running the Docker container.

### What it installs

| Resource | Source | Notes |
|----------|--------|-------|
| GRCh38 reference genome | UCSC | FASTA + `.fai` + `.dict` |
| hg19 → hg38 liftover chain | UCSC | Required for coordinate conversion |
| ClinVar | NCBI FTP | bgzip VCF + tabix index |
| CIViC | CIViC nightly | Converted to VCF by `civic_formating.py` |
| REVEL | MSSM | bgzip TSV + tabix index |
| SpliceAI | Illumina BaseSpace | **Manual download required** (see below) |
| gnomAD constraints (LOEUF) | Google Drive | bgzip TSV + tabix index |
| VEP cache | Ensembl FTP | `homo_sapiens_merged_vep_114_GRCh38` (~15 GB) |
| VEP plugins | GitHub | Clones `Ensembl/VEP_plugins` |
| Cancer Hotspots | GitHub | bgzip VCF + tabix index |
| CADD | Kircher lab | **Optional, opt-in** (`INSTALL_CADD=1`) — very large (~80 GB) |

### Requirements

- `wget` / `curl`, `git`, `gunzip`, `tabix`, `apptainer` on `PATH`
- Python 3.8+ (for the CIViC conversion step)
- At least **200 GB** of free disk space
- Internet access

### Usage

```bash
bash utils/get_ref_files.sh /path/to/install/root
# e.g.
bash utils/get_ref_files.sh /Volumes/ExternalSSD
```

Reference files are installed under `<base_dir>/refs/`. A `pipeline_config.sh`
file is written to the current directory recording all resolved paths and
resource versions for reproducibility.

The script is **idempotent** — existing files are detected and skipped, so it
is safe to re-run if a download is interrupted.

### SpliceAI (manual step)

SpliceAI requires a free Illumina BaseSpace account and cannot be downloaded
automatically. Before running the script, log into BaseSpace and download:

- SNV scores: `spliceai_scores.raw.snv.hg38.vcf.gz`
- INDEL scores: `spliceai_scores.raw.indel.hg38.vcf.gz`

Place both files in `<base_dir>/spliceai_staging/` then run the script.

### CADD (optional, opt-in)

CADD adds genome-wide deleteriousness scores (`CADD_PHRED` / `CADD_RAW`) covering
all variant types, not just missense (unlike REVEL). The GRCh38 score files are
**very large** (whole-genome SNVs ~80 GB + indels), so they are **not** downloaded
by default. To install them, re-run the downloader with the opt-in flag:

```bash
INSTALL_CADD=1 bash utils/get_ref_files.sh /path/to/install/root
```

This populates `<base_dir>/refs/CADD/` with `whole_genome_SNVs.tsv.gz` and
`gnomad.genomes.r4.0.indel.tsv.gz` (+ tabix indexes) from CADD v1.7. CADD is
optional: if the directory is present the entrypoint adds the `--plugin CADD`
track automatically; if absent, annotation runs unchanged. Source and licence:
<https://cadd.bihealth.org/> (free for non-commercial use).

---

## `civic_formating.py` — Convert CIViC TSV to VCF

Called automatically by `get_ref_files.sh`. Can also be run standalone when
CIViC releases a new nightly build and the VCF needs refreshing.

Reads the CIViC `nightly-VariantSummaries.tsv` and writes a bgzip-compressed,
tabix-indexed VCF suitable for use as a VEP `--custom` annotation track.
Uses only the Python standard library — no third-party packages required.

### Usage

```bash
python3 utils/civic_formating.py \
    --input  nightly-VariantSummaries.tsv \
    --output /path/to/refs/CIVIC/civic_grch38.vcf.gz \
    --assembly grch38
```

| Option | Default | Description |
|--------|---------|-------------|
| `--input` | required | CIViC nightly TSV file |
| `--output` | required | Output `.vcf.gz` path |
| `--assembly` | `grch38` | Genome build: `grch38` or `grch37` |
| `--chain-file` | — | Chain file for liftover (used by `get_ref_files.sh`) |
| `--bgzip` | `bgzip` | Path or wrapper for bgzip binary |
| `--tabix` | `tabix` | Path or wrapper for tabix binary |

---

## `cosmic_cmc_to_vcf.py` — Convert COSMIC Cancer Mutation Census to VCF

Builds the **optional** COSMIC reference track for SMART. COSMIC data is licence-gated
and **not redistributable**, so SMART cannot ship it — each user downloads the Cancer
Mutation Census (CMC) with their own (free for non-commercial research) COSMIC account
and runs this script to produce the VCF the pipeline consumes. No COSMIC data is bundled
with SMART; only this transformation code.

If the resulting VCF is present at `<refs>/COSMIC/cosmic_cmc_grch38.vcf.gz`, the
entrypoint automatically adds it as a VEP `--custom` track, populating `COSMIC_CNT`
(recurrence), `COSMIC_TIER` (significance tier 1/2/3), `COSMIC_ONC_TSG` (gene role) and
related columns. If the file is absent, annotation runs unchanged — COSMIC is skipped.

### Get the input

1. Register / log in at <https://cancer.sanger.ac.uk/cosmic>.
2. Downloads → project **Cancer Mutation Census** → **All Data CMC** (GRCh37 — the CMC is
   only released on GRCh37, but the file carries a GRCh38 coordinate column, which this
   script uses, so the output is GRCh38-native).
3. You get `CancerMutationCensus_AllData_Tsv_v<NN>_GRCh37.tar`.

### Usage

```bash
python3 utils/cosmic_cmc_to_vcf.py \
    --input  CancerMutationCensus_AllData_Tsv_v104_GRCh37.tar \
    --output /path/to/refs/COSMIC/cosmic_cmc_grch38.vcf.gz
```

The `.tar`, the inner `.tsv.gz`, or a plain `.tsv` are all accepted as `--input`.
Requires `bgzip`, `tabix` and `sort` on `PATH` (uses only the Python standard library).

| Option | Default | Description |
|--------|---------|-------------|
| `--input` | required | CMC `.tar` / `.tsv.gz` / `.tsv` |
| `--output` | required | Output `.vcf.gz` path |
| `--bgzip` / `--tabix` / `--sort` | from PATH | Override binary locations |

### Scope

Only substitutions where REF and ALT are equal-length A/C/G/T strings (SNVs and MNVs) are
emitted — they match VEP `exact` by position+allele with no reference genome needed. True
indels (≈9 % of the CMC) are skipped (they would need left-normalisation against the
genome FASTA) and reported in the run summary. **Do not commit or redistribute the
resulting VCF** — it is derived COSMIC data (see <https://www.cosmickb.org/terms/>).

---

## `get_oncokb_transcripts.py` — Fetch OncoKB canonical transcripts

Queries the OncoKB API (`/utils/allCuratedGenes`) and writes a transcript
whitelist file containing the GRCh38 RefSeq NM accession that OncoKB uses
internally for each curated gene. The output can be passed directly to the
SMART pipeline via `--transcripts-file`.

### Why this matters

OncoKB annotates variants based on a specific transcript per gene. If your
preferred-transcript whitelist uses a different isoform, the protein change
sent to OncoKB may not match the one it expects, leading to incorrect or
missing evidence levels. Using OncoKB's own transcripts as the whitelist
guarantees that the protein change the pipeline extracts from VEP is the same
one OncoKB will recognise.

### Usage

```bash
python3 utils/get_oncokb_transcripts.py \
    --token  "$ONCOKB_TOKEN" \
    --output oncokb_transcripts.txt \
    --tsv    oncokb_transcripts_summary.tsv
```

| Option | Default | Description |
|--------|---------|-------------|
| `--token` | required | OncoKB API token |
| `--output` | `oncokb_transcripts.txt` | One NM accession per line; pass to `--transcripts-file` |
| `--tsv` | — | Optional summary TSV with gene, NM accession, and Ensembl transcript columns |
| `--no-version` | off | Strip version suffix (e.g. `NM_005228` instead of `NM_005228.3`) |

### Output

`oncokb_transcripts.txt` — one NM accession per line, ready for use:

```
NM_005228.3
NM_000546.5
NM_004333.4
...
```

`oncokb_transcripts_summary.tsv` (optional) — tab-separated with three columns:

```
HugoSymbol    grch38RefSeq     grch38Isoform
EGFR          NM_005228.3      ENST00000275493
TP53          NM_000546.5      ENST00000269305
BRAF          NM_004333.4      ENST00000646891
...
```

As of the current OncoKB release, **1008 genes** are curated; **981** have a
GRCh38 RefSeq transcript. The 27 genes without one are mostly non-coding loci
(T-cell receptor segments `TRA`/`TRB`/`TRD`/`TRG`) or catch-all entries
(`Other Biomarkers`). Notable gap: `NTRK3` is curated for clinical evidence
but currently has no GRCh38 RefSeq assigned.

### Note on transcript versions

OncoKB's transcript versions lag slightly behind the current NCBI releases
(e.g. `NM_005228.3` for EGFR vs the current MANE Select `NM_005228.5`).
The SMART pipeline compares NM accessions **without** version suffixes, so
this version difference does not affect matching.
