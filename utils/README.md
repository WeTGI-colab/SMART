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
