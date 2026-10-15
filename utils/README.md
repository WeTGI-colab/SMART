# SMART — Utility Scripts

This directory contains standalone helper scripts that run **outside** the Docker
container: to set up the pipeline environment, to build the optional reference tables
SMART uses when they are present, to generate inputs for a run, and to produce the
comparisons and metrics that are not part of one. None of these scripts require the
Docker image to be built first.

| Script | What it is for |
|---|---|
| [`get_ref_files.sh`](#get_ref_filessh--download-reference-files) | Download and prepare every required reference dataset |
| [`civic_formating.py`](#civic_formatingpy--convert-civic-tsv-to-vcf) | CIViC nightly TSV to VCF (called by `get_ref_files.sh`) |
| [`cosmic_cmc_to_vcf.py`](#cosmic_cmc_to_vcfpy--convert-cosmic-cancer-mutation-census-to-vcf) | Optional COSMIC CMC annotation track |
| [`get_oncokb_transcripts.py`](#get_oncokb_transcriptspy--fetch-oncokb-canonical-transcripts) | OncoKB's own canonical transcripts as a whitelist |
| [`build_svig_uk_resources.py`](#build_svig_uk_resourcespy--cosmic-tables-for-b2-and-o4) | COSMIC gene roles and O4 case counts (Tier 4) |
| [`build_genie_o4.py`](#build_genie_o4py--genie-unique-patient-counts-for-o4) | GENIE unique-patient counts for O4 (Tier 4) |
| [`build_functional_assays.py`](#build_functional_assayspy--mavedb-measurements-for-o10--b6) | MaveDB functional-assay evidence for O10/B6 (Tier 4) |
| [`build_o5_references.py`](#build_o5_referencespy--extra-oncogenic-missense-references-for-o5) | Extra oncogenic missense references for O5 (Tier 4) |
| [`build_gnomad_chip_ab.py`](#build_gnomad_chip_abpy--gnomad-allele-balance-for-the-o3-chip-correction) | gnomAD allele balance for the O3 CHIP correction (Tier 4) |
| [`build_dbnsfp.py`](#build_dbnsfppy--in-silico-predictor-table-for-vep) | Optional dbNSFP predictor table for VEP |
| [`jacklist_to_vcf.py`](#jacklist_to_vcfpy--curated-variant-list-to-vcf) | Curated expert-classified variant list to VCF |
| [`spliceai_compare.py`](#spliceai_comparepy--local-spliceai-against-the-precomputed-tables) | Local SpliceAI vs the precomputed tables |
| [`freeze_spliceai_lookup.py`](#freeze_spliceai_lookuppy--record-the-broads-spliceai-lookup-answers) | Freeze the hosted SpliceAI-lookup answers for the tests |
| [`run_azurify.py`](#run_azurifypy--run-the-azurify-ml-classifier-on-tier-3-output) | Run the Azurify ML classifier on Tier 3 output |
| [`collect_github_traffic.py`](#collect_github_trafficpy--keep-a-permanent-history-of-github-traffic) | Keep a permanent history of GitHub traffic |
| [`gen_fields_doc.py`](#gen_fields_docpy--generate-the-public-field-reference) | Render `docs/fields.md` from `Config.yaml` |

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

### Verification and the reference manifest (STEP 13)

The last step of the build runs `scripts/verify_reference_files.py generate`, which inspects
the files it has just built and records `smart_reference_manifest.json` in the install root
(`<base_dir>/refs/`). That manifest is what every later reference check consumes —
`verify_reference_files.py verify` before annotating, the entrypoint's pre-run check, and
`tests/verification5/test_verify_reference_files.py`.

The generator refuses to certify a file that fails its structural checks, so a non-zero exit
from this step means the reference directory is not fit to annotate against. (A structurally
corrupt CIViC VCF survived six months in the lab reference directory precisely because the
build never looked at its own output; the manifest exists to stop that recurring.)

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

---

# Optional reference tables

These builders produce reference tables SMART uses when they are present and does without
when they are not. All of them read data that is licence-gated, large, or both, which is why
SMART cannot ship the results: each one is built once from the user's own download and kept
with the rest of the reference data. **None of their outputs may be committed or
redistributed** unless the section below says otherwise.

The five SVIG-UK builders feed the Tier 4 oncogenicity classifier
(`scripts/svig_uk_classifier.py`); without them the affected evidence codes are simply not
awarded, so the classification is weaker rather than wrong.

## `build_svig_uk_resources.py` — COSMIC tables for B2 and O4

Builds the two licence-gated lookup tables the classifier needs, from the user's own COSMIC
downloads, into `--out-dir`:

- `cgc_gene_roles.tsv` — gene mode of action from the COSMIC Cancer Gene Census
  (`ROLE_IN_CANCER`), the resource SVIG-UK names for mode of action (B2);
- `cosmic_o4_counts.tsv.gz` — somatic case counts for evidence code O4, aggregated to the
  protein level as SVIG-UK requires (missense counted per amino-acid change, regardless of the
  underlying nucleotide change). The count is `COSMIC_SAMPLE_MUTATED` across all tissues: the
  CMC `DISEASE` breakdown is empty for ~99.5 % of mutations (v104), so tissue-specific counting
  is not possible from this file.

```bash
python3 utils/build_svig_uk_resources.py \
    --cgc-tar /refs/COSMIC/Cosmic_CancerGeneCensus_Tsv_v104_GRCh38.tar \
    --cmc-tar /refs/COSMIC/CancerMutationCensus_AllData_Tsv_v104_GRCh37.tar \
    --out-dir /refs/SVIG_UK
```

The GRCh37 CMC release is fine here, because O4 counts are protein-level.

## `build_genie_o4.py` — GENIE unique-patient counts for O4

Builds the AACR Project GENIE case-count table the classifier prefers over COSMIC for O4
("enriched in a somatic variant database ... international (e.g. GENIE)"), passed as
`--genie-o4`.

For every mutation in `data_mutations_extended.txt` the builder records the sample's OncoTree
code, tissue and lineage (so counting can be limited to on-target tumour types at the level
`--o4-on-target` selects), lifts GENIE's GRCh37 coordinates to GRCh38 with the UCSC chain file
SMART already uses — checking the reference allele against the GRCh38 FASTA and left-aligning
indels, so variants match by genomic change rather than by GENIE's transcript numbering — and
counts **unique patients** per variant, as O4 requires. A sidecar `<out>.cohort.tsv` records
the unique patients per level, the denominators shown beside the counts in the output columns.

```bash
python3 utils/build_genie_o4.py \
    --genie-dir /refs/GENIE/release --oncotree /refs/GENIE/oncotree.json \
    --chain /refs/liftover/hg19ToHg38.over.chain.gz --fasta /refs/genome/hg38.fa \
    --out /refs/GENIE/genie_o4_counts.tsv.gz
```

GENIE is released under a data-use agreement that forbids redistribution.

## `build_functional_assays.py` — MaveDB measurements for O10 / B6

Builds the functional-assay reference for O10 (and B6) from MaveDB. SVIG-UK's operative
criterion is that a variant is functionally abnormal when protein activity is **<25 % of
wild-type**; the *strength* is not fixed by the guideline (it follows ClinGen SVI / CanVIG-UK
calibration of the individual assay), so it is recorded per assay with its justification.

Unlike the other builders, this one writes into the repository, because the assay it
currently carries is CC0 (the `licence` column records this per assay) —
`reference_lists/svig_uk_functional_assays.tsv` (one row per assay:
thresholds, strength, provenance) and `reference_lists/svig_uk_functional_assay_scores.tsv`
(one row per variant: assay, protein change, % of wild type).

```bash
python3 utils/build_functional_assays.py            # refresh from MaveDB
python3 utils/build_functional_assays.py --check    # verify the stored files, no network
```

## `build_o5_references.py` — extra oncogenic missense references for O5

Builds an optional table of (likely) oncogenic missense reference variants for the SVIG-UK O5
missense route, passed as `--o5-references`. The SVIG-UK canonical variants list is always
used by the classifier; this table adds the Cancer Genome Interpreter (CGI) catalogue of
validated oncogenic mutations, read from the gene:protein dictionary distributed with
[OncoVI](https://github.com/MGCarta/oncovi) (`resources/cgi_dictionary.txt`).

CGI entries carry no transcript, so the classifier uses them only when the residue number and
reference amino acid agree with the transcript SMART annotated.

```bash
python3 utils/build_o5_references.py \
    --cgi-dictionary /path/to/oncovi/resources/cgi_dictionary.txt \
    --out /refs/SVIG_UK/o5_references.tsv
```

Output columns: `gene transcript protein assessment source`. CGI data are subject to CGI's
terms of use.

## `build_gnomad_chip_ab.py` — gnomAD allele balance for the O3 CHIP correction

Builds the gnomAD v4.1 allele-balance table used for the SVIG-UK O3 CHIP correction
(`--chip-allele-balance`): in CHIP-associated genes, gnomAD's allele-balance histogram can be
used to exclude low-VAF alleles and re-calculate the frequency from alleles with VAF between
0.3 and 0.65.

For every coding RefSeq transcript of the requested genes, coding exon parts ±25 bp are
queried **remotely** (tabix over HTTPS) in the gnomAD v4.1 exome and genome sites VCFs; only
variants with `AC > 0` are written. Requires `tabix` built with HTTPS support and internet
access to `storage.googleapis.com`.

```bash
python3 utils/build_gnomad_chip_ab.py \
    --transcript-models /refs/ncbiRefSeq/ncbiRefSeq.txt.gz \
    --genes DNMT3A ASXL1 TET2 TP53 \
    --out /refs/SVIG_UK/gnomad_v4.1_chip_ab.tsv.gz
```

`--genes` defaults to the four genes named in the guideline, which adds "but not limited to".
Output columns: `chrom pos ref alt dataset filter AC AN nhomalt ab_hist_alt_bin_freq`.

## `build_dbnsfp.py` — in silico predictor table for VEP

Builds the tabix-indexed dbNSFP file SMART passes to VEP's dbNSFP plugin (optional, like
CADD). dbNSFP academic releases are CC BY-NC-ND 4.0 and must be downloaded by each user from
<https://www.dbnsfp.org>; the built file must not be redistributed.

Only predictors SMART does not already annotate are kept, one entry per tool, so no score is
repeated: SIFT, SIFT4G, PolyPhen-2, REVEL and CADD are excluded because SMART annotates them
directly, and where dbNSFP ships several flavours of a tool one is chosen (VARITY_R, BayesDel
without allele frequency, Eigen, GERP++, phyloP/phastCons on the 100-way vertebrate
alignment). Rank scores are left out — they are rescaled copies of the same score.

```bash
python3 utils/build_dbnsfp.py \
    --dbnsfp-dir /path/to/dbNSFP5.0a \
    --out <base_dir>/refs/dbNSFP/dbNSFP5.0a_grch38.smart.gz
```

`--dbnsfp-dir` is the directory holding `dbNSFP5.0a_variant.chr*.gz` from the release zip.

---

# Development, comparison and reporting tools

These are not part of a pipeline run. They exist to produce evidence for a decision, to turn
curated data into pipeline input, or to keep metrics.

## `jacklist_to_vcf.py` — curated variant list to VCF

Converts a curated, expert-classified variant list into a bgzip-compressed, tabix-indexed VCF
ready for SMART, for validating the pipeline's classification against a human one. Two export
formats are recognised automatically: the original header-less 18-column `jack_list.txt`, and
a header-based CSV/TSV export (`VariantChr`, `VariantPosition`, ... `VariantStatus`).

Alongside the VCF it writes a sidecar `<out>.classifications.tsv` holding the human
classification per variant, keyed by the same `CHROM:POS:REF:ALT` identifier the SMART output
uses, plus curation provenance (transcript, dates, analyst, reviewer, status, source file).
Several exports can be given at once: they are concatenated and deduplicated on
`CHROM:POS:REF:ALT`, keeping the most recently reviewed row, and variants curated twice with
different classifications are reported so they can be resolved before use.

```bash
python3 utils/jacklist_to_vcf.py \
    --input jack2.csv \
    --output data/jack2_validation/jack2.vcf.gz \
    --transcripts data/jack2_validation/transcripts.txt \
    --tumour-type Myeloid
```

Also takes `--status Curated` (keep only signed-off rows), `--since DD/MM/YYYY`, and
`--keep-duplicates`. The tumour type and the transcript change the result (the tumour type
affects O4); everything else is provenance.

## `spliceai_compare.py` — local SpliceAI against the precomputed tables

Scores variants with a locally run SpliceAI and compares the result with Illumina's
precomputed score files, which SMART looks variants up in. Those tables carry no transcript —
one score per gene, computed on GENCODE v24 canonical — and were precomputed within a ±50 bp
window. Running the model locally removes both limits; the question this script answers is
whether that changes anything that matters, and at what cost: how many variants change score,
how many cross 0.2 (the delta from which the classifier awards O6 and `ANNOTATION_FLAGS`
keeps a LOW/MODIFIER row), how many cross 0.5, where the strongest signal lies beyond 50 bp,
and wall-clock per variant.

No model code is reimplemented: the `spliceai` CLI reads a VCF and writes a VCF, with the gene
annotation built by `scripts/build_spliceai_annotation.py` from the panel transcript list.
Because SpliceAI reports the annotation's `NAME` column as its `SYMBOL` field, and that
annotation names rows `GENE--REFSEQ`, the output states which transcript each score belongs
to — the thing the precomputed tables never recorded.

```bash
python3 utils/spliceai_compare.py \
    --maf path/to/Final_result_tier1.maf \
    --annotation .../tso500_grch38.spliceai_annotation.txt \
    --reference .../hg38.fa \
    --precomputed .../SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz \
    --distance 500 --out comparison.tsv
```

`--spliceai` points at the executable if it is not on `PATH` (a venv, typically).

## `freeze_spliceai_lookup.py` — record the Broad's SpliceAI-lookup answers

Freezes the hosted SpliceAI-lookup service's answers for a small set of **public** variants,
so the test suite can compare against an independent implementation of the same model without
calling it on every run (the service is rate limited, asks that batch work be done locally,
and returned HTTP 500 on four of 93 rows during this work). The answers are recorded with the
date and the endpoint; if the frozen values later stop matching, either SMART's annotation
changed or the service's Gencode release did, and both are worth knowing.

```bash
python3 utils/freeze_spliceai_lookup.py tests/verification5/data/spliceai_lookup_frozen.json
```

Run it deliberately, not in CI — queries are paced at 8 s. Nothing patient-derived is sent to
a third party. The frozen file is consumed by
`tests/verification5/test_spliceai_lookup_agreement.py`.

## `run_azurify.py` — run the Azurify ML classifier on Tier 3 output

Runs Azurify (faryabiLab/Azurify, bioRxiv 2025), a CatBoost model trained on >15,000
clinically classified variants, over SMART's `Final_result_tier3.tsv`: converts the TSV to
Azurify's input format, runs Azurify (snpEff annotation plus its internal feature lookup),
merges the predictions back and writes `Final_result_tier3_ml.tsv`.

```bash
python3 utils/run_azurify.py \
    --smart-output /path/to/Output_Results \
    --azurify-dir  /path/to/Azurify \
    --snpeff-jar   /path/to/snpEff/snpEff.jar \
    --output-dir   /path/to/Output_Results/azurify
```

Installed outside Docker (`pip install catboost liftover tqdm pandas`, plus a clone of
Azurify) and needs >16 GB RAM. SNVs and small indels only: CNAs and structural variants are
skipped and preserved in the merged output with empty ML columns.

## `collect_github_traffic.py` — keep a permanent history of GitHub traffic

GitHub's Insights → Traffic keeps only 14 days and offers no download counter for a
repository. This script snapshots what GitHub does expose, merges it into CSV files under
`metrics/` that grow over time, and adds the counters that *are* cumulative (Docker Hub pulls,
release asset downloads, stars/forks).

```bash
python3 utils/collect_github_traffic.py              # collect and merge
python3 utils/collect_github_traffic.py --summary    # print what is stored, no network
```

Needs the `gh` CLI authenticated with push access (the traffic endpoints are not public). Run
it at least every 14 days, or days fall off the end and are lost. Daily rows
(`traffic_clones.csv`, `traffic_views.csv`) are merged on the date, because a day is only
final once it has passed; the snapshot files (`referrers.csv`, `paths.csv`, `totals.csv`)
append.

## `gen_fields_doc.py` — generate the public field reference

Renders `docs/fields.md`, the field reference the README links as "Full field reference",
from `Config.yaml`, so that the configuration is the single thing to edit. Output is
deterministic — no timestamps, nothing run-dependent — which is what makes `--check`
meaningful and keeps a regeneration from producing a spurious diff.

```bash
python3 utils/gen_fields_doc.py            # rewrite docs/fields.md
python3 utils/gen_fields_doc.py --check    # exit 1 if docs/fields.md is stale
```

`--config` and `--out` override the input configuration and the page written.
`tests/verification5/test_fields_doc.py` fails if the two have drifted apart.
