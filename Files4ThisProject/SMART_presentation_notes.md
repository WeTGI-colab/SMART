# SMART — Technical Presentation Notes
## Step-by-step walkthrough of what the pipeline does

---

## 0. Startup & environment check

- Docker container starts; all tools (VEP, GATK, Python, OncoKB annotator) are
  pre-installed — no dependency management needed at runtime
- Pipeline version is read from the `VERSION` file and printed in the terminal banner
- External API versions are fetched live: OncoKB API version and data version
  (e.g. v7.0, 30/03/2026)
- Tool versions are printed: VEP, GATK, bcftools, Python
- Reference file versions are printed: ClinVar filename, CIViC, SpliceAI, REVEL, LOEUF,
  OncoKB data date
- Runtime config is printed: PASS filter on/off, liftover on/off, transcript file path
- All input VCF files in the input directory are discovered and listed

---

## Per-sample loop (repeated for every `.vcf.gz` in the input directory)

---

## Step 1 — PASS filtering (`bcftools` / `awk`)

- Reads the raw VCF (may contain thousands of variants including low-quality calls)
- Keeps only variants where the `FILTER` column is exactly `PASS`
- Variant callers differ in how they mark filtered variants:
  - MuTect2/Strelka/Pisces use `PASS`
  - SomaticSniper uses `.` (no PASS) — requires `--no-pass` flag to skip this step
- Output: compressed `*.PASS_only.vcf.gz` + tabix index
- Variant count recorded: Original → PASS

---

## Step 2 — LiftOver (optional, GATK 4.6 `LiftoverVcf`)

- Converts variant coordinates from hg19/GRCh37 → GRCh38
- Uses the UCSC `hg19ToHg38.over.chain` chain file
- Uses the GRCh38 reference FASTA for sequence validation
- Variants that cannot be lifted (unmapped regions, strand flips) are written to a
  separate `*_rejected.vcf.gz` file and excluded from downstream analysis
- `--WRITE_ORIGINAL_POSITION` flag preserves the original hg19 coordinates as INFO
  fields in the output
- Skipped with `--no-liftover` if input is already GRCh38 (e.g. DRAGEN TSO500 output)
- Variant count recorded: PASS → Lifted

---

## Step 3 — VEP annotation (Ensembl VEP 114, offline, merged cache)

- The most computationally intensive step
- Runs fully offline using the pre-downloaded merged VEP cache (~25 GB)
- `--everything` flag enables all standard annotations
- `--merged` cache includes both RefSeq and Ensembl transcripts
- For each variant, VEP returns one CSQ (consequence) entry **per overlapping
  transcript** — a single variant may generate 10–50+ CSQ entries
- Plugins activated:
  - **SpliceAI** (cutoff 0.5): delta scores for splice acceptor/donor gain/loss
  - **REVEL**: missense pathogenicity score (0–1)
  - **LOEUF**: loss-of-function constraint metric from gnomAD, matched by transcript
- Custom annotation sources added via `--custom`:
  - **ClinVar** (exact match): clinical significance, review status, disease, HGVS
  - **CIViC** (exact match): variant type, consequence, clinical evidence
  - **CancerHotspots** (exact match): hotspot flags (2D, 3D, non-coding)
- All annotations are embedded in the VCF INFO/CSQ field — one annotated VCF per sample
- Variant count recorded: Lifted → VEP-annotated

---

## Step 4 — OncoKB annotation (`oncokb2.0.py`)

- Custom Python script — the clinical heart of the pipeline
- Reads the VEP-annotated VCF variant by variant
- For each variant, **transcript prioritisation** selects ONE preferred transcript
  from the VEP CSQ field using a 3-tier logic:
  - **Tier 1**: first transcript whose NM accession matches the preferred whitelist
  - **Tier 2**: transcript tagged as MANE Select or MANE Plus Clinical by VEP
  - **Tier 3**: fallback — first transcript VEP reports
- The selected transcript determines:
  - For **SNVs/indels**: the protein change (HGVSp_Short, e.g. `L858R`) sent to OncoKB
  - For **CNAs**: the gene name (Hugo Symbol) sent to OncoKB
- Routes each variant to the correct OncoKB REST API endpoint:
  - SNVs, indels, MantaINS, MantaBND → `/annotate/mutations/byProteinChange`
  - MantaDUP, MantaDEL, GAIN, LOSS → `/annotate/copyNumberAlterations`
- OncoKB returns per-variant: oncogenicity, mutation effect, all therapeutic levels
  (LEVEL_1 through LEVEL_R2), FDA levels, diagnostic levels, prognostic levels,
  gene and variant summaries, PMIDs
- Results written back into the VCF as INFO fields (`ONCOKB_*`)
- Variant count recorded: VEP → OncoKB-annotated

> **Known limitation (v0.1.0):** when a variant overlaps multiple transcripts in the
> preferred whitelist, only the first match is used. The other isoforms are discarded.
> Most affected gene: CDKN2A (p16 vs p14ARF). Will be fixed in v0.2.0.

---

## Step 5 — VCF to Table (`vcf2table.py`)

- Converts the annotated VCF into a flat CSV table — one row per variant
- Parses the multi-value CSQ field, applying the same 3-tier transcript prioritisation
  as Step 4 to ensure the transcript shown in the table is identical to the one used
  for the OncoKB query
- All VEP fields are expanded into individual columns
- Populates `NM_Transcript` column with the selected transcript ID (version-stripped)
- Output: one CSV per sample in `Table/`

---

## Step 6 — MafAnnotator (OncoKB MAF annotator)

- Runs the official OncoKB `MafAnnotator.py` on the per-sample CSV
- Re-queries OncoKB using the standard MAF format fields (Hugo_Symbol, HGVSp_Short)
- Produces a MAF-format file per sample in `FINAL_Table/`
- Note: MafAnnotator does not correctly handle CNA rows — the pipeline's post-analysis
  step overrides MafAnnotator's CNA output with values from the OncoKB API (Step 4)

---

## Step 7 — Post-analysis (`post_analysis.py`)

Runs once after all samples are processed. Performs:

- **Merge**: combines all per-sample MAF files from `FINAL_Table/` into a single
  merged DataFrame
- **JSON expansion**: OncoKB returns treatment, diagnostic and prognostic implications
  as JSON arrays; these are expanded into individual numbered columns
  (`ONCOKB_TX_0_level`, `ONCOKB_TX_0_drugs`, `ONCOKB_DIAG_0_level`, etc.)
  — the number of columns varies per run depending on how many implications the
  variants carry (this is why per-sample MAFs have different column counts)
- **CNA correction**: overrides MafAnnotator's incorrect CNA annotation with values
  from the API stored in Step 4
- **Tier splitting**: applies the Config.yaml field configuration to split columns into
  three audience-targeted output files
- **MAF header**: writes `#SMART_VERSION x.y.z` as the first line of the MAF file

### Output files

| File | Audience | Description |
|---|---|---|
| `Final_result_tier1.maf` | Downstream tools / bioinformatics | All non-dropped columns, standard MAF format |
| `Final_result_tier2.tsv` | Bioinformaticians | Selected fields, two-row header (names + metadata) |
| `Final_result_tier3.tsv` | Clinical scientists | Clinically relevant fields, two-row header |

---

## Summary — data flow

```
Input VCFs
    │
    ▼
[Step 1] PASS filter         → removes low-quality calls
    │
    ▼
[Step 2] LiftOver (optional) → hg19 → GRCh38
    │
    ▼
[Step 3] VEP annotation      → functional consequence, population freq,
    │                           splice, pathogenicity, ClinVar, CIViC,
    │                           CancerHotspots (all transcripts in CSQ)
    ▼
[Step 4] OncoKB annotation   → transcript prioritisation → clinical evidence
    │                           (oncogenicity, therapeutic levels, summaries)
    ▼
[Step 5] VCF → Table         → flat CSV, same transcript as Step 4
    │
    ▼
[Step 6] MafAnnotator        → MAF format per sample
    │
    ▼
[Step 7] Post-analysis       → merge all samples, expand JSON, tier split
    │
    ▼
Final_result_tier1.maf  /  tier2.tsv  /  tier3.tsv
```

---

## Key design decisions worth highlighting

- **Single transcript per variant**: ensures the transcript used for VEP display,
  the transcript used to build the OncoKB query, and the transcript shown in the
  clinical report are always the same — no inconsistency between annotation layers
- **Offline VEP**: no internet required for VEP; only OncoKB requires live API access
- **Docker**: all tools and versions are fixed inside the image — fully reproducible
  across institutions and machines
- **Version tracking**: `VERSION` file → banner → MAF header; every output file is
  self-describing with the pipeline version that produced it
- **Config.yaml**: controls which fields appear in which output tier; can be updated
  without touching the pipeline code