# Field Configuration System

## Overview

This directory contains a YAML-based field metadata and tier classification system used to document, filter, and prioritise columns in the variant annotation pipeline output. The system is split into two files:

- `ConfigORIGINAL.yaml` — the human-maintained template defining field metadata
- `Config.yaml` — the tier-annotated output consumed by downstream tools, produced by running `update_confg.py`

> The build script is `create_config/update_confg.py`. The misspelling is in the real filename; there is no `build_config.py` in this repository, although the script's own header and parts of its documentation used to call it that.

---

## File Descriptions

### `ConfigORIGINAL.yaml` — The Template

This is the source of truth for field definitions. It contains one top-level section, `fields`.

#### `fields`

A dictionary where each key is a column name as it appears in the output TSV, and each value is a metadata block with three required properties:

```yaml
fields:

  HGVSp:
    description: The HGVS protein sequence name
    source: VEP
    version: 114.0
```

| Property | Description |
|---|---|
| `description` | Human-readable explanation of the field, including links to documentation where relevant |
| `source` | The tool, database, or standard that produces this field (e.g. VEP, ClinVar, OncoKB, Manta) |
| `version` | The version of the source used during annotation |

Note that `description`, `source` and `version` here are a *fallback*. Each run harvests the provenance of every field it can from its own VCF headers (`scripts/vcf_metadata.py`), and that takes precedence; `Config.yaml` supplies the tier, the MAF standard columns, what SMART computes, and the OncoKB MafAnnotator columns, which carry no header of their own.

#### `field_patterns` — in `Config.yaml`, not in the template

A list of wildcard pattern entries used for fields that are generated dynamically in groups — specifically the expanded OncoKB diagnostic and treatment implication columns. Instead of defining hundreds of near-identical entries, a single pattern covers all indexed variants of a field:

```yaml
field_patterns:

  - pattern: "ONCOKB_DIAG_*_levelOfEvidence"
    description: "OncoKB diagnostic implication level of evidence"
    source: OncoKB
    version: "Refer to the ONCOKB_DATA_VERSION column and datetime of the last update of the variant"
    tier: 2

  - pattern: "ONCOKB_TX_*_drugs"
    description: "OncoKB treatment drugs (combined as Drug1 + Drug2)"
    source: OncoKB
    version: "Refer to the ONCOKB_DATA_VERSION column and datetime of the last update of the variant"
    tier: 2
```

The `*` wildcard matches any index (e.g. `ONCOKB_DIAG_0_`, `ONCOKB_DIAG_1_`, ..., `ONCOKB_TX_29_`). **The patterns are never expanded into concrete field entries.** They stay a list of patterns in `Config.yaml` — 37 entries in the root copy — and are matched against the run's actual column names at run time by `fnmatch`: `resolve_field_meta()` in `scripts/post_analysis.py` and the pattern fallback in `scripts/vcf_metadata.py`. `update_confg.py` only assigns each pattern a tier.

Currently defined pattern groups are:

- `ONCOKB_DIAG_*` — OncoKB diagnostic implication fields (level of evidence, tumor type, PMIDs, description, etc.)
- `ONCOKB_TX_*` — OncoKB treatment implication fields (drugs, FDA level, associated cancer type, etc.)

The committed `ConfigORIGINAL.yaml` has drifted from the generated file: it carries no `field_patterns` section and 305 fields, against 427 fields and 37 patterns in the root `Config.yaml`. Re-running `update_confg.py` against the template as committed would therefore *not* reproduce the config the pipeline uses. In practice the root `Config.yaml` is now edited directly, and `update_confg.py` is the tier-assignment tool for when the template is brought back up to date.

---

### `Config.yaml` — The Generated Output

This file is produced by running `update_confg.py`. It is identical in structure to `ConfigORIGINAL.yaml` but with two additions applied to every field:

1. A `tier` property assigned based on the tier classification rules (see below)
2. Updated `description` values for `ClinVar_*` fields, sourced automatically from the ClinVar VCF header

### Where the generated file has to end up

`update_confg.py` writes `Config.yaml` into the directory it is run from (`OUTPUT_YAML = "Config.yaml"`, a relative path), so by default it lands here in `create_config/`. **Nothing in the pipeline reads that copy.** The config is passed in explicitly:

- `entrypoint.sh` takes `--config`, and the documented container invocation passes `--config /data/Config.yaml` — i.e. a copy that sits in the mounted input directory;
- the repository root `Config.yaml` is the authoritative copy. `tests/verification5/test_config_copies_agree.py` compares every other `Config.yaml` in the repository against it, so a correction made only in `create_config/` is a correction the pipeline never sees.

A regenerated config must therefore be copied to the repository root (and to the input directory that is mounted as `/data`) before it has any effect.

---

## Build Script: `update_confg.py`

### Usage

```bash
python3 create_config/update_confg.py        # reads ./ConfigORIGINAL.yaml, writes ./Config.yaml
```

Input and output paths are module-level constants near the top of the script:

```python
INPUT_YAML   = "ConfigORIGINAL.yaml"
OUTPUT_YAML  = "Config.yaml"
CLINVAR_VCF  = "/mnt/data1/vep_test/ClinVar/clinvar.vcf.gz"
HOTSPOTS_VCF = "/mnt/data1/vep_test/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz"
```

Both VCF paths are absolute paths on the machine the script was written on and have to be edited before a run. `HOTSPOTS_VCF` is defined but not read anywhere in the script.

### What the script does

1. **Loads** `ConfigORIGINAL.yaml`
2. **Assigns tiers** to every field, from the sets defined in the script (see Tier Definitions below)
3. **Assigns tiers** to every entry in `field_patterns`, by matching the pattern string against `TIER2_ONCOKB_TX_PATTERNS` and `TIER2_ONCOKB_DIAG_PATTERNS`
4. **Updates ClinVar descriptions** by reading the `##INFO` header lines from the ClinVar VCF and replacing placeholder descriptions for `ClinVar_*` fields with the authoritative text from the source file
5. **Writes** the result to `Config.yaml`

---

## Tier Definitions

Every field in `Config.yaml` receives a `tier` value. The tier decides which output files the column reaches — the routing is in `scripts/post_analysis.py`, where the final column sets are resolved:

| Tier | Value in YAML | Where the column goes |
|---|---|---|
| Tier 3 | `3` | **Both** `Final_result_tier3.tsv` (clinical scientists) **and** `Final_result_tier2.tsv` |
| Tier 2 | `2` | `Final_result_tier2.tsv` (bioinformaticians) only |
| Tier 1 | `1` | `Final_result_tier1.maf` only — **neither TSV** |
| Drop | `"drop"` | Excluded from every output file |

**Tier 3 is the clinical sheet, and Tier 1 is not a sheet at all.** Tier 1 means "kept in the MAF, in no TSV": a field assigned tier 1 is silently absent from both tiered TSVs. To put a column in front of a clinical scientist, give it tier 3 (which also puts it in the Tier 2 sheet, a superset); to give it to a bioinformatician only, tier 2. The sets are named for their audience in `update_confg.py` itself — `# --- Tier 3 (clinician) ---` and `# --- Tier 2 (bioinformatician) ---`.

Tier assignment priority (highest to lowest), as implemented:

```
DROP > TIER3 > TIER2 > ONCOKB_TX_*/ONCOKB_DIAG_* patterns > TIER1 (default)
```

Fields not explicitly listed in `DROP_FIELDS`, `TIER3_FIELDS`, or `TIER2_FIELDS`, and not matching an OncoKB pattern group, default to Tier 1 — and so appear in the MAF only.

The sets `DROP_FIELDS`, `TIER2_FIELDS`, and `TIER3_FIELDS` are defined directly in `update_confg.py` and should be updated there when the field list changes. For reference, the root `Config.yaml` currently assigns 273 fields to tier 2 and 92 to tier 3; a run writes more columns than that, because the `ONCOKB_TX_*`/`ONCOKB_DIAG_*` patterns resolve to one column per indexed entry present in the data.

---

## Adding New Fields

### Static fields (single column)

Add an entry to the `fields` section of `ConfigORIGINAL.yaml`:

```yaml
  MY_NEW_FIELD:
    description: What this field contains
    source: ToolName
    version: tool_v1.0
```

Then assign its tier by adding it to the appropriate set in `update_confg.py` (`TIER3_FIELDS` for the clinical sheet, `TIER2_FIELDS` for the bioinformatician sheet, `DROP_FIELDS` to exclude it). Leaving it unset defaults it to Tier 1, which keeps it out of both TSVs.

### Dynamic fields (indexed groups)

If the new fields follow a numbered pattern (e.g. `MY_TOOL_0_score`, `MY_TOOL_1_score`, ...), add a pattern entry to the `field_patterns` section of `Config.yaml`:

```yaml
  - pattern: "MY_TOOL_*_score"
    description: "Score reported by MyTool for implication entry N"
    source: MyTool
    version: mytool_v1.0
    tier: 2
```

No column list has to be supplied: `fnmatch` resolves the pattern against whatever columns the run produces. Patterns outside the two OncoKB groups are not covered by the script's pattern-tier rules, so give such an entry its `tier` directly.

---

## Sources and Versions

| Source | Version used | Fields prefixed with |
|---|---|---|
| VEP (Ensembl Variant Effect Predictor) | 114.0 | *(VEP standard fields)* |
| ClinVar | ClinVar_2024-12 | `ClinVar_` |
| CIViC | CIViC_v3.6 | `CIViC_` |
| OncoKB | See `ONCOKB_DATA_VERSION` column | `ONCOKB_`, `ONCOKB_DIAG_*`, `ONCOKB_TX_*` |
| gnomAD Exomes | gnomAD_v4.1 | `gnomADe_` |
| gnomAD Genomes | gnomAD_v4.1 | `gnomADg_` |
| SpliceAI | 1.3 | `SpliceAI_` |
| REVEL | 1.3 | `REVEL` |
| Cancer Hotspots | changv2 | `CancerHotspots_` |
| Manta | Manta_v1.6 | *(the SV INFO fields: `SVTYPE`, `SVLEN`, `END`, `MATEID`, `BND_DEPTH`, `CIPOS`/`CIEND`, `HOMLEN`/`HOMSEQ`, ...)* |
| VCF standard | VCF_v4.2 | `FORMAT`, `AD`, `GT`, etc. |

---

## Notes

- OncoKB field versions are not static — they reflect the version of the OncoKB dataset at the time of annotation. Always refer to the `ONCOKB_DATA_VERSION` column in the output TSV for the exact version used per variant.
- The Manta row above is a *fallback attribution*, and names the caller the laboratory's inputs were called with. The SV INFO fields come from whichever caller wrote the input VCF; a run harvests the caller's own `##INFO` descriptions and `##source` from the VCF headers, which is what the output reports. Fields the caller writes but SMART does not add — `hotspot`, `FractionInformativeReads` — are attributed to "Input VCF (variant caller)" for that reason, not to a named tool.
- `ClinVar_GENEINFO_` is in the drop list because it is a composite field not suitable for direct output.
- `ONCOKB_JSON`, `ONCOKB_treatments`, `ONCOKB_diagnosticImplications`, and `ONCOKB_prognosticImplications` are dropped because their content is expanded into dedicated indexed columns (`ONCOKB_TX_*`, `ONCOKB_DIAG_*`) during post-processing.
