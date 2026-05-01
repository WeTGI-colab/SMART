# Field Configuration System

## Overview

This directory contains a YAML-based field metadata and tier classification system used to document, filter, and prioritise columns in the variant annotation pipeline output. The system is split into two files:

- `ConfigORIGINAL.yaml` — the human-maintained template defining field metadata and wildcard patterns
- `Config.yaml` — the generated output consumed by downstream tools, produced by running `build_config.py`

---

## File Descriptions

### `ConfigORIGINAL.yaml` — The Template

This is the source of truth for field definitions. It contains two top-level sections:

#### `fields`

A dictionary where each key is a column name as it appears in the output TSV, and each value is a metadata block with three required properties:

```yaml
fields:

  HGVSp:
    description: The HGVS protein sequence name
    source: VEP
    version: 114.2
```

| Property | Description |
|---|---|
| `description` | Human-readable explanation of the field, including links to documentation where relevant |
| `source` | The tool, database, or standard that produces this field (e.g. VEP, ClinVar, OncoKB, Manta) |
| `version` | The version of the source used during annotation |

#### `field_patterns`

A list of wildcard pattern entries used for fields that are generated dynamically in groups — specifically the expanded OncoKB diagnostic and treatment implication columns. Instead of defining hundreds of near-identical entries, a single pattern covers all indexed variants of a field:

```yaml
field_patterns:

  - pattern: "ONCOKB_DIAG_*_levelOfEvidence"
    description: "OncoKB diagnostic implication level of evidence"
    source: OncoKB
    version: "Refer to the ONCOKB_DATA_VERSION column and datetime of the last update of the variant"

  - pattern: "ONCOKB_TX_*_drugs"
    description: "OncoKB treatment drugs (combined as Drug1 + Drug2)"
    source: OncoKB
    version: "Refer to the ONCOKB_DATA_VERSION column and datetime of the last update of the variant"
```

The `*` wildcard matches any numeric index (e.g. `ONCOKB_DIAG_0_`, `ONCOKB_DIAG_1_`, ..., `ONCOKB_TX_29_`). Pattern matching uses Python's `fnmatch` module. At build time, the script expands these patterns against the actual columns present in the data, injecting concrete field entries into the `fields` block before tier assignment.

Currently defined pattern groups are:

- `ONCOKB_DIAG_*` — OncoKB diagnostic implication fields (level of evidence, tumor type, PMIDs, description, etc.)
- `ONCOKB_TX_*` — OncoKB treatment implication fields (drugs, FDA level, associated cancer type, etc.)

---

### `Config.yaml` — The Generated Output

This file is produced by running `build_config.py` and should **not be edited manually**. It is identical in structure to `ConfigORIGINAL.yaml` but with two additions applied to every field:

1. A `tier` property assigned based on the tier classification rules (see below)
2. Updated `description` values for `ClinVar_*` fields, sourced automatically from the ClinVar VCF header

---

## Build Script: `build_config.py`

### Usage

```bash
python3 build_config.py
```

Input and output paths are defined at the top of the script:

```python
INPUT_YAML  = "ConfigORIGINAL.yaml"
OUTPUT_YAML = "Config.yaml"
CLINVAR_VCF = "/path/to/clinvar.vcf.gz"
HOTSPOTS_VCF = "/path/to/hotspots.vcf.gz"
```

### What the script does

1. **Loads** `ConfigORIGINAL.yaml`
2. **Expands** `field_patterns` wildcards against actual column names (requires passing actual columns — see note below)
3. **Assigns tiers** to every field based on the sets defined in the script (see Tier Definitions below)
4. **Updates ClinVar descriptions** by reading the `##INFO` header lines from the ClinVar VCF and replacing placeholder descriptions for `ClinVar_*` fields with the authoritative text from the source file
5. **Writes** the result to `Config.yaml`

---

## Tier Definitions

Every field in `Config.yaml` receives a `tier` value. Tiers are used downstream to control which columns are shown, exported, or filtered.

| Tier | Value in YAML | Meaning |
|---|---|---|
| Tier 1 | `1` | Primary clinical or analytical fields — shown by default |
| Tier 2 | `2` | Supporting or technical fields — available but not shown by default |
| Tier 3 | `3` | Reserved for future use — defined by the team |
| Drop | `"drop"` | Fields to be excluded from output entirely |

Tier assignment priority (highest to lowest):

```
DROP > TIER3 > TIER2 > TIER1 (default)
```

Fields not explicitly listed in `DROP_FIELDS`, `TIER3_FIELDS`, or `TIER2_FIELDS` default to Tier 1.

The sets `DROP_FIELDS`, `TIER2_FIELDS`, and `TIER3_FIELDS` are defined directly in the script and should be updated there when the field list changes. **TODO**

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

Then assign its tier by adding it to the appropriate set in `build_config.py` (`TIER2_FIELDS`, `DROP_FIELDS`, etc.), or leave it unset to default to Tier 1.

### Dynamic fields (indexed groups)

If the new fields follow a numbered pattern (e.g. `MY_TOOL_0_score`, `MY_TOOL_1_score`, ...), add a pattern entry to the `field_patterns` section of `ConfigORIGINAL.yaml`:

```yaml
  - pattern: "MY_TOOL_*_score"
    description: "Score reported by MyTool for implication entry N"
    source: MyTool
    version: mytool_v1.0
```

And update `expand_field_patterns()` in the script to pass the actual column list from your data.

---

## Sources and Versions

| Source | Version used | Fields prefixed with |
|---|---|---|
| VEP (Ensembl Variant Effect Predictor) | 114.2 | *(VEP standard fields)* |
| ClinVar | ClinVar_2024-12 | `ClinVar_` |
| CIViC | CIViC_v3.6 | `CIViC_` |
| OncoKB | See `ONCOKB_DATA_VERSION` column | `ONCOKB_`, `ONCOKB_DIAG_*`, `ONCOKB_TX_*` |
| gnomAD Exomes | gnomAD_v4.0 | `gnomADe_` |
| gnomAD Genomes | gnomAD_v4.0 | `gnomADg_` |
| SpliceAI | 1.3 | `SpliceAI_` |
| REVEL | 1.3 | `REVEL` |
| Cancer Hotspots | changv2 | `CancerHotspots_` |
| Manta (structural variants) | Manta_v1.6 | *(SV-specific fields)* |
| VCF standard | VCF_v4.2 | `FORMAT`, `AD`, `GT`, etc. |

---

## Notes

- OncoKB field versions are not static — they reflect the version of the OncoKB dataset at the time of annotation. Always refer to the `ONCOKB_DATA_VERSION` column in the output TSV for the exact version used per variant.
- `ClinVar_GENEINFO_` is in the drop list because it is a composite field not suitable for direct output.
- `ONCOKB_JSON`, `ONCOKB_treatments`, `ONCOKB_diagnosticImplications`, and `ONCOKB_prognosticImplications` are dropped because their content is expanded into dedicated indexed columns (`ONCOKB_TX_*`, `ONCOKB_DIAG_*`) during post-processing.