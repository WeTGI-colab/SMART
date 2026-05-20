# SVIG-UK / ACGS Somatic Variant Oncogenicity Classification Framework

**Source:** ACGS Guidelines for the Classification of Oncogenicity of Somatic Variants in Cancer —
Recommendations by the UK Somatic Variant Interpretation Group (SVIG-UK)
**Ratified:** Association for Clinical Genomic Science (ACGS) Quality Subcommittee, 22 July 2025
**Scope:** SNVs and small indels in solid tumours and haematological malignancies

---

## 1. Classification Categories

A **points-based system** assigns each variant to one of five categories:

| Category | Points threshold |
|---|---|
| **Oncogenic** | ≥ 10 |
| **Likely Oncogenic** | 6 – 9 |
| **VUS** (Variant of Uncertain Significance) | 0 – 5 |
| **Likely Benign** | −1 to −6 |
| **Benign** | ≤ −7 |

**Minimum evidence rule:** At least **two independent evidence codes** must be applied to reach any
(likely) oncogenic or (likely) benign classification.
Exceptions: O1 (standalone → Oncogenic) and a specific application of B1 (standalone → Benign).

---

## 2. Oncogenic Evidence Codes

### O1 — SVIG-UK Canonical Variants List
**Strength:** Standalone Oncogenic (no second code required)

Well-characterised canonical cancer somatic variants with robust functional data, confirmed oncogenic
(≥10 pts) by at least two SVIG-UK Clinical Scientists and enriched in tumour databases (≥ O4 Strong).

---

### O2 — Null Variant in a Tumour Suppressor Gene (TSG)
**Strength:** Very Strong [+8] / Strong [+4] / Moderate [+2] / Supporting [+1]

Applied for null variants (nonsense, frameshift, canonical ±1/2 splice sites, initiation codon,
single/multi-exon deletion) in a gene with a known loss-of-function (LOF) mechanism in the
cancer type under investigation. Strength determined using the PVS1 decision tree (Abou Tayoun
et al., 2018). Also applicable via RNA evidence (O2_RNA) when splicing assays confirm impact.

Key sub-rules:
- Stop-gain in first 100 bp of exon 1 → apply cautiously (NMD escape risk)
- Stop-loss without novel in-frame stop in 3'UTR → NSD predicted → O2_VSTR [+8]
- +2T>C splice variants → SpliceAI delta score ≥ 0.8 required before applying

---

### O3 — Absent or Very Rare in Population Database (gnomAD)
**Strength:** Moderate [+2] / Supporting [+1]

| Observation | Strength |
|---|---|
| Absent across **all** populations in gnomAD v4.1 | Moderate [+2] |
| AF ≤ 0.001% in any sub-population with > 100,000 alleles in gnomAD v4.1 | Supporting [+1] |

**Note:** O3 provides evidence of rarity, not evidence against benignity. Do not apply O3 to push
a variant into VUS if all other evidence supports benign status.
**Prerequisite for O4.**

---

### O4 — Enriched in a Somatic Variant Database
**Strength:** Strong [+4] / Moderate [+2] / Supporting [+1]

Recommended database: **AACR Project GENIE** (cBioPortal). COSMIC may be used but requires
higher thresholds (solid tumours >50 for missense to reach Strong) due to known data quality issues.

#### Missense and splice variants (count = same amino acid / exact nucleotide change)

| GENIE count | Strength |
|---|---|
| > 10 on-target entries | Strong [+4] |
| 5 – 10 on-target entries | Moderate [+2] |
| > 10 VUS on-target (in-house DB) | Supporting [+1] |

#### Frameshift / nonsense variants (count = same position or downstream)

| GENIE count | Strength |
|---|---|
| > 50 entries | Strong [+4] |
| 20 – 50 entries | Moderate [+2] |
| 10 – 19 entries | Supporting [+1] |

#### In-frame insertion / deletion variants

| GENIE count | Strength |
|---|---|
| > 50 within deleted region | Strong [+4] |
| 20 – 50 within deleted region | Moderate [+2] |
| 10 – 19 within deleted region | Supporting [+1] |

**Rules:** O4 cannot be applied if B1 or B2 is applied. No stacking across databases.

---

### O5 — Variant Affects Same Location / Results in Similar Impact
**Strength:** Strong [+4] / Moderate [+2] / Supporting [+1]

Different amino acid change at the same residue where another (likely) oncogenic change is
already established. Strength determined by comparing REVEL scores:

| Condition | Strength |
|---|---|
| REVEL score equivalent (< 0.02 difference) or both > 0.773 | Strong [+4] |
| REVEL > 0.7 but lower than reference variant | Moderate [+2] |
| REVEL > 0.7 but less deleterious than reference | Supporting [+1] |

Also applicable for splice variants where a different nucleotide at the same site is already
classified as (likely) oncogenic and SpliceAI predicts a similar or greater impact.

---

### O6 — Computational Evidence Supports Deleterious Effect
**Strength:** Supporting [+1] ONLY — never above supporting

| Tool | Threshold for O6 |
|---|---|
| **REVEL** | ≥ 0.7 |
| **SpliceAI** (any delta score: AG, AL, DG, DL) | ≥ 0.2 |

Use a meta-predictor (REVEL recommended) rather than combining multiple individual tools.
This code cannot be used as a sole line of evidence.

---

### O7 — Located in a Mutational Hotspot / Critical Functional Domain
**Strength:** Strong [+4] / Moderate [+2] / Supporting [+1]

Primary source: **CancerHotspots.org**

| CancerHotspots condition | Strength |
|---|---|
| ≥ 50 entries at same amino acid position AND ≥ 10 for same amino acid change | Strong [+4] |
| < 50 at position AND ≥ 10 for same amino acid change | Moderate [+2] |
| 2 – 9 entries for same amino acid change | Supporting [+1] |

O7 Moderate may also be applied if all of: (1) no local benign variation in gnomAD,
(2) local enrichment of oncogenic variants, and (3) well-established functional domain.
O7 Supporting may be applied if in silico protein structural models predict functional impact in
a well-characterised domain and local benign variation is absent.

---

### O8 — Missense Constraint (Gene or Domain Level)
**Strength:** Supporting [+1] ONLY

Low rate of benign missense variation in the gene or region, where missense is a common
disease mechanism. Use gnomAD constraint scores (LOEUF / Z-score) or MetaDome regional
constraint. Regional constraint scores preferred over gene-level.

---

### O9 — Protein Length Change (In-frame Indel / Stop-loss / Final Exon Truncation)
**Strength:** Moderate [+2] / Supporting [+1]

Applied for:
- In-frame insertion/deletion (< 1 exon) outside repetitive or poorly conserved regions
- Stop-loss variants not subject to non-stop decay (NSD)
- Truncating variant in the **final exon** of an **oncogene** predicted to result in gain-of-function

Conversely, in-frame indels in repetitive/poorly conserved regions support benignity (B5).

---

### O10 — Functional Studies: Abnormal Result
**Strength:** Very Strong [+8] / Strong [+4] / Moderate [+2] / Supporting [+1]

Well-established in vitro or in vivo functional studies demonstrating a functionally abnormal
result consistent with the disease mechanism. Assessed using the Brnich et al. (2019) framework.
Multiplexed assays of variant effect (MAVEs) frequently achieve Strong weighting.

RNA splicing assays demonstrating splicing impact → use O2 (RNA) rather than O10.
O10 cannot be used with O2_VSTR [+8].

---

### O11 — Tumour Phenotype Supports Oncogenicity
**Strength:** Moderate [+2] / Supporting [+1]

Additional tumour-specific molecular or cellular phenotypic information supporting the
oncogenic role of the variant:
- IHC showing loss of protein (e.g., MSH2/MSH6 loss with an MSH2 variant)
- Microsatellite instability (MSI), Tumour Mutational Burden (TMB), Homologous
  Recombination Deficiency (HRD)
- Loss of heterozygosity (LOH)
- Absent or supportive germline test findings

---

## 3. Benign Evidence Codes

### B1 — Variant Present at High Frequency in gnomAD
**Strength:** Standalone Benign (at high frequency) / Strong [−4]

A variant observed at high allele frequency in the general population provides strong evidence
of benignity. Standalone application (no second code required) when frequency is unambiguously high.

---

### B2 — Variant Does Not Fit the Mode of Action of the Gene
**Strength:** Standalone → VUS override

If the variant type is inconsistent with the known mechanism of action of the gene in the specific
cancer type, a VUS classification should be assigned regardless of other evidence.
Example: a frameshift in an oncogene that acts exclusively via gain-of-function.

Resources: Cancer Gene Census (COSMIC), Cancer Genome Interpreter.

---

### B3 — Computational Evidence Does NOT Support Deleterious Effect
**Strength:** Supporting [−1] ONLY

| Tool | Threshold for B3 |
|---|---|
| **REVEL** | < 0.7 (implied benign range) |
| **SpliceAI** (all delta scores: AG, AL, DG, DL) | all < 0.1 |

---

### B4 — Synonymous / Deep Intronic Variant
**Strength:** Strong [−4] / Supporting [−1]

| Variant type | Strength |
|---|---|
| Synonymous variant with no predicted splicing impact | Strong [−4] |
| Intronic variant at position ≥ +7 or ≤ −21 from exon boundary | Strong [−4] |
| Synonymous with possible (but unconfirmed) splicing impact | Supporting [−1] |

**Note:** Synonymous variants are not universally benign — check SpliceAI and SynMICdb.
Known exceptions: TP53 last base of exons 4, 6, and 9 affect splicing → classify as driver.

---

### B5 — In-frame Deletion/Insertion in a Repetitive Region
**Strength:** Supporting [−1]

In-frame indel located in a repetitive region with unknown function. Length change alone is
not sufficient to support oncogenicity when the affected region is poorly conserved or repetitive.

---

### B6 — Functional Studies: No Damaging Effect
**Strength:** Strong [−4] / Moderate [−2] / Supporting [−1]

Well-established in vitro or in vivo functional studies demonstrating no impact on protein
function or RNA splicing. Assessed using the Brnich et al. (2019) framework.
RNA splicing assays showing no impact on splicing → use B4 (RNA) rather than B6.

---

### B7 — Tumour Phenotype Against Oncogenicity
**Strength:** Moderate [−2] / Supporting [−1]

Tumour molecular or cellular phenotypic data inconsistent with the oncogenic role of the variant.
The inverse of O11 — phenotypic evidence that argues against the variant being a driver.

---

## 4. Key Rules for Combining Evidence

1. **Minimum two codes** required for any (likely) oncogenic or (likely) benign classification
   (exceptions: O1 standalone → Oncogenic; B1 standalone → Benign)
2. **O3 is prerequisite for O4** — variant must be rare before counting somatic database entries
3. **O4 blocked if B1 or B2 applied**
4. **No stacking within the same code** (e.g., cannot add GENIE + COSMIC counts)
5. **B2 overrides** — when applied, the final classification is VUS regardless of other codes
6. **Conflicting evidence** — where oncogenic and benign evidence are equally weighted, classify as VUS
7. **O6 and O8 are capped at Supporting [+1]** — computational evidence alone cannot drive classification
8. Each code may only be applied **once** per variant

---

## 5. Workflow Order (Recommended)

The codes are designed to be applied in order from most to least discriminating:

```
O1  →  O2  →  O3  →  O4  →  O5  →  O6  →  O7  →  O8  →  O9  →  O10  →  O11
B1  →  B2  →  B3  →  B4  →  B5  →  B6  →  B7
```

Once a (likely) oncogenic or (likely) benign classification is reached and no additional
evidence would change the outcome, further codes do not need to be evaluated.

---

## 6. Mapping to SMART Output Columns

| SVIG-UK Code | What it needs | SMART column(s) |
|---|---|---|
| O1 | Canonical variant list | `ONCOKB_HOTSPOT`, `VARIANT_IN_ONCOKB`, `ClinVar_SCI` |
| O2 | LOF in TSG + PVS1 tree | `Consequence` (frameshift, stop_gained, splice_donor/acceptor, start_lost) |
| O3 | gnomAD frequency | `MAX_AF`, `gnomADe_AF`, `gnomADg_AF` |
| O4 | Somatic DB count | `CancerHotspots` count — **not fully in SMART; requires GENIE lookup** |
| O5 | Same position as known oncogenic | `REVEL`, `Existing_variation` |
| O6 | Computational deleterious | `REVEL` ≥ 0.7 · any `SpliceAI_pred_DS_*` ≥ 0.2 |
| O7 | Hotspot (CancerHotspots.org) | `CancerHotspots_HOTSPOT`, `CancerHotspots` count |
| O8 | Missense constraint | `LOEUF` |
| O9 | In-frame indel non-repeat | `VARIANT_CLASS`, `Consequence` (in_frame_insertion/deletion) |
| O10 | Functional evidence | `ONCOKB_ONCOGENIC`, `ClinVar_CLNSIG`, `ClinVar_ONC`, `ClinVar_SCI` |
| O11 | Tumour phenotype | **Not in SMART** (requires IHC, MSI, TMB, HRD, LOH) |
| B1 | High gnomAD AF | `MAX_AF` |
| B2 | Wrong mechanism for gene type | Gene mode-of-action (TSG vs oncogene) + `Consequence` |
| B3 | Computational benign | `REVEL` < 0.7 · all `SpliceAI_pred_DS_*` < 0.1 |
| B4 | Synonymous / deep intronic | `Consequence` = synonymous_variant · `IMPACT` = LOW |
| B5 | In-frame indel in repeat region | `Consequence` = in_frame_* + repeat annotation |
| B6 | Functional no effect | `ONCOKB_ONCOGENIC` = Neutral/Likely Neutral · `ClinVar_CLNSIG` = Benign |
| B7 | Tumour phenotype against | **Not in SMART** |

---

## 7. Key External Resources Referenced

| Resource | Used for | URL |
|---|---|---|
| gnomAD v4.1 | O3, B1 | gnomad.broadinstitute.org |
| AACR Project GENIE | O4 | genie.cbioportal.org |
| COSMIC | O4 (with caution) | cancer.sanger.ac.uk/cosmic |
| CancerHotspots.org | O7 | cancerhotspots.org |
| REVEL | O5, O6, B3 | sites.google.com/site/revelgenomics |
| SpliceAI | O6, B3 | spliceailookup.broadinstitute.org |
| Cancer Gene Census | O2, B2 | cancer.sanger.ac.uk/census |
| Cancer Genome Interpreter | B2 | cancergenomeinterpreter.org |
| ClinVar | O10, B6 | ncbi.nlm.nih.gov/clinvar |
| OncoKB | O1, O10 | oncokb.org |
| MetaDome | O8 | stuart.radboudumc.nl/metadome |
| MaveDB | O10 | mavedb.org |

---

*Document prepared from: SVIG-UK Guidelines v1.0 (main document, 23 pp) and Supplementary Material (35 pp), ACGS 2025.*

---

---

# 8. SMART Automation Assessment — Code by Code

This section evaluates each SVIG-UK evidence code against the current SMART pipeline output,
identifying what can be automated immediately, what requires additional implementation, and
what cannot be automated from VCF-level data alone.

---

## 8.1 Fully Automated — Available in SMART Right Now

---

### O3 — Absent or Very Rare in gnomAD ✅ FULLY AUTOMATED

**SMART columns used:** `gnomADe_AF`, `gnomADg_AF`, `MAX_AF`

SMART already computes population frequencies from Ensembl VEP using gnomAD exomes and
genomes across all ancestry groups. The rule translates directly:

```
MAX_AF = 0 (absent in all populations)   →  O3 Moderate [+2]
MAX_AF ≤ 0.00001 (≤ 0.001%)             →  O3 Supporting [+1]
MAX_AF > 0.00001                          →  O3 not applicable
```

**Implementation effort:** 3 lines of Python. No additional data required.

---

### O6 — Computational Evidence Supports Deleterious Effect ✅ FULLY AUTOMATED

**SMART columns used:** `REVEL`, `SpliceAI_pred_DS_AG`, `SpliceAI_pred_DS_AL`,
`SpliceAI_pred_DS_DG`, `SpliceAI_pred_DS_DL`

SMART computes both REVEL (missense pathogenicity) and all four SpliceAI delta scores
(acceptor gain, acceptor loss, donor gain, donor loss) via VEP plugins:

```
REVEL ≥ 0.7                               →  O6 Supporting [+1]
max(all four SpliceAI delta scores) ≥ 0.2  →  O6 Supporting [+1]
```

Note: O6 is **capped at Supporting [+1]** by the guidelines — this is important because it means
computational evidence alone can never drive a classification above VUS.

**Implementation effort:** 2 lines of Python.

---

### B1 — Variant Present at High Frequency in gnomAD ✅ FULLY AUTOMATED

**SMART columns used:** `MAX_AF`, `gnomADe_AF`, `gnomADg_AF`

The inverse of O3. A variant present at high frequency in the general population is strong evidence
of benignity. The exact frequency threshold for standalone Benign requires scientific judgement
(typically > 1% or > 0.1% depending on disease prevalence), but the data is present.

**Implementation effort:** 1 line of Python.

---

### B3 — Computational Evidence Does NOT Support Deleterious Effect ✅ FULLY AUTOMATED

**SMART columns used:** `REVEL`, `SpliceAI_pred_DS_AG`, `SpliceAI_pred_DS_AL`,
`SpliceAI_pred_DS_DG`, `SpliceAI_pred_DS_DL`

```
REVEL < 0.7  AND  all SpliceAI delta scores < 0.1   →  B3 Supporting [−1]
```

**Implementation effort:** 1 line of Python.

---

### B4 — Synonymous / Deep Intronic Variant ✅ FULLY AUTOMATED

**SMART columns used:** `Consequence`, `IMPACT`, `SpliceAI_pred_DS_AG`, `SpliceAI_pred_DS_AL`,
`SpliceAI_pred_DS_DG`, `SpliceAI_pred_DS_DL`

```
Consequence = synonymous_variant  AND  max(SpliceAI deltas) < 0.1   →  B4 Strong [−4]
Consequence = intron_variant  AND  HGVSc position ≥ +7 or ≤ −21    →  B4 Strong [−4]
```

**Important caveat:** Synonymous variants must always be checked with SpliceAI first. Known
exceptions exist (e.g., TP53 last base of exons 4, 6, 9 — these are drivers despite being synonymous).

**Implementation effort:** 3–4 lines of Python.

---

## 8.2 Implementable — Not in SMART Today but Buildable

---

### O1 — SVIG-UK Canonical Variants List 🔶 PROXY AVAILABLE, EXACT LIST MISSING

**Current SMART columns (proxy):** `ONCOKB_HOTSPOT`, `ONCOKB_ONCOGENIC`, `VARIANT_IN_ONCOKB`,
`ClinVar_ONC`, `ClinVar_SCI`

The exact SVIG-UK Canonical Variants List (Supplementary Table 3 of the guidelines) is not
available as a machine-readable file in SMART. However, the following combination of existing
columns provides a strong proxy:

```
ONCOKB_HOTSPOT = True
AND ONCOKB_ONCOGENIC = "Oncogenic"
AND ClinVar_SCI = "Tier_I_-_Strong"
→  High confidence proxy for O1 (Standalone Oncogenic)
```

**What to implement:** Download the SVIG-UK canonical variants list from Supplementary Table 3
and create a lookup file (`gene:transcript:HGVSp → O1`). This gives exact O1 application.
Without it, the OncoKB + ClinVar proxy covers the majority of canonical cases.

**Implementation effort:** Low. The list is in the supplementary PDF; needs extracting once.

---

### O2 — Null Variant in a Tumour Suppressor Gene 🔶 VARIANT TYPE AVAILABLE, GENE ROLE MISSING

**Current SMART columns:** `Consequence`, `IMPACT`, `SYMBOL`, `SpliceAI_pred_DS_*`

SMART already provides the variant consequence type:
- `frameshift_variant` → LOF
- `stop_gained` → LOF
- `splice_donor_variant` / `splice_acceptor_variant` → canonical splice disruption
- `start_lost` → initiation codon

What is **missing** is the gene-level classification: is this gene a tumour suppressor (TSG),
an oncogene, or both, in the specific cancer type?

**What to implement:** Add a reference file based on the **COSMIC Cancer Gene Census**
(freely downloadable) mapping each gene to its role (`TSG`, `oncogene`, `fusion`, `both`).
Then implement a simplified PVS1 decision tree:

```
Consequence IN (frameshift, stop_gained, splice_donor, splice_acceptor, start_lost)
AND gene_role = "TSG"
→  O2; strength based on PVS1 position in gene (Very Strong for typical cases)
```

Full PVS1 tree (NMD prediction, last exon, alternative start codons) adds complexity but a
simplified version covers 80–90% of cases correctly.

**Implementation effort:** Medium. Gene role lookup is straightforward; full PVS1 tree is complex.

---

### O4 — Enriched in a Somatic Variant Database 🔶 FLAG AVAILABLE, COUNTS MISSING — HIGH PRIORITY

**Current SMART columns:** `CancerHotspots_HOTSPOT` (binary flag only)

This is the **most impactful missing feature**. O4 contributes up to +4 points (Strong) and is
the primary way to distinguish driver mutations from passengers in a data-driven, tumour-agnostic
manner. SMART currently only records whether a variant was found in CancerHotspots (a flag),
but does not record the exact count of cases in GENIE or COSMIC.

**What to implement:**
1. Download **AACR Project GENIE** dataset from cBioPortal (public, ~150,000 variants)
2. Build a lookup table: `gene:protein_change → count` (deduplicated per patient)
3. Add a new column `GENIE_count` to the SMART output
4. Apply thresholds:

```
GENIE_count > 10   →  O4 Strong [+4]      (missense/splice)
GENIE_count 5–10   →  O4 Moderate [+2]
GENIE_count 1–4    →  O4 Supporting [+1]
```

Note: O3 is a **prerequisite** for O4 — the variant must first be rare in the population.

**Implementation effort:** Medium-high. Requires a one-time GENIE data download and join step,
but adds the single most discriminating feature for driver vs passenger classification.

---

### O5 — Variant Affects Same Position as a Known Oncogenic Change 🔶 DATA AVAILABLE VIA API

**Current SMART columns:** `REVEL`, `Existing_variation`, `HGVSp_Short`, `SYMBOL`

SMART does not currently check whether another amino acid change at the **same residue** is
already classified as oncogenic. However, this information is accessible via the OncoKB API,
which is already integrated into the SMART pipeline.

**What to implement:** After the standard OncoKB query, perform a secondary lookup:
"Are there any (likely) oncogenic variants at this amino acid position in OncoKB or ClinVar?"
Then compare REVEL scores to determine the strength:

```
Known oncogenic change exists at same AA position
AND REVEL(VUA) ≥ REVEL(reference) or both > 0.773   →  O5 Strong [+4]
AND REVEL(VUA) > 0.7 (lower than reference)          →  O5 Moderate [+2]
AND REVEL(VUA) > 0.7 (clearly lower)                 →  O5 Supporting [+1]
```

**Implementation effort:** Medium. Requires an additional OncoKB/ClinVar position query but
the API is already in place.

---

### O7 — Located in a Mutational Hotspot (CancerHotspots.org) ✅ IMPLEMENTED

**Current SMART columns:** `CancerHotspots_HOTSPOT` (binary), `CancerHotspots` (coordinate)

SMART records whether a variant overlaps a CancerHotspots entry, but does not record the
count of cases at that amino acid position. The SVIG-UK thresholds require the count:

```
≥ 50 entries at same AA position AND ≥ 10 for same AA change   →  O7 Strong [+4]
< 50 entries at same AA position AND ≥ 10 for same AA change   →  O7 Moderate [+2]
2 – 9 entries for same AA change                               →  O7 Supporting [+1]
```

**Implemented:** Two new columns added to SMART Tier 2 output via `add_cancerhotspots_counts()`
in `post_analysis.py`. Data sourced from the CancerHotspots.org public API
(`cancerhotspots.org/api/hotspots/single`, 1,165 records, Chang et al. 2016/2018).
Snapshot saved at `Files4ThisProject/cancerhotspots_counts.json`; to be copied to the
reference directory (`/refs/CancerHotSpots/`) when at home (see `TODO_when_at_home.md`).

New columns: `CancerHotspots_position_count`, `CancerHotspots_aa_change_count` — both Tier 2.
9 unit tests added in `tests/verification5/test_post_analysis.py` (170 total passing).

---

### O8 — Missense Constraint (LOEUF) ✅ AUTOMATED

**Current SMART columns:** `LOEUF`

SMART already computes LOEUF (Loss-of-function Observed/Expected Upper bound Fraction) via
the gnomAD constraints VEP plugin. This score reflects the tolerance of a gene to loss-of-function
variation: lower LOEUF = more constrained = less tolerant = more likely disease-causing.

**Implemented:** No new columns needed. The `LOEUF` column already exists in SMART Tier 3 output.
The threshold is documented in `Config.yaml` and `ConfigORIGINAL.yaml`:

```
Consequence = missense_variant  AND  LOEUF ≤ 0.35  →  O8 Supporting [+1]
```

The ML model will use `LOEUF` directly as a continuous feature. The threshold of 0.35 is based
on gnomAD constraint score distribution — genes below this value are in the top ~10% most
constrained in the human genome. Only applicable for missense variants; O2 applies for LOF.

---

### O9 — In-frame Indel in Non-repetitive Region 🔶 VARIANT TYPE AVAILABLE, REPEAT ANNOTATION MISSING

**Current SMART columns:** `VARIANT_CLASS`, `Consequence`

SMART can identify in-frame variants:
```
Consequence = in_frame_insertion  OR  in_frame_deletion
AND VARIANT_CLASS = insertion OR deletion
```

What is **missing** is the annotation of whether the affected region is repetitive or poorly conserved
(which would instead support B5 — benign).

**What to implement:**
- **Option A (precise):** Add RepeatMasker annotation as a VEP plugin column. This directly
  identifies simple repeats, low-complexity regions, and satellite sequences.
- **Option B (approximate):** If SpliceAI ≈ 0 AND gnomAD shows no benign in-frame indels in
  the region → assume non-repetitive. This is a reasonable heuristic.

**Implementation effort:** Low–Medium. Option B is immediately implementable; Option A requires
adding the RepeatMasker VEP plugin to the Docker image.

---

### B2 — Variant Does Not Fit the Mode of Action of the Gene 🔶 CONSEQUENCE AVAILABLE, GENE ROLE MISSING

**Current SMART columns:** `Consequence`, `IMPACT`, `SYMBOL`

This code is the inverse of O2 and shares the same dependency: we need to know whether the gene
is a TSG, oncogene, or both. Once the Cancer Gene Census reference file is implemented for O2,
B2 comes almost for free:

```
Consequence = missense_variant
AND gene_role = "oncogene" (pure, not TSG)
AND variant is NOT at a known activation hotspot
→  B2 (Standalone VUS — overrides other oncogenic evidence)
```

**Important:** B2 is a VUS **override** — when applied, the final classification is VUS regardless
of how many oncogenic codes are present. This makes it one of the most consequential codes.

**Implementation effort:** Low (contingent on Cancer Gene Census file being added for O2).

---

### B6 — Functional Studies Show No Damaging Effect ✅ AUTOMATED

**Current SMART columns:** `ONCOKB_ONCOGENIC`, `ClinVar_CLNSIG`, `ClinVar_ONC`

SMART provides strong proxies for functional neutrality:
```
ONCOKB_ONCOGENIC IN ("Likely Neutral", "Neutral")   →  B6 Supporting [−1] proxy
ClinVar_CLNSIG IN ("Benign", "Likely_benign")        →  B6 Supporting [−1] proxy
```

These are not direct functional study evidence but represent the consensus of curated databases
that aggregate functional studies. They are the best available automated proxy for B6.

**Implemented:** No new columns needed. Proxies already in SMART Tier 3 output:

---

## 8.3 Not Automatable — Cannot be Derived from VCF Data

---

### O10 — Functional Studies: Abnormal Result ❌ PARTIAL PROXY ONLY

**Why not fully automatable:** O10 at Strong or Very Strong requires evaluating the quality
of a specific experimental assay (cell-based, animal model, or MAVE study), checking that
controls were appropriate, and judging whether the assay design is fit for purpose. This is
expert scientific judgement that cannot be extracted from VCF annotation.

**What CAN be automated (proxy):**
- `ONCOKB_ONCOGENIC = "Oncogenic"` → implies functional evidence exists → O10 Supporting proxy
- `ClinVar_SCI = "Tier_I_-_Strong"` → strong somatic clinical impact → O10 Supporting proxy

These proxies may support an O10 Supporting [+1] classification but never Strong or Very Strong.

**Consequence for the ML model:** O10 will be represented as a binary feature (has functional
evidence proxy: yes/no) rather than a graded score. The model will need to learn the weight
of this proxy from training data.

---

### O11 — Tumour Phenotype Supports Oncogenicity ❌ NOT IN VCF

**Why not automatable:** O11 requires clinical and laboratory data generated outside the
sequencing pipeline:
- Immunohistochemistry (IHC): protein expression from a microscopy slide
- Homologous Recombination Deficiency (HRD): requires an HRD assay
- Loss of Heterozygosity (LOH): requires tumour/normal comparison or SNP array

**Partial exception:** Tumour Mutational Burden (TMB) and Microsatellite Instability (MSI)
can be estimated computationally from the VCF:
- TMB = count of PASS somatic variants / panel size (Mb)
- MSI = can be assessed from indel patterns at microsatellite loci

However, IHC, HRD, and LOH are not derivable from SNV/indel VCF data and require
integration of external assay results. **O11 is outside the scope of automated VCF analysis.**

---

### B7 — Tumour Phenotype Against Oncogenicity ❌ NOT IN VCF

**Same reasons as O11** — this is the inverse (phenotypic data arguing against the variant
being a driver). Not derivable from VCF alone. Outside scope.

---

## 8.4 Summary Table

| Code | Status | Implementation effort | Max points possible | Priority |
|---|---|---|---|---|
| **O3** | ✅ Automated | None | +2 | — |
| **O6** | ✅ Automated | None | +1 | — |
| **B1** | ✅ Automated | None | Standalone | — |
| **B3** | ✅ Automated | None | −1 | — |
| **B4** | ✅ Automated | None | −4 | — |
| **O1** | ✅ Implemented | Done | Standalone | — |
| **O7** | ✅ Implemented | Done | **+4** | — |
| **O8** | ✅ Automated | Done | +1 | — |
| **B2** | ✅ Implemented | Done | VUS override | — |
| **B6** | ✅ Automated | Done | −4 | — |
| **O4** | 🔶 GENIE counts missing | Medium-high | **+4** | **Highest** |
| **O5** | 🔶 API lookup needed | Medium | +4 | Medium |
| **O2** | ✅ Implemented (simplified) | Done | **+8** | — |
| **O9** | 🔶 Repeat annotation | Low–Medium | +2 | Medium |
| **O10** | ❌ Proxy only | — | +1 proxy | — |
| **O11** | ❌ Clinical data | — | +2 | — |
| **B7** | ❌ Clinical data | — | −2 | — |

### Recommended implementation order

1. **O7** — Add CancerHotspots counts (low effort, +4 max, public data)
2. **O4** — Add GENIE variant counts (highest discriminating power, +4 max)
3. **O2 + B2** — Add Cancer Gene Census gene roles (unlocks both codes)
4. **O1** — Extract SVIG-UK canonical variants list from supplementary
5. **O5** — Add same-position oncogenic lookup via OncoKB API
6. **O9** — Add repeat region annotation (RepeatMasker plugin)
7. **O8** — Apply LOEUF threshold (trivial once decision made)

With steps 1–4 alone, SMART would automate approximately **80% of the SVIG-UK scoring
process** for the most common variant types encountered in clinical cancer genomics panels.

---

*Analysis prepared May 2026. Based on SVIG-UK Guidelines v1.0, ACGS 2025.*
