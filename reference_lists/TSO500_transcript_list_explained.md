# TSO500 Transcript List — Explained

## What TSO500 officially is

The **Illumina TruSight Oncology 500 (TSO500)** panel covers **523 genes**. Despite its
name, it is not 500 genes — the number is the product series name. The panel detects
SNVs, indels, CNAs, fusions, and TMB/MSI across those 523 genes using a targeted
capture approach.

There is **no single official "TSO500 transcript list"** — Illumina defines the panel
by genomic regions (a BED file of capture targets), not by transcript IDs. The
transcript list received from the Birmingham clinical team is their local
interpretation: for each captured gene, which transcript should be used for annotation.

---

## Why 522 transcripts for 515 genes?

The file `TSO500_transcript_MANE.txt` contains **522 transcripts covering 515 unique
genes**. The extra entries come from two sources:

### 1. CDKN2A has two transcripts (intentional)

| Transcript | Isoform | MANE status |
|---|---|---|
| `NM_000077.5` | CDKN2A p16/INK4A | MANE Select |
| `NM_058195.4` | CDKN2A p14ARF | MANE Plus Clinical |

CDKN2A is biologically unique — it encodes **two completely different proteins** from
the same genomic locus using different promoters and reading frames (alternative
splicing of exon 1). Both are clinically relevant tumour suppressors:

- **p16/INK4A** inhibits CDK4/6, regulating the Rb pathway
- **p14ARF** stabilises TP53 via MDM2 inhibition

Both transcripts are included deliberately. This is not an error.

### 2. Four transcripts not found in MANE (gaps)

| TSO500 transcript | Gene | MANE Select (v1.5) |
|---|---|---|
| `NM_001014431.2` | AKT1 | `NM_001382430.1` |
| `NM_001127500.3` | MET | `NM_000245.4` |
| `NM_001130442.3` | HRAS | `NM_005343.4` |
| `NM_001174067.2` | FGFR1 | `NM_023110.3` |

These four TSO500 transcripts are **alternative isoforms not designated as MANE Select
or MANE Plus Clinical**. This is clinically significant:

- Using a non-MANE transcript can produce a **different amino acid position** for the
  same genomic variant, which changes the protein change string sent to OncoKB.
- In the worst case this produces a **false-positive** result: e.g. for FGFR1, the
  non-MANE transcript `NM_001174067.2` maps the same variant to p.Asn577Lys
  (Likely Oncogenic, LEVEL_4 in OncoKB), while the MANE Select `NM_023110.3` gives
  p.Asn546Lys (Unknown in OncoKB). The therapeutic signal is an artefact of
  transcript choice, not biology.
- **MET** is a special case: `NM_001127500.3` (non-MANE) is in the TSO500 list, but
  the MANE Select `NM_000245.4` is also present as a separate entry. The pipeline
  will therefore correctly prioritise the MANE transcript for MET.

---

## MANE coverage summary

| Category | Count |
|---|---|
| Total transcripts in list | 522 |
| Unique genes | 515 |
| Transcripts matching MANE Select | majority |
| Transcripts matching MANE Plus Clinical | minority (e.g. CDKN2A p14ARF) |
| Transcripts not in MANE at all | 4 (AKT1, MET, HRAS, FGFR1) |

Analysis performed against **MANE release 1.5**
(`MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt`),
downloaded from NCBI FTP on 2026-04-22.