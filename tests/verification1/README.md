# verification1

Field-level verification of the SMART pipeline against live external APIs.

## Purpose

After running the pipeline on a known set of variants, this verification checks
that the annotation fields in the output MAF match what the source APIs return
when queried directly for the same variants. It catches regressions in OncoKB
annotation logic, VEP field extraction, and post-processing transformations.

## Input

| File | Description |
|---|---|
| `verification1.vcf.gz` | 18 synthetic variants (14 SNV/indel + 4 CNA) |
| `verification1_transcripts.txt` | 15 preferred NM transcript IDs |
| `Config.yaml` | Pipeline config (tiers, column metadata) |

## Variants

| Gene | Variant ID | Change | Type |
|---|---|---|---|
| NRAS | NRAS_Q61R | p.Q61R | SNV |
| IDH1 | IDH1_R132H | p.R132H | SNV |
| PIK3CA | PIK3CA_H1047R | p.H1047R | SNV |
| EGFR | EGFR_del19 | p.Glu746_Ala750del | Indel |
| EGFR | EGFR_L858R | p.L858R | SNV |
| BRAF | BRAF_V600E | p.V600E | SNV |
| GNAQ | GNAQ_Q209L | p.Gln209= | SNV (synonymous) |
| PTEN | PTEN_R130stop | p.R130* | SNV |
| KRAS | KRAS_G13D | p.G13D | SNV |
| KRAS | KRAS_G12D | p.G12D | SNV |
| BRCA2 | BRCA2_L3101R | p.L3101R | SNV |
| DICER1 | DICER1_E1705K | p.E1705* | SNV |
| TP53 | TP53_R175H | p.R175H | SNV |
| ERBB2 | ERBB2_S310F | p.S310F | SNV |
| MET | MantaDUP:MET_AMP | Amplification | CNA |
| CDKN2A | MantaDEL:CDKN2A_DEL | Deletion | CNA |
| CDK4 | MantaDUP:CDK4_AMP | Amplification | CNA |
| ERBB2 | MantaDUP:ERBB2_AMP | Amplification | CNA |

## Verification script

`verify.py` runs two independent modules against the pipeline MAF output:

**OncoKB module** — queries `oncokb.org/api/v1` for each variant and compares:
- Gene/variant presence (`GENE_IN_ONCOKB`, `VARIANT_IN_ONCOKB`, `ONCOKB_ALLELE_EXIST`)
- Oncogenicity and mutation effect (`ONCOGENIC`, `MUTATION_EFFECT`, `ONCOKB_VUS`, `ONCOKB_HOTSPOT`)
- Treatment levels (`LEVEL_1` through `LEVEL_R2`, `HIGHEST_SENSITIVE_LEVEL`, `HIGHEST_RESISTANCE_LEVEL`)
- FDA level and implication levels (`ONCOKB_highestFdaLevel`, `ONCOKB_DIAG_LVL`, `ONCOKB_PROG_LVL`)
- Summaries (`ONCOKB_GENE_SUMMARY`, `ONCOKB_VARIANT_SUMMARY`, `ONCOKB_diagnosticSummary`, `ONCOKB_prognosticSummary`)

**VEP module** — queries `rest.ensembl.org` for each SNV/indel (CNAs skipped) and compares:
- Consequence, impact, biotype, exon/intron number
- HGVSc, HGVSp, protein position, amino acids, codons
- SIFT, PolyPhen, strand, canonical flag

Fields not verified (require local plugin databases unavailable via REST API):
`gnomAD*`, `ClinVar*`, `SpliceAI*`, `REVEL`, `LOEUF`, `CancerHotspots*`, `CIViC*`

## Usage

```bash
# Build and run the pipeline first (something like)
docker run --rm -v /path/to/SMART/tests/verification1:/data -v /path/to/SMART/tests/verification1/output:/output -v /Volumes/ExternalSSD/refs:/refs:ro monkiky/smart:1.0.0 "$ONCOKB_TOKEN" --transcripts-file /data/verification1_transcripts.txt --ref-dir /refs --config /data/Config.yaml --no-liftover --keep-tmp --input-dir /data --keep-tables

# Run both verification modules
python3 tests/verification1/verify.py \
    --maf  tests/verification1/output/output/Final_result_tier1.maf \
    --token $ONCOKB_TOKEN \
    --output tests/verification1/results.tsv
```

Exit code `0` = all checks passed. Exit code `1` = one or more mismatches or errors.

It will generate a report like this

```
############################################################
MODULE: OncoKB
############################################################

────────────────────────────────────────────────────────────
[1/18]  TUMOR  |  NRAS  |  NRAS_Q61R  |  p.Q61R
  → Mutation query: NRAS Q61R
  ✓ ONCOKB_ALLELE_EXIST                 MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='NRAS, a GTPase, is mutated in a diverse range of cancer'  API='NRAS, a GTPase, is mutated in a diverse range of cancer'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda2'  API='LEVEL_Fda2'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='The NRAS Q61R mutation is known to be oncogenic.'  API='The NRAS Q61R mutation is known to be oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Gain-of-function'  API='Gain-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF=''  API=''  [PASS]
  ✓ LEVEL_2                             MAF='Cobimetinib,Trametinib'  API='Cobimetinib,Trametinib'  [PASS]
  ✓ LEVEL_3A                            MAF='Binimetinib,Cobimetinib,Selumetinib+Iodine I 131-6-Beta'  API='Binimetinib,Cobimetinib,Selumetinib,Iodine I 131-6-Beta'  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF='Binimetinib+Ribociclib'  API='Binimetinib,Ribociclib'  [PASS]
  ✓ LEVEL_R1                            MAF='Cetuximab,Tucatinib+Trastuzumab,Panitumumab'  API='Cetuximab,Tucatinib,Trastuzumab,Panitumumab'  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_2'  API='LEVEL_2'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF='LEVEL_R1'  API='LEVEL_R1'  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF='LEVEL_Dx2'  API='LEVEL_Dx2'  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF='LEVEL_Px1'  API='LEVEL_Px1'  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[2/18]  TUMOR  |  IDH1  |  IDH1_R132H  |  p.R132H
  → Mutation query: IDH1 R132H
  ✓ ONCOKB_ALLELE_EXIST                 MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='IDH1, a cell metabolism enzyme, is recurrently mutated '  API='IDH1, a cell metabolism enzyme, is recurrently mutated '  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda2'  API='LEVEL_Fda2'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='The IDH1 R132H mutation is known to be oncogenic.'  API='The IDH1 R132H mutation is known to be oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Switch-of-function'  API='Switch-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF='Ivosidenib,Olutasidenib,Vorasidenib'  API='Ivosidenib,Olutasidenib,Vorasidenib'  [PASS]
  ✓ LEVEL_2                             MAF='Ivosidenib'  API='Ivosidenib'  [PASS]
  ✓ LEVEL_3A                            MAF='Ivosidenib'  API='Ivosidenib'  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF=''  API=''  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_1'  API='LEVEL_1'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF=''  API=''  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF='LEVEL_Dx2'  API='LEVEL_Dx2'  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF='LEVEL_Px2'  API='LEVEL_Px2'  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[3/18]  TUMOR  |  PIK3CA  |  PIK3CA_H1047R  |  p.H1047R
  → Mutation query: PIK3CA H1047R
  ✓ ONCOKB_ALLELE_EXIST                 MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='PIK3CA, the catalytic subunit of PI3-kinase, is frequen'  API='PIK3CA, the catalytic subunit of PI3-kinase, is frequen'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda2'  API='LEVEL_Fda2'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='The PIK3CA H1047R mutation is known to be oncogenic.'  API='The PIK3CA H1047R mutation is known to be oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Gain-of-function'  API='Gain-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF='Alpelisib+Fulvestrant,Capivasertib+Fulvestrant,Inavolis'  API='Alpelisib,Fulvestrant,Capivasertib,Inavolisib,Palbocicl'  [PASS]
  ✓ LEVEL_2                             MAF=''  API=''  [PASS]
  ✓ LEVEL_3A                            MAF=''  API=''  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF='RLY-2608,RLY-2608+Fulvestrant'  API='RLY-2608,Fulvestrant'  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_1'  API='LEVEL_1'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF=''  API=''  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[4/18]  TUMOR  |  EGFR  |  EGFR_del19  |  p.Glu746_Ala750del
  → Mutation query: EGFR Glu746_Ala750del
  ✓ ONCOKB_ALLELE_EXIST                 MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='EGFR, a receptor tyrosine kinase, is altered by amplifi'  API='EGFR, a receptor tyrosine kinase, is altered by amplifi'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda2'  API='LEVEL_Fda2'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='The EGFR E746_A750del alteration is known to be oncogen'  API='The EGFR E746_A750del alteration is known to be oncogen'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Gain-of-function'  API='Gain-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF='Afatinib,Amivantamab+Chemotherapy,Amivantamab+Lazertini'  API='Afatinib,Amivantamab,Chemotherapy,Lazertinib,Dacomitini'  [PASS]
  ✓ LEVEL_2                             MAF=''  API=''  [PASS]
  ✓ LEVEL_3A                            MAF='Patritumab Deruxtecan'  API='Patritumab Deruxtecan'  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF=''  API=''  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_1'  API='LEVEL_1'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF=''  API=''  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[5/18]  TUMOR  |  EGFR  |  EGFR_L858R  |  p.L858R
  → Mutation query: EGFR L858R
  ✓ ONCOKB_ALLELE_EXIST                 MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='EGFR, a receptor tyrosine kinase, is altered by amplifi'  API='EGFR, a receptor tyrosine kinase, is altered by amplifi'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda2'  API='LEVEL_Fda2'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='The EGFR L858R mutation is known to be oncogenic.'  API='The EGFR L858R mutation is known to be oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Gain-of-function'  API='Gain-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF='Afatinib,Amivantamab+Chemotherapy,Amivantamab+Lazertini'  API='Afatinib,Amivantamab,Chemotherapy,Lazertinib,Dacomitini'  [PASS]
  ✓ LEVEL_2                             MAF=''  API=''  [PASS]
  ✓ LEVEL_3A                            MAF='Patritumab Deruxtecan'  API='Patritumab Deruxtecan'  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF=''  API=''  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_1'  API='LEVEL_1'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF=''  API=''  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[6/18]  TUMOR  |  MET  |  MantaDUP:MET_AMP  |  
  → CNA query: MET AMPLIFICATION
  ✓ ONCOKB_ALLELE_EXIST                 MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='MET, a receptor tyrosine kinase, is recurrently altered'  API='MET, a receptor tyrosine kinase, is recurrently altered'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda2'  API='LEVEL_Fda2'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='MET amplification is known to be oncogenic.'  API='MET amplification is known to be oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Gain-of-function'  API='Gain-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF=''  API=''  [PASS]
  ✓ LEVEL_2                             MAF='Capmatinib,Crizotinib,Tepotinib'  API='Capmatinib,Crizotinib,Tepotinib'  [PASS]
  ✓ LEVEL_3A                            MAF='Telisotuzumab Vedotin'  API='Telisotuzumab Vedotin'  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF=''  API=''  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF='Erlotinib,Gefitinib,Osimertinib'  API='Erlotinib,Gefitinib,Osimertinib'  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_2'  API='LEVEL_2'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF='LEVEL_R2'  API='LEVEL_R2'  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[7/18]  TUMOR  |  BRAF  |  BRAF_V600E  |  p.V600E
  → Mutation query: BRAF V600E
  ✓ ONCOKB_ALLELE_EXIST                 MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='BRAF, an intracellular kinase, is frequently mutated in'  API='BRAF, an intracellular kinase, is frequently mutated in'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda2'  API='LEVEL_Fda2'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='The BRAF V600E mutation is known to be oncogenic.'  API='The BRAF V600E mutation is known to be oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Gain-of-function'  API='Gain-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF='Dabrafenib,Dabrafenib+Trametinib,Encorafenib+Binimetini'  API='Dabrafenib,Trametinib,Encorafenib,Binimetinib,Cetuximab'  [PASS]
  ✓ LEVEL_2                             MAF='Encorafenib+Panitumumab,Selumetinib,Vemurafenib,Vemuraf'  API='Encorafenib,Panitumumab,Selumetinib,Vemurafenib,Cobimet'  [PASS]
  ✓ LEVEL_3A                            MAF='Vemurafenib,Dabrafenib'  API='Vemurafenib,Dabrafenib'  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF=''  API=''  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_1'  API='LEVEL_1'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF=''  API=''  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF='LEVEL_Dx2'  API='LEVEL_Dx2'  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[8/18]  TUMOR  |  CDKN2A  |  MantaDEL:CDKN2A_DEL  |  
  → CNA query: CDKN2A DELETION
  ✓ ONCOKB_ALLELE_EXIST                 MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='CDKN2A, which encodes both a cyclin-dependent kinase in'  API='CDKN2A, which encodes both a cyclin-dependent kinase in'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda3'  API='LEVEL_Fda3'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='CDKN2A deletion is known to be oncogenic.'  API='CDKN2A deletion is known to be oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Loss-of-function'  API='Loss-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF=''  API=''  [PASS]
  ✓ LEVEL_2                             MAF=''  API=''  [PASS]
  ✓ LEVEL_3A                            MAF=''  API=''  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF='Palbociclib,Ribociclib,Abemaciclib'  API='Palbociclib,Ribociclib,Abemaciclib'  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_4'  API='LEVEL_4'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF=''  API=''  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF='LEVEL_Dx2'  API='LEVEL_Dx2'  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[9/18]  TUMOR  |  GNAQ  |  GNAQ_Q209L  |  p.Gln209=
  → Mutation query: GNAQ Gln209=
  ✓ ONCOKB_ALLELE_EXIST                 MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='GNAQ, a G protein subunit, is recurrently mutated in uv'  API='GNAQ, a G protein subunit, is recurrently mutated in uv'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF=''  API=''  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='This is a synonymous mutation and is not annotated by O'  API='This is a synonymous mutation and is not annotated by O'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='False'  API='False'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Unknown'  API='Unknown'  [PASS]
  ✓ ONCOGENIC                           MAF='Unknown'  API='Unknown'  [PASS]
  ✓ LEVEL_1                             MAF=''  API=''  [PASS]
  ✓ LEVEL_2                             MAF=''  API=''  [PASS]
  ✓ LEVEL_3A                            MAF=''  API=''  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF=''  API=''  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF=''  API=''  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF=''  API=''  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[10/18]  TUMOR  |  PTEN  |  PTEN_R130stop  |  p.R130*
  → Mutation query: PTEN R130*
  ✓ ONCOKB_ALLELE_EXIST                 MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='PTEN, a lipid and protein phosphatase, is one of the mo'  API='PTEN, a lipid and protein phosphatase, is one of the mo'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda2'  API='LEVEL_Fda2'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='The PTEN R130* mutation is known to be oncogenic.'  API='The PTEN R130* mutation is known to be oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Loss-of-function'  API='Loss-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF='Capivasertib+Fulvestrant'  API='Capivasertib,Fulvestrant'  [PASS]
  ✓ LEVEL_2                             MAF=''  API=''  [PASS]
  ✓ LEVEL_3A                            MAF=''  API=''  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF='GSK2636771,AZD8186'  API='GSK2636771,AZD8186'  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_1'  API='LEVEL_1'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF=''  API=''  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF='LEVEL_Dx3'  API='LEVEL_Dx3'  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[11/18]  TUMOR  |  KRAS  |  KRAS_G13D  |  p.G13D
  → Mutation query: KRAS G13D
  ✓ ONCOKB_ALLELE_EXIST                 MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='KRAS, a GTPase which functions as an upstream regulator'  API='KRAS, a GTPase which functions as an upstream regulator'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda2'  API='LEVEL_Fda2'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='The KRAS G13D mutation is known to be oncogenic.'  API='The KRAS G13D mutation is known to be oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Gain-of-function'  API='Gain-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF='Avutometinib+Defactinib'  API='Avutometinib,Defactinib'  [PASS]
  ✓ LEVEL_2                             MAF='Cobimetinib,Trametinib'  API='Cobimetinib,Trametinib'  [PASS]
  ✓ LEVEL_3A                            MAF='Cobimetinib,Trametinib'  API='Cobimetinib,Trametinib'  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF='Trametinib,Cobimetinib,Binimetinib'  API='Trametinib,Cobimetinib,Binimetinib'  [PASS]
  ✓ LEVEL_R1                            MAF='Cetuximab,Tucatinib+Trastuzumab,Panitumumab'  API='Cetuximab,Tucatinib,Trastuzumab,Panitumumab'  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_1'  API='LEVEL_1'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF='LEVEL_R1'  API='LEVEL_R1'  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF='LEVEL_Dx2'  API='LEVEL_Dx2'  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[12/18]  TUMOR  |  KRAS  |  KRAS_G12D  |  p.G12D
  → Mutation query: KRAS G12D
  ✓ ONCOKB_ALLELE_EXIST                 MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='KRAS, a GTPase which functions as an upstream regulator'  API='KRAS, a GTPase which functions as an upstream regulator'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda2'  API='LEVEL_Fda2'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='The KRAS G12D mutation is known to be oncogenic.'  API='The KRAS G12D mutation is known to be oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Gain-of-function'  API='Gain-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF='Avutometinib+Defactinib'  API='Avutometinib,Defactinib'  [PASS]
  ✓ LEVEL_2                             MAF='Cobimetinib,Trametinib'  API='Cobimetinib,Trametinib'  [PASS]
  ✓ LEVEL_3A                            MAF='Cobimetinib,Daraxonrasib,Setidegrasib,Trametinib'  API='Cobimetinib,Daraxonrasib,Setidegrasib,Trametinib'  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF='Daraxonrasib,MRTX1133,Setidegrasib,Trametinib,Cobimetin'  API='Daraxonrasib,MRTX1133,Setidegrasib,Trametinib,Cobimetin'  [PASS]
  ✓ LEVEL_R1                            MAF='Cetuximab,Tucatinib+Trastuzumab,Panitumumab'  API='Cetuximab,Tucatinib,Trastuzumab,Panitumumab'  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_1'  API='LEVEL_1'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF='LEVEL_R1'  API='LEVEL_R1'  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF='LEVEL_Dx2'  API='LEVEL_Dx2'  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[13/18]  TUMOR  |  CDK4  |  MantaDUP:CDK4_AMP  |  
  → CNA query: CDK4 AMPLIFICATION
  ✓ ONCOKB_ALLELE_EXIST                 MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='CDK4, an intracellular kinase, is altered by amplificat'  API='CDK4, an intracellular kinase, is altered by amplificat'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda3'  API='LEVEL_Fda3'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='CDK4 amplification is known to be oncogenic.'  API='CDK4 amplification is known to be oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Gain-of-function'  API='Gain-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF=''  API=''  [PASS]
  ✓ LEVEL_2                             MAF=''  API=''  [PASS]
  ✓ LEVEL_3A                            MAF=''  API=''  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF='Palbociclib,Abemaciclib'  API='Palbociclib,Abemaciclib'  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_4'  API='LEVEL_4'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF=''  API=''  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[14/18]  TUMOR  |  BRCA2  |  BRCA2_L3101R  |  p.L3101R
  → Mutation query: BRCA2 L3101R
  ✓ ONCOKB_ALLELE_EXIST                 MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='BRCA2, a tumor suppressor involved in the DNA damage re'  API='BRCA2, a tumor suppressor involved in the DNA damage re'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda2'  API='LEVEL_Fda2'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='The BRCA2 L3101R mutation is likely oncogenic.'  API='The BRCA2 L3101R mutation is likely oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Likely Loss-of-function'  API='Likely Loss-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Likely Oncogenic'  API='Likely Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF='Niraparib,Niraparib+Abiraterone Acetate+Prednisone,Olap'  API='Niraparib,Abiraterone Acetate,Prednisone,Olaparib,Abira'  [PASS]
  ✓ LEVEL_2                             MAF='Olaparib,Rucaparib,Niraparib'  API='Olaparib,Rucaparib,Niraparib'  [PASS]
  ✓ LEVEL_3A                            MAF='Olaparib,Talazoparib'  API='Olaparib,Talazoparib'  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF=''  API=''  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_1'  API='LEVEL_1'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF=''  API=''  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[15/18]  TUMOR  |  DICER1  |  DICER1_E1705K  |  p.E1705*
  → Mutation query: DICER1 E1705*
  ✓ ONCOKB_ALLELE_EXIST                 MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='DICER1, an endoribonuclease, is altered in various canc'  API='DICER1, an endoribonuclease, is altered in various canc'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF=''  API=''  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='The DICER1 E1705* is a truncating mutation in a tumor s'  API='The DICER1 E1705* is a truncating mutation in a tumor s'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='False'  API='False'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Likely Loss-of-function'  API='Likely Loss-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Likely Oncogenic'  API='Likely Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF=''  API=''  [PASS]
  ✓ LEVEL_2                             MAF=''  API=''  [PASS]
  ✓ LEVEL_3A                            MAF=''  API=''  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF=''  API=''  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF=''  API=''  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF=''  API=''  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[16/18]  TUMOR  |  TP53  |  TP53_R175H  |  p.R175H
  → Mutation query: TP53 R175H
  ✓ ONCOKB_ALLELE_EXIST                 MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='TP53, a tumor suppressor in the DNA damage pathway, is '  API='TP53, a tumor suppressor in the DNA damage pathway, is '  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF=''  API=''  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='The TP53 R175H mutation is known to be oncogenic.'  API='The TP53 R175H mutation is known to be oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Loss-of-function'  API='Loss-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF=''  API=''  [PASS]
  ✓ LEVEL_2                             MAF=''  API=''  [PASS]
  ✓ LEVEL_3A                            MAF=''  API=''  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF=''  API=''  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF=''  API=''  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF=''  API=''  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF='LEVEL_Px1'  API='LEVEL_Px1'  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[17/18]  TUMOR  |  ERBB2  |  MantaDUP:ERBB2_AMP  |  
  → CNA query: ERBB2 AMPLIFICATION
  ✓ ONCOKB_ALLELE_EXIST                 MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='ERBB2, a receptor tyrosine kinase, is altered by mutati'  API='ERBB2, a receptor tyrosine kinase, is altered by mutati'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda2'  API='LEVEL_Fda2'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='ERBB2 amplification is known to be oncogenic.'  API='ERBB2 amplification is known to be oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Gain-of-function'  API='Gain-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF='Ado-Trastuzumab Emtansine,Lapatinib,Capecitabine,Marget'  API='Ado-Trastuzumab Emtansine,Lapatinib,Capecitabine,Marget'  [PASS]
  ✓ LEVEL_2                             MAF='Ado-Trastuzumab Emtansine,Lapatinib,Trastuzumab,Trastuz'  API='Ado-Trastuzumab Emtansine,Lapatinib,Trastuzumab,Trastuz'  [PASS]
  ✓ LEVEL_3A                            MAF='Trastuzumab Deruxtecan,Zanidatamab'  API='Trastuzumab Deruxtecan,Zanidatamab'  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF=''  API=''  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_1'  API='LEVEL_1'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF=''  API=''  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[18/18]  TUMOR  |  ERBB2  |  ERBB2_S310F  |  p.S310F
  → Mutation query: ERBB2 S310F
  ✓ ONCOKB_ALLELE_EXIST                 MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_GENE_SUMMARY                 MAF='ERBB2, a receptor tyrosine kinase, is altered by mutati'  API='ERBB2, a receptor tyrosine kinase, is altered by mutati'  [PASS]
  ✓ ONCOKB_HOTSPOT                      MAF='True'  API='True'  [PASS]
  ✓ ONCOKB_VUS                          MAF='False'  API='False'  [PASS]
  ✓ ONCOKB_highestFdaLevel              MAF='LEVEL_Fda2'  API='LEVEL_Fda2'  [PASS]
  ✓ ONCOKB_VARIANT_SUMMARY              MAF='The ERBB2 S310F mutation is known to be oncogenic.'  API='The ERBB2 S310F mutation is known to be oncogenic.'  [PASS]
  ✓ GENE_IN_ONCOKB                      MAF='True'  API='True'  [PASS]
  ✓ VARIANT_IN_ONCOKB                   MAF='True'  API='True'  [PASS]
  ✓ MUTATION_EFFECT                     MAF='Gain-of-function'  API='Gain-of-function'  [PASS]
  ✓ ONCOGENIC                           MAF='Oncogenic'  API='Oncogenic'  [PASS]
  ✓ LEVEL_1                             MAF='Trastuzumab Deruxtecan'  API='Trastuzumab Deruxtecan'  [PASS]
  ✓ LEVEL_2                             MAF='Ado-Trastuzumab Emtansine,Neratinib,Neratinib+Trastuzum'  API='Ado-Trastuzumab Emtansine,Neratinib,Trastuzumab,Fulvest'  [PASS]
  ✓ LEVEL_3A                            MAF='Neratinib,Sevabertinib,Trastuzumab+Pertuzumab+Docetaxel'  API='Neratinib,Sevabertinib,Trastuzumab,Pertuzumab,Docetaxel'  [PASS]
  ✓ LEVEL_3B                            MAF=''  API=''  [PASS]
  ✓ LEVEL_4                             MAF='Neratinib,Pertuzumab+Trastuzumab,Trastuzumab Deruxtecan'  API='Neratinib,Pertuzumab,Trastuzumab,Trastuzumab Deruxtecan'  [PASS]
  ✓ LEVEL_R1                            MAF=''  API=''  [PASS]
  ✓ LEVEL_R2                            MAF=''  API=''  [PASS]
  ✓ HIGHEST_SENSITIVE_LEVEL             MAF='LEVEL_1'  API='LEVEL_1'  [PASS]
  ✓ HIGHEST_RESISTANCE_LEVEL            MAF=''  API=''  [PASS]
  ✓ ONCOKB_DIAG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_PROG_LVL                     MAF=''  API=''  [PASS]
  ✓ ONCOKB_diagnosticSummary            MAF=''  API=''  [PASS]
  ✓ ONCOKB_prognosticSummary            MAF=''  API=''  [PASS]

============================================================
OncoKB SUMMARY
============================================================
  Variants checked:  18
  Field checks:      414
  PASS:              414
  MISMATCH:          0
  ERROR:             0
============================================================

────────────────────────────────────────────────────────────
FIELD COVERAGE REPORT
────────────────────────────────────────────────────────────
  Tested (≥1 variant had data):  20
    ✓ ONCOKB_ALLELE_EXIST
    ✓ ONCOKB_GENE_SUMMARY
    ✓ ONCOKB_HOTSPOT
    ✓ ONCOKB_VUS
    ✓ ONCOKB_highestFdaLevel
    ✓ ONCOKB_VARIANT_SUMMARY
    ✓ GENE_IN_ONCOKB
    ✓ VARIANT_IN_ONCOKB
    ✓ MUTATION_EFFECT
    ✓ ONCOGENIC
    ✓ LEVEL_1
    ✓ LEVEL_2
    ✓ LEVEL_3A
    ✓ LEVEL_4
    ✓ LEVEL_R1
    ✓ LEVEL_R2
    ✓ HIGHEST_SENSITIVE_LEVEL
    ✓ HIGHEST_RESISTANCE_LEVEL
    ✓ ONCOKB_DIAG_LVL
    ✓ ONCOKB_PROG_LVL

  UNTESTABLE (always empty in MAF + API): 3
    ○ LEVEL_3B
        → No hint available.
    ○ ONCOKB_diagnosticSummary
        → Always empty when querying without a tumorType. The API returns highestDiagnosticImplicationLevel (ONCOKB_DIAG_LVL) pan-cancer, but diagnosticSummary text is only populated for a specific tumor type. Variants in this dataset that HAVE a DIAG level: IDH1 R132H (Dx2), BRAF V600E (Dx2), NRAS Q61R (Dx2), KRAS G12D/G13D (Dx2), CDKN2A_DEL (Dx2), PTEN R130* (Dx3). To test, re-query with e.g. tumorType=Glioma for IDH1 R132H.
    ○ ONCOKB_prognosticSummary
        → Always empty when querying without a tumorType. Same reason as diagnosticSummary. Variants with PROG level: IDH1 R132H (Px2), NRAS Q61R (Px1), TP53 R175H (Px1). To test, re-query with a specific tumor type.
────────────────────────────────────────────────────────────

############################################################
MODULE: VEP
############################################################

────────────────────────────────────────────────────────────
[1/18]  TUMOR  |  NRAS  |  NRAS_Q61R  |  p.Q61R
  → VEP query: ENST00000369535:c.182A>G  (refseq=False)
  ✓ Consequence          MAF='missense_variant'  API='missense_variant'  [PASS]
  ✓ IMPACT               MAF='MODERATE'  API='MODERATE'  [PASS]
  ✓ BIOTYPE              MAF='protein_coding'  API='protein_coding'  [PASS]
  ✓ EXON                 MAF='3/7'  API='3/7'  [PASS]
  ✓ INTRON               MAF=''  API=''  [PASS]
  ✓ HGVSc                MAF='c.182A>G'  API='c.182A>G'  [PASS]
  ✓ HGVSp                MAF='ENSP00000358548.4:p.Gln61Arg'  API='ENSP00000358548:p.Gln61Arg'  [PASS]
  ✓ Protein_position     MAF='61'  API='61'  [PASS]
  ✓ Amino_acids          MAF='Q/R'  API='Q/R'  [PASS]
  ✓ Codons               MAF='cAa/cGa'  API='cAa/cGa'  [PASS]
  ✓ SIFT                 MAF='deleterious_low_confidence(0.02)'  API='deleterious_low_confidence(0.02)'  [PASS]
  ✓ PolyPhen             MAF=''  API=''  [PASS]
  ✓ STRAND               MAF='-1'  API='-1'  [PASS]
  ✓ CANONICAL            MAF='YES'  API='YES'  [PASS]

────────────────────────────────────────────────────────────
[2/18]  TUMOR  |  IDH1  |  IDH1_R132H  |  p.R132H
  → VEP query: ENST00000345146:c.395G>A  (refseq=False)
  ✓ Consequence          MAF='missense_variant'  API='missense_variant'  [PASS]
  ✓ IMPACT               MAF='MODERATE'  API='MODERATE'  [PASS]
  ✓ BIOTYPE              MAF='protein_coding'  API='protein_coding'  [PASS]
  ✓ EXON                 MAF='4/10'  API='4/10'  [PASS]
  ✓ INTRON               MAF=''  API=''  [PASS]
  ✓ HGVSc                MAF='c.395G>A'  API='c.395G>A'  [PASS]
  ✓ HGVSp                MAF='ENSP00000260985.2:p.Arg132His'  API='ENSP00000260985:p.Arg132His'  [PASS]
  ✓ Protein_position     MAF='132'  API='132'  [PASS]
  ✓ Amino_acids          MAF='R/H'  API='R/H'  [PASS]
  ✓ Codons               MAF='cGt/cAt'  API='cGt/cAt'  [PASS]
  ✓ SIFT                 MAF='tolerated_low_confidence(0.07)'  API='tolerated_low_confidence(0.07)'  [PASS]
  ✓ PolyPhen             MAF='benign(0.009)'  API='benign(0.009)'  [PASS]
  ✓ STRAND               MAF='-1'  API='-1'  [PASS]
  ✓ CANONICAL            MAF='YES'  API='YES'  [PASS]

────────────────────────────────────────────────────────────
[3/18]  TUMOR  |  PIK3CA  |  PIK3CA_H1047R  |  p.H1047R
  → VEP query: ENST00000263967:c.3140A>G  (refseq=False)
  ✓ Consequence          MAF='missense_variant'  API='missense_variant'  [PASS]
  ✓ IMPACT               MAF='MODERATE'  API='MODERATE'  [PASS]
  ✓ BIOTYPE              MAF='protein_coding'  API='protein_coding'  [PASS]
  ✓ EXON                 MAF='21/21'  API='21/21'  [PASS]
  ✓ INTRON               MAF=''  API=''  [PASS]
  ✓ HGVSc                MAF='c.3140A>G'  API='c.3140A>G'  [PASS]
  ✓ HGVSp                MAF='ENSP00000263967.3:p.His1047Arg'  API='ENSP00000263967:p.His1047Arg'  [PASS]
  ✓ Protein_position     MAF='1047'  API='1047'  [PASS]
  ✓ Amino_acids          MAF='H/R'  API='H/R'  [PASS]
  ✓ Codons               MAF='cAt/cGt'  API='cAt/cGt'  [PASS]
  ✓ SIFT                 MAF='deleterious(0)'  API='deleterious(0)'  [PASS]
  ✓ PolyPhen             MAF='benign(0.088)'  API='benign(0.088)'  [PASS]
  ✓ STRAND               MAF='1'  API='1'  [PASS]
  ✓ CANONICAL            MAF='YES'  API='YES'  [PASS]

────────────────────────────────────────────────────────────
[4/18]  TUMOR  |  EGFR  |  EGFR_del19  |  p.Glu746_Ala750del
  → VEP query: ENST00000275493:c.2236_2250del  (refseq=False)
  ✓ Consequence          MAF='inframe_deletion'  API='inframe_deletion'  [PASS]
  ✓ IMPACT               MAF='MODERATE'  API='MODERATE'  [PASS]
  ✓ BIOTYPE              MAF='protein_coding'  API='protein_coding'  [PASS]
  ✓ EXON                 MAF='19/28'  API='19/28'  [PASS]
  ✓ INTRON               MAF=''  API=''  [PASS]
  ✓ HGVSc                MAF='c.2236_2250del'  API='c.2236_2250del'  [PASS]
  ✓ HGVSp                MAF='ENSP00000275493.2:p.Glu746_Ala750del'  API='ENSP00000275493:p.Glu746_Ala750del'  [PASS]
  ✓ Protein_position     MAF='746-750'  API='746-750'  [PASS]
  ✓ Amino_acids          MAF='ELREA/-'  API='ELREA/-'  [PASS]
  ✓ Codons               MAF='GAATTAAGAGAAGCA/-'  API='GAATTAAGAGAAGCA/-'  [PASS]
  ✓ SIFT                 MAF=''  API=''  [PASS]
  ✓ PolyPhen             MAF=''  API=''  [PASS]
  ✓ STRAND               MAF='1'  API='1'  [PASS]
  ✓ CANONICAL            MAF='YES'  API='YES'  [PASS]

────────────────────────────────────────────────────────────
[5/18]  TUMOR  |  EGFR  |  EGFR_L858R  |  p.L858R
  → VEP query: ENST00000275493:c.2573T>G  (refseq=False)
  ✓ Consequence          MAF='missense_variant'  API='missense_variant'  [PASS]
  ✓ IMPACT               MAF='MODERATE'  API='MODERATE'  [PASS]
  ✓ BIOTYPE              MAF='protein_coding'  API='protein_coding'  [PASS]
  ✓ EXON                 MAF='21/28'  API='21/28'  [PASS]
  ✓ INTRON               MAF=''  API=''  [PASS]
  ✓ HGVSc                MAF='c.2573T>G'  API='c.2573T>G'  [PASS]
  ✓ HGVSp                MAF='ENSP00000275493.2:p.Leu858Arg'  API='ENSP00000275493:p.Leu858Arg'  [PASS]
  ✓ Protein_position     MAF='858'  API='858'  [PASS]
  ✓ Amino_acids          MAF='L/R'  API='L/R'  [PASS]
  ✓ Codons               MAF='cTg/cGg'  API='cTg/cGg'  [PASS]
  ✓ SIFT                 MAF='deleterious_low_confidence(0)'  API='deleterious_low_confidence(0)'  [PASS]
  ✓ PolyPhen             MAF='probably_damaging(0.997)'  API='probably_damaging(0.997)'  [PASS]
  ✓ STRAND               MAF='1'  API='1'  [PASS]
  ✓ CANONICAL            MAF='YES'  API='YES'  [PASS]

────────────────────────────────────────────────────────────
[6/18]  TUMOR  |  MET  |  MantaDUP:MET_AMP  |  
  SKIP: CNA — VEP REST API does not annotate structural variants

────────────────────────────────────────────────────────────
[7/18]  TUMOR  |  BRAF  |  BRAF_V600E  |  p.V600E
  → VEP query: ENST00000646891:c.1799T>A  (refseq=False)
  ✓ Consequence          MAF='missense_variant'  API='missense_variant'  [PASS]
  ✓ IMPACT               MAF='MODERATE'  API='MODERATE'  [PASS]
  ✓ BIOTYPE              MAF='protein_coding'  API='protein_coding'  [PASS]
  ✓ EXON                 MAF='15/18'  API='15/18'  [PASS]
  ✓ INTRON               MAF=''  API=''  [PASS]
  ✓ HGVSc                MAF='c.1799T>A'  API='c.1799T>A'  [PASS]
  ✓ HGVSp                MAF='ENSP00000493543.1:p.Val600Glu'  API='ENSP00000493543:p.Val600Glu'  [PASS]
  ✓ Protein_position     MAF='600'  API='600'  [PASS]
  ✓ Amino_acids          MAF='V/E'  API='V/E'  [PASS]
  ✓ Codons               MAF='gTg/gAg'  API='gTg/gAg'  [PASS]
  ✓ SIFT                 MAF='deleterious_low_confidence(0)'  API='deleterious_low_confidence(0)'  [PASS]
  ✓ PolyPhen             MAF='probably_damaging(0.935)'  API='probably_damaging(0.935)'  [PASS]
  ✓ STRAND               MAF='-1'  API='-1'  [PASS]
  ✓ CANONICAL            MAF='YES'  API='YES'  [PASS]

────────────────────────────────────────────────────────────
[8/18]  TUMOR  |  CDKN2A  |  MantaDEL:CDKN2A_DEL  |  
  SKIP: CNA — VEP REST API does not annotate structural variants

────────────────────────────────────────────────────────────
[9/18]  TUMOR  |  GNAQ  |  GNAQ_Q209L  |  p.Gln209=
  SKIP: missing HGVSc ('') or Feature ('ENST00000286548')

────────────────────────────────────────────────────────────
[10/18]  TUMOR  |  PTEN  |  PTEN_R130stop  |  p.R130*
  → VEP query: ENST00000371953:c.388C>T  (refseq=False)
  ✓ Consequence          MAF='stop_gained'  API='stop_gained'  [PASS]
  ✓ IMPACT               MAF='HIGH'  API='HIGH'  [PASS]
  ✓ BIOTYPE              MAF='protein_coding'  API='protein_coding'  [PASS]
  ✓ EXON                 MAF='5/9'  API='5/9'  [PASS]
  ✓ INTRON               MAF=''  API=''  [PASS]
  ✓ HGVSc                MAF='c.388C>T'  API='c.388C>T'  [PASS]
  ✓ HGVSp                MAF='ENSP00000361021.3:p.Arg130Ter'  API='ENSP00000361021:p.Arg130Ter'  [PASS]
  ✓ Protein_position     MAF='130'  API='130'  [PASS]
  ✓ Amino_acids          MAF='R/*'  API='R/*'  [PASS]
  ✓ Codons               MAF='Cga/Tga'  API='Cga/Tga'  [PASS]
  ✓ SIFT                 MAF=''  API=''  [PASS]
  ✓ PolyPhen             MAF=''  API=''  [PASS]
  ✓ STRAND               MAF='1'  API='1'  [PASS]
  ✓ CANONICAL            MAF='YES'  API='YES'  [PASS]

────────────────────────────────────────────────────────────
[11/18]  TUMOR  |  KRAS  |  KRAS_G13D  |  p.G13D
  → VEP query: ENST00000311936:c.38G>A  (refseq=False)
  ✓ Consequence          MAF='missense_variant'  API='missense_variant'  [PASS]
  ✓ IMPACT               MAF='MODERATE'  API='MODERATE'  [PASS]
  ✓ BIOTYPE              MAF='protein_coding'  API='protein_coding'  [PASS]
  ✓ EXON                 MAF='2/5'  API='2/5'  [PASS]
  ✓ INTRON               MAF=''  API=''  [PASS]
  ✓ HGVSc                MAF='c.38G>A'  API='c.38G>A'  [PASS]
  ✓ HGVSp                MAF='ENSP00000308495.3:p.Gly13Asp'  API='ENSP00000308495:p.Gly13Asp'  [PASS]
  ✓ Protein_position     MAF='13'  API='13'  [PASS]
  ✓ Amino_acids          MAF='G/D'  API='G/D'  [PASS]
  ✓ Codons               MAF='gGc/gAc'  API='gGc/gAc'  [PASS]
  ✓ SIFT                 MAF='deleterious_low_confidence(0.04)'  API='deleterious_low_confidence(0.04)'  [PASS]
  ✓ PolyPhen             MAF=''  API=''  [PASS]
  ✓ STRAND               MAF='-1'  API='-1'  [PASS]
  ✓ CANONICAL            MAF='YES'  API='YES'  [PASS]

────────────────────────────────────────────────────────────
[12/18]  TUMOR  |  KRAS  |  KRAS_G12D  |  p.G12D
  → VEP query: ENST00000311936:c.35G>A  (refseq=False)
  ✓ Consequence          MAF='missense_variant'  API='missense_variant'  [PASS]
  ✓ IMPACT               MAF='MODERATE'  API='MODERATE'  [PASS]
  ✓ BIOTYPE              MAF='protein_coding'  API='protein_coding'  [PASS]
  ✓ EXON                 MAF='2/5'  API='2/5'  [PASS]
  ✓ INTRON               MAF=''  API=''  [PASS]
  ✓ HGVSc                MAF='c.35G>A'  API='c.35G>A'  [PASS]
  ✓ HGVSp                MAF='ENSP00000308495.3:p.Gly12Asp'  API='ENSP00000308495:p.Gly12Asp'  [PASS]
  ✓ Protein_position     MAF='12'  API='12'  [PASS]
  ✓ Amino_acids          MAF='G/D'  API='G/D'  [PASS]
  ✓ Codons               MAF='gGt/gAt'  API='gGt/gAt'  [PASS]
  ✓ SIFT                 MAF='deleterious_low_confidence(0.04)'  API='deleterious_low_confidence(0.04)'  [PASS]
  ✓ PolyPhen             MAF=''  API=''  [PASS]
  ✓ STRAND               MAF='-1'  API='-1'  [PASS]
  ✓ CANONICAL            MAF='YES'  API='YES'  [PASS]

────────────────────────────────────────────────────────────
[13/18]  TUMOR  |  CDK4  |  MantaDUP:CDK4_AMP  |  
  SKIP: CNA — VEP REST API does not annotate structural variants

────────────────────────────────────────────────────────────
[14/18]  TUMOR  |  BRCA2  |  BRCA2_L3101R  |  p.L3101R
  → VEP query: ENST00000380152:c.9302T>G  (refseq=False)
  ✓ Consequence          MAF='missense_variant'  API='missense_variant'  [PASS]
  ✓ IMPACT               MAF='MODERATE'  API='MODERATE'  [PASS]
  ✓ BIOTYPE              MAF='protein_coding'  API='protein_coding'  [PASS]
  ✓ EXON                 MAF='25/27'  API='25/27'  [PASS]
  ✓ INTRON               MAF=''  API=''  [PASS]
  ✓ HGVSc                MAF='c.9302T>G'  API='c.9302T>G'  [PASS]
  ✓ HGVSp                MAF='ENSP00000369497.3:p.Leu3101Arg'  API='ENSP00000369497:p.Leu3101Arg'  [PASS]
  ✓ Protein_position     MAF='3101'  API='3101'  [PASS]
  ✓ Amino_acids          MAF='L/R'  API='L/R'  [PASS]
  ✓ Codons               MAF='cTg/cGg'  API='cTg/cGg'  [PASS]
  ✓ SIFT                 MAF='deleterious(0)'  API='deleterious(0)'  [PASS]
  ✓ PolyPhen             MAF='possibly_damaging(0.888)'  API='possibly_damaging(0.888)'  [PASS]
  ✓ STRAND               MAF='1'  API='1'  [PASS]
  ✓ CANONICAL            MAF='YES'  API='YES'  [PASS]

────────────────────────────────────────────────────────────
[15/18]  TUMOR  |  DICER1  |  DICER1_E1705K  |  p.E1705*
  → VEP query: NM_030621.4:c.5113G>T  (refseq=True)
  ✓ Consequence          MAF='stop_gained'  API='stop_gained'  [PASS]
  ✓ IMPACT               MAF='HIGH'  API='HIGH'  [PASS]
  ✓ BIOTYPE              MAF='protein_coding'  API='protein_coding'  [PASS]
  ✓ EXON                 MAF='26/29'  API='26/29'  [PASS]
  ✓ INTRON               MAF=''  API=''  [PASS]
  ✓ HGVSc                MAF='c.5113G>T'  API='c.5113G>T'  [PASS]
  ✓ HGVSp                MAF='NP_085124.2:p.Glu1705Ter'  API='NP_085124:p.Glu1705Ter'  [PASS]
  ✓ Protein_position     MAF='1705'  API='1705'  [PASS]
  ✓ Amino_acids          MAF='E/*'  API='E/*'  [PASS]
  ✓ Codons               MAF='Gaa/Taa'  API='Gaa/Taa'  [PASS]
  ✓ SIFT                 MAF=''  API=''  [PASS]
  ✓ PolyPhen             MAF=''  API=''  [PASS]
  ✓ STRAND               MAF='-1'  API='-1'  [PASS]
  ✓ CANONICAL            MAF=''  API=''  [PASS]

────────────────────────────────────────────────────────────
[16/18]  TUMOR  |  TP53  |  TP53_R175H  |  p.R175H
  → VEP query: ENST00000269305:c.524G>A  (refseq=False)
  ✓ Consequence          MAF='missense_variant'  API='missense_variant'  [PASS]
  ✓ IMPACT               MAF='MODERATE'  API='MODERATE'  [PASS]
  ✓ BIOTYPE              MAF='protein_coding'  API='protein_coding'  [PASS]
  ✓ EXON                 MAF='5/11'  API='5/11'  [PASS]
  ✓ INTRON               MAF=''  API=''  [PASS]
  ✓ HGVSc                MAF='c.524G>A'  API='c.524G>A'  [PASS]
  ✓ HGVSp                MAF='ENSP00000269305.4:p.Arg175His'  API='ENSP00000269305:p.Arg175His'  [PASS]
  ✓ Protein_position     MAF='175'  API='175'  [PASS]
  ✓ Amino_acids          MAF='R/H'  API='R/H'  [PASS]
  ✓ Codons               MAF='cGc/cAc'  API='cGc/cAc'  [PASS]
  ✓ SIFT                 MAF='tolerated(0.08)'  API='tolerated(0.08)'  [PASS]
  ✓ PolyPhen             MAF=''  API=''  [PASS]
  ✓ STRAND               MAF='-1'  API='-1'  [PASS]
  ✓ CANONICAL            MAF='YES'  API='YES'  [PASS]

────────────────────────────────────────────────────────────
[17/18]  TUMOR  |  ERBB2  |  MantaDUP:ERBB2_AMP  |  
  SKIP: CNA — VEP REST API does not annotate structural variants

────────────────────────────────────────────────────────────
[18/18]  TUMOR  |  ERBB2  |  ERBB2_S310F  |  p.S310F
  → VEP query: ENST00000269571:c.929C>T  (refseq=False)
  ✓ Consequence          MAF='missense_variant'  API='missense_variant'  [PASS]
  ✓ IMPACT               MAF='MODERATE'  API='MODERATE'  [PASS]
  ✓ BIOTYPE              MAF='protein_coding'  API='protein_coding'  [PASS]
  ✓ EXON                 MAF='8/27'  API='8/27'  [PASS]
  ✓ INTRON               MAF=''  API=''  [PASS]
  ✓ HGVSc                MAF='c.929C>T'  API='c.929C>T'  [PASS]
  ✓ HGVSp                MAF='ENSP00000269571.4:p.Ser310Phe'  API='ENSP00000269571:p.Ser310Phe'  [PASS]
  ✓ Protein_position     MAF='310'  API='310'  [PASS]
  ✓ Amino_acids          MAF='S/F'  API='S/F'  [PASS]
  ✓ Codons               MAF='tCc/tTc'  API='tCc/tTc'  [PASS]
  ✓ SIFT                 MAF='deleterious(0)'  API='deleterious(0)'  [PASS]
  ✓ PolyPhen             MAF='probably_damaging(0.984)'  API='probably_damaging(0.984)'  [PASS]
  ✓ STRAND               MAF='1'  API='1'  [PASS]
  ✓ CANONICAL            MAF='YES'  API='YES'  [PASS]

============================================================
VEP SUMMARY
============================================================
  Rows in MAF:       18
  Skipped (CNA/etc): 5
  Variants checked:  13
  Field checks:      182
  PASS:              182
  MISMATCH:          0
  ERROR:             0
============================================================

────────────────────────────────────────────────────────────
FIELD COVERAGE REPORT
────────────────────────────────────────────────────────────
  Tested (≥1 variant had data):  13
    ✓ Consequence
    ✓ IMPACT
    ✓ BIOTYPE
    ✓ EXON
    ✓ HGVSc
    ✓ HGVSp
    ✓ Protein_position
    ✓ Amino_acids
    ✓ Codons
    ✓ SIFT
    ✓ PolyPhen
    ✓ STRAND
    ✓ CANONICAL

  UNTESTABLE (always empty in MAF + API): 1
    ○ INTRON
        → Only populated for intronic variants; no intronic variants in this dataset.
────────────────────────────────────────────────────────────

############################################################
COMBINED SUMMARY
############################################################
  ONCOKB    PASS=414  MISMATCH=0  ERROR=0
  VEP       PASS=182  MISMATCH=0  ERROR=0
  TOTAL     PASS=596  MISMATCH=0  ERROR=0
############################################################

Results written to: results.tsv
```