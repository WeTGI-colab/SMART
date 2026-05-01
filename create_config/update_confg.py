#!/usr/bin/env python3

# =============================================================================
# build_config.py
# =============================================================================
#
# DESCRIPTION:
#   Reads ConfigORIGINAL.yaml (the human-maintained field metadata template)
#   and produces Config.yaml (the generated config consumed by merge_maf.py).
#
#   The script performs three tasks:
#     1. Assigns a tier to every field based on the sets defined below
#     2. Updates ClinVar_* field descriptions from the ClinVar VCF header
#     3. Writes the result to Config.yaml
#
# TIER DEFINITIONS:
#   drop   : Redundant, duplicated, or internal fields. Excluded from all outputs.
#   tier 1 : All non-drop fields. Present in Final_result_tier1.maf.
#   tier 2 : ~400 fields for bioinformaticians. Includes technical depth,
#            population frequencies, full ClinVar/CIViC/OncoKB actionability,
#            SV fields, and a curated subset of ONCOKB_TX_*/ONCOKB_DIAG_* columns.
#   tier 3 : ~200 fields for clinical scientists. Core variant identity,
#            consequence, top-line actionability, ClinVar significance,
#            hotspot flags, population frequency summaries, and key
#            OncoKB/CIViC clinical fields. No expanded TX/DIAG detail.
#
# USAGE:
#   python build_config.py
#
# RELATED FILES:
#   - ConfigORIGINAL.yaml : human-maintained template
#   - Config.yaml         : generated output
#   - merge_maf.py        : consumes Config.yaml
# =============================================================================

import yaml
import gzip
import re
import fnmatch

INPUT_YAML   = "ConfigORIGINAL.yaml"
OUTPUT_YAML  = "Config.yaml"
CLINVAR_VCF  = "/mnt/data1/vep_test/ClinVar/clinvar.vcf.gz"
HOTSPOTS_VCF = "/mnt/data1/vep_test/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz"


# =============================================================================
# DROP FIELDS
# Redundant, duplicated, raw JSON, or internal query fields.
# Excluded from every output file.
# =============================================================================

DROP_FIELDS = {
    # Raw/redundant CHROM duplicate
    "CHROM",

    # OncoKB fields superseded by expanded ONCOKB_TX_*/ONCOKB_DIAG_* columns
    "ONCOKB_highestSensitiveLevel",
    "ONCOKB_highestResistanceLevel",
    "ONCOKB_highestDiagnosticImplicationLevel",
    "ONCOKB_highestPrognosticImplicationLevel",
    "ONCOKB_geneExist",
    "ONCOKB_GENE_EXIST",
    "ONCOKB_VARIANT_EXIST",
    "ONCOKB_variantExist",
    "ONCOKB_alleleExist",
    "ONCOKB_ONCOGENIC",
    "ONCOKB_oncogenic",
    "ONCOKB_EFFECT",
    "ONCOKB_mutationEffect.knownEffect",
    "ONCOKB_treatments",
    "ONCOKB_diagnosticImplications",
    "ONCOKB_prognosticImplications",
    "ONCOKB_dataVersion",
    "ONCOKB_hotspot",
    "ONCOKB_geneSummary",
    "ONCOKB_variantSummary",
    "ONCOKB_tumorTypeSummary",
    "ONCOKB_lastUpdate",
    "ONCOKB_VUS.1",

    # OncoKB query internals — not useful for interpretation
    "ONCOKB_query.referenceGenome",
    "ONCOKB_query.hugoSymbol",
    "ONCOKB_query.entrezGeneId",
    "ONCOKB_query.alteration",
    "ONCOKB_query.tumorType",
    "ONCOKB_query.id",
    "ONCOKB_query.alterationType",
    "ONCOKB_query.svType",
    "ONCOKB_query.consequence",
    "ONCOKB_query.proteinStart",
    "ONCOKB_query.proteinEnd",
    "ONCOKB_query.hgvs",
    "ONCOKB_query.hgvsInfo",
    "ONCOKB_query.canonicalTranscript",

    # Truncated legacy OncoKB column names (superseded by uppercase equivalents)
    "ONCOKB_genefx",
    "ONCOKB_variant",
    "ONCOKB_allelef",
    "ONCOKB_oncog",
    "ONCOKB_highes",
    "ONCOKB_others",
    "ONCOKB_geneS",
    "ONCOKB_variantf",
    "ONCOKB_tumor1",
    "ONCOKB_progno",
    "ONCOKB_diagno",
    "ONCOKB_treatm",
    "ONCOKB_dataVe",
    "ONCOKB_lastUp",
    "ONCOKB_exon",

    # Raw unexpanded annotation columns (content moved to _CSQ_ subfields)
    "CIViC_CSQ",
    "ONCOKB_JSON",
    "ONCOKB_TREATMENTS",
    "ONCOKB_DIAGNOSTIC_IMPLICATIONS",
    "ONCOKB_PROGNOSTIC_IMPLICATIONS",

    # Misc
    "ClinVar_GENEINFO_",
    "ONCOKB_VUS.1",
}
print("Drop:",len(DROP_FIELDS))

# =============================================================================
# TIER 2 FIELDS — Bioinformatician friendly (~400 columns)
#
# Rationale:
#   - Core variant identity and consequence fields
#   - Full transcript annotation detail
#   - Population frequencies (gnomAD exomes + genomes, per-population)
#   - Complete ClinVar classification and review fields
#   - Complete CIViC evidence fields
#   - OncoKB actionability levels and summaries
#   - SpliceAI scores and positions
#   - Structural variant fields
#   - FORMAT/genotype fields
#   - A curated subset of ONCOKB_TX_*/ONCOKB_DIAG_* columns
#     (see TIER2_ONCOKB_TX_PATTERNS and TIER2_ONCOKB_DIAG_PATTERNS below)
# =============================================================================

TIER2_FIELDS = {
    # --- VCF core ---
    "Tumor_Sample_Barcode",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "NCBI_Build",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "VAF",

    # --- VEP consequence ---
    "HGVSp_Short",
    "HGVSp",
    "HGVSc",
    "NM_Transcript",
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "DISTANCE",
    "STRAND",
    "FLAGS",
    "VARIANT_CLASS",
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "CANONICAL",
    "MANE",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "TSL",
    "APPRIS",
    "CCDS",
    "ENSP",
    "SWISSPROT",
    "TREMBL",
    "UNIPARC",
    "UNIPROT_ISOFORM",
    "SOURCE",
    "GENE_PHENO",
    "DOMAINS",
    "miRNA",
    "HGVS_OFFSET",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",

    # --- In silico predictors ---
    "SIFT",
    "PolyPhen",
    "REVEL",
    "LOEUF",

    # --- SpliceAI ---
    "SpliceAI_cutoff",
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
    "SpliceAI_pred_DP_AG",
    "SpliceAI_pred_DP_AL",
    "SpliceAI_pred_DP_DG",
    "SpliceAI_pred_DP_DL",
    "SpliceAI_pred_SYMBOL",

    # --- Population frequencies ---
    "AF",
    "AFR_AF",
    "AMR_AF",
    "EAS_AF",
    "EUR_AF",
    "SAS_AF",
    "MAX_AF",
    "MAX_AF_POPS",
    "gnomADe_AF",
    "gnomADe_AFR_AF",
    "gnomADe_AMR_AF",
    "gnomADe_ASJ_AF",
    "gnomADe_EAS_AF",
    "gnomADe_FIN_AF",
    "gnomADe_MID_AF",
    "gnomADe_NFE_AF",
    "gnomADe_REMAINING_AF",
    "gnomADe_SAS_AF",
    "gnomADg_AF",
    "gnomADg_AFR_AF",
    "gnomADg_AMI_AF",
    "gnomADg_AMR_AF",
    "gnomADg_ASJ_AF",
    "gnomADg_EAS_AF",
    "gnomADg_FIN_AF",
    "gnomADg_MID_AF",
    "gnomADg_NFE_AF",
    "gnomADg_REMAINING_AF",
    "gnomADg_SAS_AF",

    # --- ClinVar ---
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
    "ClinVar",
    "ClinVar_AF_ESP",
    "ClinVar_AF_EXAC",
    "ClinVar_AF_TGP",
    "ClinVar_ALLELEID",
    "ClinVar_CLNDN",
    "ClinVar_CLNDNINCL",
    "ClinVar_CLNDISDB",
    "ClinVar_CLNDISDBINCL",
    "ClinVar_CLNHGVS",
    "ClinVar_CLNREVSTAT",
    "ClinVar_CLNSIG",
    "ClinVar_CLNSIGCONF",
    "ClinVar_CLNSIGINCL",
    "ClinVar_CLNSIGSCV",
    "ClinVar_CLNVC",
    "ClinVar_CLNVCSO",
    "ClinVar_CLNVI",
    "ClinVar_DBVARID",
    "ClinVar_GENEINFO",
    "ClinVar_MC",
    "ClinVar_ONCDN",
    "ClinVar_ONCDNINCL",
    "ClinVar_ONCDISDB",
    "ClinVar_ONCDISDBINCL",
    "ClinVar_ONC",
    "ClinVar_ONCINCL",
    "ClinVar_ONCREVSTAT",
    "ClinVar_ONCSCV",
    "ClinVar_ONCCONF",
    "ClinVar_ORIGIN",
    "ClinVar_RS",
    "ClinVar_SCIDN",
    "ClinVar_SCIDNINCL",
    "ClinVar_SCIDISDB",
    "ClinVar_SCIDISDBINCL",
    "ClinVar_SCIREVSTAT",
    "ClinVar_SCI",
    "ClinVar_SCIINCL",
    "ClinVar_SCISCV",

    # --- CIViC ---
    "CIViC",
    "CIViC_GN",
    "CIViC_VT",
    "CIViC_Entity_ID",
    "CIViC_CSQ_Allele",
    "CIViC_CSQ_Consequence",
    "CIViC_CSQ_SYMBOL",
    "CIViC_CSQ_HGVSc",
    "CIViC_CSQ_HGVSp",
    "CIViC_CSQ_CIViC Variant Name",
    "CIViC_CSQ_CIViC Variant ID",
    "CIViC_CSQ_CIViC Variant URL",
    "CIViC_CSQ_CIViC Molecular Profile Name",
    "CIViC_CSQ_CIViC Molecular Profile Score",
    "CIViC_CSQ_CIViC Entity Type",
    "CIViC_CSQ_CIViC Entity Disease",
    "CIViC_CSQ_CIViC Entity Significance",
    "CIViC_CSQ_CIViC Entity Direction",
    "CIViC_CSQ_CIViC Entity Therapies",
    "CIViC_CSQ_CIViC Entity Therapy Interaction Type",
    "CIViC_CSQ_CIViC Evidence Level",
    "CIViC_CSQ_CIViC Evidence Rating",
    "CIViC_CSQ_CIViC Assertion AMP Category",
    "CIViC_CSQ_CIViC Assertion NCCN Guideline",
    "CIViC_CSQ_CIViC Assertion Regulatory Approval",
    "CIViC_CSQ_CIViC Assertion FDA Companion Test",
    "CIViC_CSQ_CIViC Assertion ACMG Codes",
    "CIViC_CSQ_CIViC Entity Status",
    "CIViC_CSQ_CIViC Entity Variant Origin",
    "CIViC_CSQ_CIViC Entity Source",
    "CIViC_CSQ_CIViC Entity URL",
    "CIViC_CSQ_ClinVar IDs",
    "CIViC_CSQ_Allele Registry ID",
    "CIViC_CSQ_CIViC HGVS",
    "CIViC_CSQ_CIViC Evidence Phenotypes",
    "CIViC_CSQ_CIViC Molecular Profile Aliases",
    "CIViC_CSQ_CIViC Molecular Profile URL",
    "CIViC_CSQ_CIViC Molecular Profile ID",
    "CIViC_CSQ_CIViC Variant Aliases",
    "CIViC_CSQ_Entrez Gene ID",
    "CIViC_CSQ_Feature_type",
    "CIViC_CSQ_Feature",

    # --- Cancer Hotspots ---
    "CancerHotspots",
    "CancerHotspots_HOTSPOT",
    "CancerHotspots_HOTSPOT_GENE",
    "CancerHotspots_HOTSPOT_HGVSp",
    "CancerHotspots_HOTSPOT3D",
    "CancerHotspots_HOTSPOT3D_GENE",
    "CancerHotspots_HOTSPOT3D_HGVSp",
    "CancerHotspots_HOTSPOTNC",
    "CancerHotspots_HOTSPOTNC_GENE",
    "CancerHotspots_HOTSPOTNC_HGVSc",
    "hotspot",

    # --- OncoKB actionability ---
    "ANNOTATED",
    "GENE_IN_ONCOKB",
    "VARIANT_IN_ONCOKB",
    "ONCOGENIC",
    "MUTATION_EFFECT",
    "MUTATION_EFFECT_CITATIONS",
    "ONCOKB_ONCOGENIC",
    "ONCOKB_HOTSPOT",
    "ONCOKB_SENS_LVL",
    "ONCOKB_FDA_LVL",
    "ONCOKB_DIAG_LVL",
    "ONCOKB_EFFECT_DESC",
    "ONCOKB_PMIDS",
    "ONCOKB_DATA_VERSION",
    "ONCOKB_LAST_UPDATE",
    "ONCOKB_GENE_SUMMARY",
    "ONCOKB_TUMOR_TYPE_SUMMARY",
    "ONCOKB_VARIANT_SUMMARY",
    "ONCOKB_ALLELE_EXIST",
    "ONCOKB_QUERY_ALTERATION",
    "ONCOKB_QUERY_TUMOR_TYPE",
    "ONCOKB_QUERY_TYPE",
    "ONCOKB_QUERY_REF_GENOME",
    "ONCOKB_QUERY_HUGO_SYMBOL",
    "ONCOKB_QUERY_ENTREZ_GENE_ID",
    "ONCOKB_VUS",
    "ONCOKB_otherSignificantSensitiveLevels",
    "ONCOKB_otherSignificantResistanceLevels",
    "ONCOKB_mutationEffect.description",
    "ONCOKB_mutationEffect.citations.pmids",
    "ONCOKB_mutationEffect.citations.abstracts",
    "ONCOKB_prognosticSummary",
    "ONCOKB_diagnosticSummary",
    "ONCOKB_highestFdaLevel",
    "ONCOKB_ALLELE_EXIST",
    "Hugo_Symbol",
    "NCBI_Build",
    "LEVEL_1",
    "LEVEL_2",
    "LEVEL_3A",
    "LEVEL_3B",
    "LEVEL_4",
    "LEVEL_R1",
    "LEVEL_R2",
    "HIGHEST_LEVEL",
    "HIGHEST_SENSITIVE_LEVEL",
    "HIGHEST_RESISTANCE_LEVEL",
    "TX_CITATIONS",
    "LEVEL_Dx1",
    "LEVEL_Dx2",
    "LEVEL_Dx3",
    "HIGHEST_DX_LEVEL",
    "DX_CITATIONS",
    "LEVEL_Px1",
    "LEVEL_Px2",
    "LEVEL_Px3",
    "HIGHEST_PX_LEVEL",
    "PX_CITATIONS",

    # --- FORMAT / genotype ---
    "FORMAT",
    "FORMAT_DATA",
    "GT",
    "AD",
    "VAF",
    "F1R2",
    "F2R1",
    "SB",
    "DP",
    "MB",
    "PS",
    "PR",
    "SR",
    "SQ",
    "OBC",
    "OBPa",
    "OBParc",
    "OBPsnp",

    # --- Structural variants ---
    "SVTYPE",
    "SVLEN",
    "END",
    "BND_DEPTH",
    "MATE_BND_DEPTH",
    "MATEID",
    "CIPOS",
    "CIEND",
    "CIGAR",
    "HOMLEN",
    "HOMSEQ",
    "IMPRECISE",
    "LEFT_SVINSSEQ",
    "RIGHT_SVINSSEQ",
    "DUPSVLEN",
    "FractionInformativeReads",
    "MQ",
    "OriginalContig",
    "OriginalStart",
    "ReverseComplementedAlleles",
}
print("TIER2_FIELDS:",len(TIER2_FIELDS))
# =============================================================================
# TIER 2 — ONCOKB_TX_* and ONCOKB_DIAG_* pattern subsets
#
# Rationale:
#   Only the clinically and analytically meaningful subfields are included.
#   Excluded: alterations (redundant with HGVSp), abstracts, color/code/id
#   internal taxonomy fields, and parent/level hierarchy fields.
# =============================================================================

TIER2_ONCOKB_TX_PATTERNS = [
    "ONCOKB_TX_*_level",
    "ONCOKB_TX_*_fdaLevel",
    "ONCOKB_TX_*_drugs",
    "ONCOKB_TX_*_approvedIndications",
    "ONCOKB_TX_*_pmids",
    "ONCOKB_TX_*_description",
    "ONCOKB_TX_*_levelAssociatedCancerType.name",
    "ONCOKB_TX_*_levelAssociatedCancerType.mainType.name",
    "ONCOKB_TX_*_levelAssociatedCancerType.mainType.tumorForm",
    "ONCOKB_TX_*_levelAssociatedCancerType.tissue",
    "ONCOKB_TX_*_levelExcludedCancerTypes",
]
print("TIER2_ONCOKB_TX_PATTERNS:",len(TIER2_ONCOKB_TX_PATTERNS))

TIER2_ONCOKB_DIAG_PATTERNS = [
    "ONCOKB_DIAG_*_levelOfEvidence",
    "ONCOKB_DIAG_*_description",
    "ONCOKB_DIAG_*_pmids",
    "ONCOKB_DIAG_*_tumorType.name",
    "ONCOKB_DIAG_*_tumorType.mainType.name",
    "ONCOKB_DIAG_*_tumorType.mainType.tumorForm",
    "ONCOKB_DIAG_*_tumorType.tissue",
]
print("TIER2_ONCOKB_DIAG_PATTERNS:",len(TIER2_ONCOKB_DIAG_PATTERNS))

# =============================================================================
# TIER 3 FIELDS — Clinician friendly (~200 columns)
#
# Rationale:
#   Core variant identity, top-line consequence, key population frequencies,
#   ClinVar significance and review status, hotspot flags, top-line OncoKB
#   oncogenicity and mutation effect, SpliceAI delta scores, SIFT/PolyPhen/REVEL,
#   and the most actionable CIViC fields.
#   No expanded TX/DIAG detail, no per-population gnomAD breakdown,
#   no internal query fields, no SV technical detail.
# =============================================================================

TIER3_FIELDS = {
    # --- Variant identity ---
    "Tumor_Sample_Barcode",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "REF",
    "ALT",
    "VAF",
    "FILTER",
    "VARIANT_CLASS",
    "Existing_variation",

    # --- Gene and transcript ---
    "SYMBOL",
    "HGVSp_Short",
    "HGVSp",
    "HGVSc",
    "NM_Transcript",
    "Consequence",
    "IMPACT",
    "EXON",
    "INTRON",

    # --- In silico predictors ---
    "SIFT",
    "PolyPhen",
    "REVEL",
    "LOEUF",

    # --- SpliceAI (scores only, not positions) ---
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",

    # --- Population frequency (summary only) ---
    "AF",
    "gnomADe_AF",
    "gnomADg_AF",
    "MAX_AF",
    "MAX_AF_POPS",

    # --- ClinVar (classification and review) ---
    "ClinVar_CLNSIG",
    "ClinVar_CLNREVSTAT",
    "ClinVar_CLNDN",
    "ClinVar_CLNHGVS",
    "ClinVar_ONC",
    "ClinVar_ONCREVSTAT",
    "ClinVar_SCI",
    "ClinVar_SCIREVSTAT",
    "ClinVar_ORIGIN",

    # --- Cancer Hotspots ---
    "CancerHotspots",
    "CancerHotspots_HOTSPOT",
    "hotspot",

    # --- OncoKB top-line actionability ---
    "ANNOTATED",
    "GENE_IN_ONCOKB",
    "VARIANT_IN_ONCOKB",
    "ONCOGENIC",
    "MUTATION_EFFECT",
    "ONCOKB_HOTSPOT",
    "ONCOKB_SENS_LVL",
    "ONCOKB_FDA_LVL",
    "ONCOKB_DIAG_LVL",
    "ONCOKB_DATA_VERSION",
    "ONCOKB_LAST_UPDATE",
    "HIGHEST_LEVEL",
    "HIGHEST_SENSITIVE_LEVEL",
    "HIGHEST_RESISTANCE_LEVEL",
    "HIGHEST_DX_LEVEL",
    "HIGHEST_PX_LEVEL",
    "LEVEL_1",
    "LEVEL_2",
    "LEVEL_3A",
    "LEVEL_R1",

    # --- CIViC (top-line clinical fields only) ---
    "CIViC_GN",
    "CIViC_VT",
    "CIViC_CSQ_CIViC Variant Name",
    "CIViC_CSQ_CIViC Entity Significance",
    "CIViC_CSQ_CIViC Entity Direction",
    "CIViC_CSQ_CIViC Entity Therapies",
    "CIViC_CSQ_CIViC Evidence Level",
    "CIViC_CSQ_CIViC Assertion AMP Category",
    "CIViC_CSQ_CIViC Assertion Regulatory Approval",
    "CIViC_CSQ_CIViC Entity Disease",

    # --- Genotype essentials ---
    "GT",
    "AD",

    # --- Structural variant summary ---
    "SVTYPE",
    "SVLEN",
}
print("TIER3_FIELDS:",len(TIER3_FIELDS))

# =============================================================================
# Extract ClinVar INFO descriptions from VCF header
# =============================================================================

def load_clinvar_descriptions(vcf_path):

    descriptions = {}

    with gzip.open(vcf_path, "rt") as f:
        for line in f:

            if not line.startswith("##INFO="):
                continue

            id_match   = re.search(r'ID=([^,]+)', line)
            desc_match = re.search(r'Description="([^"]+)"', line)

            if id_match and desc_match:
                descriptions[id_match.group(1)] = desc_match.group(1)

    print(f"Loaded {len(descriptions)} ClinVar INFO fields")

    return descriptions


# =============================================================================
# Update YAML descriptions using ClinVar VCF header
# =============================================================================

def update_clinvar_descriptions(fields, clinvar_info):

    updated = 0
    matched = 0

    for field_name, meta in fields.items():

        if not field_name.startswith("ClinVar_"):
            continue

        clinvar_key = field_name.replace("ClinVar_", "")

        if clinvar_key in clinvar_info:
            matched += 1
            new_description = clinvar_info[clinvar_key]
            if meta.get("description") != new_description:
                meta["description"] = new_description
                updated += 1

    print(f"ClinVar fields matched: {matched}")
    print(f"ClinVar descriptions updated: {updated}")


# =============================================================================
# Check whether a field name matches any pattern in a list
# =============================================================================

def matches_any_pattern(field, patterns):

    for pattern in patterns:
        if fnmatch.fnmatch(field, pattern):
            return True

    return False


# =============================================================================
# Main: assign tiers and write Config.yaml
# =============================================================================

def apply_tiers(input_yaml, output_yaml):

    with open(input_yaml) as f:
        data = yaml.safe_load(f)

    fields = data.get("fields", {})

    counts = {"tier1": 0, "tier2": 0, "tier3": 0, "drop": 0}

    for k, v in fields.items():
        if v is None:
            print(f"WARNING: field '{k}' has no metadata")

    for field, meta in fields.items():

        field_clean = field.strip()

        # --- Drop ---
        if field_clean in DROP_FIELDS:
            meta["tier"] = "drop"
            counts["drop"] += 1
            continue

        # --- Tier 3 (clinician) ---
        if field_clean in TIER3_FIELDS:
            meta["tier"] = 3
            counts["tier3"] += 1
            continue

        # --- Tier 2 (bioinformatician) ---
        if field_clean in TIER2_FIELDS:
            meta["tier"] = 2
            counts["tier2"] += 1
            continue

        # --- ONCOKB_TX_* and ONCOKB_DIAG_* pattern-based tier assignment ---
        if matches_any_pattern(field_clean, TIER2_ONCOKB_TX_PATTERNS + TIER2_ONCOKB_DIAG_PATTERNS):
            meta["tier"] = 2
            counts["tier2"] += 1
            continue

        # --- Default: Tier 1 ---
        meta["tier"] = 1
        counts["tier1"] += 1

    # --- Update ClinVar descriptions ---
    clinvar_info = load_clinvar_descriptions(CLINVAR_VCF)
    update_clinvar_descriptions(fields, clinvar_info)

    # =============================================================================
    # Assign tier to field_patterns
    # =============================================================================

    patterns = data.get("field_patterns", [])

    for pattern_def in patterns:

        pattern = pattern_def.get("pattern", "").strip()

        # --- Drop (rare, but keep for completeness) ---
        if pattern in DROP_FIELDS:
            pattern_def["tier"] = "drop"
            continue

        # --- Tier 3 (if you ever define pattern-based Tier3) ---
        if pattern in TIER3_FIELDS:
            pattern_def["tier"] = 3
            continue

        # --- Tier 2 via explicit pattern groups ---
        if matches_any_pattern(pattern, TIER2_ONCOKB_TX_PATTERNS + TIER2_ONCOKB_DIAG_PATTERNS):
            pattern_def["tier"] = 2
            continue

        # --- Tier 2 fallback for ONCOKB ---
        # if pattern.startswith("ONCOKB_TX_") or pattern.startswith("ONCOKB_DIAG_"):
        #     pattern_def["tier"] = 2
        #     continue

        # --- Default Tier 1 ---
        pattern_def["tier"] = 1

    
    # =============================================================================
    # ✅ VALIDATION (PUT IT HERE)
    # =============================================================================
    for p in patterns:
        if "tier" not in p:
            print(f"WARNING: pattern without tier → {p.get('pattern')}")


    # --- Write output ---
    with open(output_yaml, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False)

    # --- Write output ---
    with open(output_yaml, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False)

    print("\nTier assignment complete\n")
    print("Summary:")
    #print(f"  Tier 1 : {counts['tier1']}")
    print(f"  Tier 2 : {counts['tier2']}")
    print(f"  Tier 3 : {counts['tier3']}")
    print(f"  Dropped: {counts['drop']}")

    yaml_fields = set(fields.keys())

    print("\n--- DEBUG ---")

    print(f"TIER2 missing from YAML: {len(TIER2_FIELDS - yaml_fields)}")
    print(f"TIER3 missing from YAML: {len(TIER3_FIELDS - yaml_fields)}")
    print(f"DROP missing from YAML : {len(DROP_FIELDS - yaml_fields)}")

    print("\nExample missing Tier2 fields:")
    for f in list(TIER2_FIELDS - yaml_fields)[:10]:
        print("  ", f)

if __name__ == "__main__":
    apply_tiers(INPUT_YAML, OUTPUT_YAML)