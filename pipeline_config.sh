# --------------------------------------------------
# Clinical Genomics Pipeline configuration
# Generated automatically by install_references.sh
# --------------------------------------------------
#
# Installed by : manolodominguez
# Install date : 2026-03-29
# Install host : Manolos-MacBook-Pro.local
# --------------------------------------------------

INSTALL_DATE="2026-03-29"
REF_ROOT="/Volumes/ExternalSSD/refs"

REF_GENOME="${REF_ROOT}/reference/hg38.fa"
REF_GENOME_VERSION="UCSC_hg38"

CHAIN_FILE="${REF_ROOT}/liftover/hg19ToHg38.over.chain"
CHAIN_FILE_VERSION="UCSC_hg19ToHg38"

CLINVAR_VCF="${REF_ROOT}/ClinVar/clinvar.vcf.gz"
CLINVAR_VERSION="2026-03-21"

CIVIC_FILE="${REF_ROOT}/CIVIC/nightly-VariantSummaries.tsv"
CIVIC_VCF="${REF_ROOT}/CIVIC/civic_grch38.vcf.gz"
CIVIC_VCF_INDEX="${REF_ROOT}/CIVIC/civic_grch38.vcf.gz.tbi"
CIVIC_VERSION="civic_nightly_2026-03-29"

REVEL_FILE="${REF_ROOT}/REVEL/revel_all_chromosomes.zip"
REVEL_VERSION="v1.3"

SPLICEAI_SNV="${REF_ROOT}/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="${REF_ROOT}/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"
SPLICEAI_VERSION="SpliceAIv1.3"

LOEUF_FILE="${REF_ROOT}/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz"
LOEUF_FILE_INDEX="${REF_ROOT}/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz.tbi"
LOEUF_VERSION="custom_grch38"

VEP_CACHE_DIR="${REF_ROOT}/homo_sapiens"
VEP_CACHE_VERSION="GRCh38"

VEP_PLUGINS_DIR="${REF_ROOT}/Plugins"
VEP_PLUGINS_VERSION="git_clone_2026-03-29"

CANCER_HOTSPOTS_VCF="${REF_ROOT}/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz"
CANCER_HOTSPOTS_VERSION="charlottekyng/cancer_hotspots_master"
