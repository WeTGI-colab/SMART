#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Clinical Genomics Pipeline – Reference Data Installer
# ==============================================================================
#
# This script downloads and prepares the reference datasets required to run
# the variant annotation pipeline.
#
# Features:
#   - Idempotent: existing files are skipped
#   - Creates indexes where needed
#   - Creates a pipeline_config.sh file in the current directory
#   - Records paths and resource versions where possible
#   - Verifies checksums where available
#   - Checks disk space before starting
#
# Notes:
#   - SpliceAI requires an Illumina BaseSpace account; files must be manually
#     staged to $SPLICEAI_STAGING_DIR before running this script
#   - This script assumes wget, git, gunzip, zgrep, realpath, apptainer,
#     and tabix are available on the host
#
# ==============================================================================

# ==============================================================================
# CONFIGURATION — edit these variables to match your environment
# ==============================================================================

# Directory where this script lives (used to locate companion scripts)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Base install directory — set via first argument or fallback to current directory
if [[ $# -lt 1 ]]; then
  echo "[ERROR] Usage: $0 <base_dir>"
  echo "        Example: $0 /Volumes/ExternalSSD"
  exit 1
fi
BASE_DIR="$1"

if [[ ! -d "$BASE_DIR" ]]; then
  echo "[ERROR] Directory does not exist: $BASE_DIR"
  exit 1
fi

# Where reference files will be stored
REF_ROOT="${BASE_DIR}/refs"

# SpliceAI files must be manually downloaded from BaseSpace and placed here
# before running this script (BaseSpace requires an Illumina login)
SPLICEAI_STAGING_DIR="${BASE_DIR}/spliceai_staging"

# Minimum free disk space required to proceed (in GB)
MIN_DISK_GB=200

# ==============================================================================

echo "================================================================================"
echo " Clinical Genomics Pipeline – Reference Data Installer"
echo ""
echo " Running as user : $(whoami)"
echo " Home directory  : ${HOME}"
echo " Install root    : ${REF_ROOT}"
echo " SpliceAI source : ${SPLICEAI_STAGING_DIR}"
echo ""
echo "Resources installed include:"
echo "  • Human reference genome (GRCh38)"
echo "  • Genome indexes required by GATK and samtools"
echo "  • hg19 → GRCh38 liftover chain file"
echo "  • Public annotation databases (ClinVar, CIViC, gnomAD constraints, etc.)"
echo "  • VEP cache and plugins for offline annotation"
echo ""
echo "The installer is idempotent:"
echo "  • Existing files are detected and skipped"
echo "  • Only missing resources are downloaded"
echo "  • The script can safely be re-run if interrupted"
echo "================================================================================"

# -------------------------
# Preflight checks
# -------------------------

echo ""
echo "[PREFLIGHT] Checking available disk space on ${BASE_DIR}..."
AVAIL_GB=$(df -k "${BASE_DIR}" | awk 'NR==2 {printf "%d", $4/1048576}')
echo "[PREFLIGHT] Available: ${AVAIL_GB} GB  |  Required: ${MIN_DISK_GB} GB"

if (( AVAIL_GB < MIN_DISK_GB )); then
  echo "[ERROR] Insufficient disk space."
  echo "        Available : ${AVAIL_GB} GB"
  echo "        Required  : ${MIN_DISK_GB} GB"
  echo "        Free up space under ${BASE_DIR} and re-run."
  exit 1
fi

echo "[PREFLIGHT] Disk space OK."

# -------------------------
# Container wrappers
# -------------------------
echo ""
echo "[INIT] Setting container wrappers"

SAMTOOLS_IMG="docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
GATK_IMG="docker://broadinstitute/gatk:4.6.0.0"
VEP_IMG="docker://ensemblorg/ensembl-vep:release_114.0"
HTSLIB_IMG="docker://quay.io/biocontainers/htslib:1.21--h566b1c6_1"

SAMTOOLS="apptainer exec --unsquash --bind ${BASE_DIR}:${BASE_DIR} --bind $(pwd):$(pwd) ${SAMTOOLS_IMG} samtools"
TABIX="apptainer exec --unsquash --bind ${BASE_DIR}:${BASE_DIR} --bind $(pwd):$(pwd) ${SAMTOOLS_IMG} tabix"
BGZIP="apptainer exec --unsquash --bind ${BASE_DIR}:${BASE_DIR} --bind $(pwd):$(pwd) ${HTSLIB_IMG} bgzip"
GATK="apptainer exec --unsquash --bind ${BASE_DIR}:${BASE_DIR} --bind $(pwd):$(pwd) ${GATK_IMG} gatk"
VEP="apptainer exec --unsquash --bind ${BASE_DIR}:${BASE_DIR} --bind $(pwd):$(pwd) ${VEP_IMG}"


# -------------------------
# Helper functions
# -------------------------

log_step () {
  echo ""
  echo "================================================================================"
  echo "$1"
  echo "================================================================================"
}

download_if_missing () {
  # Usage: download_if_missing <URL> <OUTPUT_FILE> [EXPECTED_SHA256]
  local url="$1"
  local out="$2"
  local expected_sha="${3:-}"

  if [[ -f "$out" ]]; then
    echo "[SKIP] Already exists: $out"
    # Still verify checksum if one was provided
    if [[ -n "$expected_sha" ]]; then
      local actual_sha
      actual_sha=$(shasum -a 256 "$out" | awk '{print $1}')
      if [[ "$actual_sha" != "$expected_sha" ]]; then
        echo "[ERROR] Checksum mismatch for existing file: $out"
        echo "        Expected : $expected_sha"
        echo "        Got      : $actual_sha"
        exit 1
      fi
      echo "[OK] Checksum verified: $out"
    fi
    return 0
  fi

  echo "[DL] $url"
  echo "     -> $out"

  local tmp_out="${out}.tmp"

  if ! curl -L --progress-bar -o "$tmp_out" "$url"; then
    echo "[ERROR] Download failed: $url"
    rm -f "$tmp_out"
    exit 1
  fi

  # Verify checksum before moving into place
  if [[ -n "$expected_sha" ]]; then
    local actual_sha
    actual_sha=$(shasum -a 256 "$tmp_out" | awk '{print $1}')
    if [[ "$actual_sha" != "$expected_sha" ]]; then
      echo "[ERROR] Checksum mismatch after download: $url"
      echo "        Expected : $expected_sha"
      echo "        Got      : $actual_sha"
      rm -f "$tmp_out"
      exit 1
    fi
    echo "[OK] Checksum verified."
  fi

  mv "$tmp_out" "$out"
}

mkdir_if_missing () {
  local d="$1"
  if [[ -d "$d" ]]; then
    echo "[OK] Directory exists: $d"
  else
    echo "[MKDIR] Creating: $d"
    mkdir -p "$d"
  fi
}

# -------------------------
# STEP 1: Liftover chain
# -------------------------
log_step "[STEP 1] Liftover chain (hg19 → hg38)"

mkdir_if_missing "${REF_ROOT}/liftover"

CHAIN_GZ="${REF_ROOT}/liftover/hg19ToHg38.over.chain.gz"
CHAIN="${REF_ROOT}/liftover/hg19ToHg38.over.chain"

if [[ -f "$CHAIN" ]]; then
  echo "[SKIP] Liftover chain already exists: $CHAIN"
else
  download_if_missing \
    "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz" \
    "$CHAIN_GZ"

  echo "[RUN] Decompressing chain file -> $CHAIN"
  gunzip -c "$CHAIN_GZ" > "$CHAIN"
  echo "[OK] Chain file ready: $CHAIN"
fi

# -------------------------
# STEP 2: GRCh38 reference + indexes
# -------------------------
log_step "[STEP 2] GRCh38 reference genome (hg38.fa + .fai + .dict)"

mkdir_if_missing "${REF_ROOT}/reference"

REF_FA_GZ="${REF_ROOT}/liftover/hg38.fa.gz"
REF_FA="${REF_ROOT}/liftover/hg38.fa"
REF_FAI="${REF_ROOT}/liftover/hg38.fa.fai"
REF_DICT="${REF_ROOT}/liftover/hg38.dict"

if [[ -f "$REF_FA" ]]; then
  echo "[SKIP] Reference FASTA already exists: $REF_FA"
else
  download_if_missing \
    "https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz" \
    "$REF_FA_GZ"

  echo "[RUN] Decompressing reference FASTA -> $REF_FA"
  # Use -c (stdout) to preserve the .gz file, protecting against interruption
  gunzip -c "$REF_FA_GZ" > "$REF_FA"
  echo "[OK] Reference FASTA ready: $REF_FA"
fi

if [[ -f "$REF_FAI" ]]; then
  echo "[SKIP] FASTA index already exists: $REF_FAI"
else
  echo "[RUN] Creating FASTA index (.fai) with samtools"
  $SAMTOOLS faidx "$REF_FA"
  echo "[OK] FASTA index created: $REF_FAI"
fi

if [[ -f "$REF_DICT" ]]; then
  echo "[SKIP] Sequence dictionary already exists: $REF_DICT"
else
  echo "[RUN] Creating sequence dictionary (.dict) with GATK"
  $GATK CreateSequenceDictionary \
      -R "$REF_FA" \
      -O "$REF_DICT"
  echo "[OK] Sequence dictionary created: $REF_DICT"
fi

# -------------------------
# STEP 3: ClinVar
# -------------------------
log_step "[STEP 3] ClinVar (GRCh38 VCF + tabix index)"

mkdir_if_missing "${REF_ROOT}/ClinVar"
CLINVAR_VCF="${REF_ROOT}/ClinVar/clinvar.vcf.gz"
CLINVAR_TBI="${REF_ROOT}/ClinVar/clinvar.vcf.gz.tbi"

download_if_missing \
  "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" \
  "$CLINVAR_VCF"

if [[ -f "$CLINVAR_TBI" ]]; then
  echo "[SKIP] ClinVar tabix index already exists: $CLINVAR_TBI"
else
  echo "[RUN] Indexing ClinVar with tabix"
  $TABIX -p vcf "$CLINVAR_VCF"
  echo "[OK] ClinVar indexed: $CLINVAR_TBI"
fi

# -------------------------
# STEP 4: CIViC
# -------------------------
log_step "[STEP 4] CIViC (nightly download)"

mkdir_if_missing "${REF_ROOT}/CIVIC"
CIVIC_FILE="${REF_ROOT}/CIVIC/nightly-VariantSummaries.tsv"

download_if_missing \
  "https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv" \
  "$CIVIC_FILE"

# Capture version from download date since CIViC nightly has no embedded version
CIVIC_DOWNLOAD_DATE="$(date +%Y-%m-%d)"

# Convert the TSV to a bgzip-compressed, tabix-indexed VCF for VEP --custom
CIVIC_VCF="${REF_ROOT}/CIVIC/civic_grch38.vcf.gz"

if [[ -f "$CIVIC_VCF" ]]; then
  echo "[SKIP] CIViC VCF already exists: $CIVIC_VCF"
else
  echo "[RUN] Converting CIViC TSV → VCF using civic_formating.py"

  # Create temporary Apptainer-backed wrappers for bgzip and tabix
  # (the Python script calls them via subprocess, so they must be executable files)
  _WRAPPER_DIR=$(mktemp -d)
  cat > "${_WRAPPER_DIR}/bgzip" <<BGZIP_WRAPPER
#!/bin/bash
apptainer exec --unsquash --bind "${BASE_DIR}:${BASE_DIR}" ${SAMTOOLS_IMG} bgzip "\$@"
BGZIP_WRAPPER
  cat > "${_WRAPPER_DIR}/tabix" <<TABIX_WRAPPER
#!/bin/bash
apptainer exec --unsquash --bind "${BASE_DIR}:${BASE_DIR}" ${SAMTOOLS_IMG} tabix "\$@"
TABIX_WRAPPER
  chmod +x "${_WRAPPER_DIR}/bgzip" "${_WRAPPER_DIR}/tabix"

  python3 "${SCRIPT_DIR}/civic_formating.py" \
    --input      "$CIVIC_FILE" \
    --output     "$CIVIC_VCF" \
    --assembly   grch38 \
    --chain-file "$CHAIN" \
    --bgzip      "${_WRAPPER_DIR}/bgzip" \
    --tabix      "${_WRAPPER_DIR}/tabix"

  rm -rf "${_WRAPPER_DIR}"
  echo "[OK] CIViC VCF ready: $CIVIC_VCF"
  echo "[OK] CIViC index   : ${CIVIC_VCF}.tbi"
fi

# -------------------------
# STEP 5: REVEL
# -------------------------
log_step "[STEP 5] REVEL (zip download)"

mkdir_if_missing "${REF_ROOT}/REVEL"
REVEL_ZIP="${REF_ROOT}/REVEL/revel_all_chromosomes.zip"

download_if_missing \
  "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip" \
  "$REVEL_ZIP"

# Process REVEL zip → bgzip-compressed TSV + tabix index
REVEL_TSV="${REF_ROOT}/REVEL/new_tabbed_revel_grch38.tsv.gz"

if [[ -f "$REVEL_TSV" && -s "$REVEL_TSV" ]]; then
  echo "[SKIP] REVEL TSV already processed: $REVEL_TSV"
else
  # Remove empty file if it exists (failed previous run)
  [[ -f "$REVEL_TSV" ]] && rm -f "$REVEL_TSV" && echo "[INFO] Removed empty REVEL TSV, reprocessing..."
  echo "[RUN] Extracting REVEL zip..."
  unzip -o "$REVEL_ZIP" -d "${REF_ROOT}/REVEL/"

  # The zip may contain a .csv or .csv.gz — find it
  REVEL_CSV=$(find "${REF_ROOT}/REVEL/" \( -name "*.csv" -o -name "*.csv.gz" -o -name "revel_with_transcript_ids" \) | head -1)
  if [[ -z "$REVEL_CSV" ]]; then
    echo "[ERROR] No CSV found in REVEL zip contents under ${REF_ROOT}/REVEL/"
    echo "        Zip contents:"
    unzip -l "$REVEL_ZIP" | tail -n +4 | head -20
    exit 1
  fi
  echo "[RUN] Converting REVEL CSV → tab-separated + bgzip: $REVEL_CSV"

  # Add chr prefix, convert commas to tabs, sort by chrom+pos, bgzip
  # Handle both plain .csv and .csv.gz inputs
  if [[ "$REVEL_CSV" == *.gz ]]; then
    gzip -dc "$REVEL_CSV"
  else
    cat "$REVEL_CSV"
  fi \
    | awk 'NR==1{
        n=split($0,a,",");
        printf "%s", a[1];
        for(i=2;i<=n;i++) printf "\t%s", a[i];
        printf "\n"; next
      } {
        n=split($0,a,",");
        chrom = (a[1] ~ /^chr/) ? a[1] : "chr"a[1];
        printf "%s", chrom;
        for(i=2;i<=n;i++) printf "\t%s", a[i];
        printf "\n"
      }' \
    | (read header; echo "$header"; sort -k1,1 -k2,2n) \
    | $BGZIP > "$REVEL_TSV"

  echo "[RUN] Indexing REVEL with tabix..."
  $TABIX -s 1 -b 2 -e 2 --skip-lines 1 "$REVEL_TSV"
  echo "[OK] REVEL ready: $REVEL_TSV"

  # Clean up extracted CSV to save space
  rm -f "$REVEL_CSV"
  echo "[OK] Removed intermediate CSV."
fi

# Safety check: ensure BGZF compression and index
if [[ -f "$REVEL_TSV" ]]; then
  # Recompress with bgzip if not already BGZF
  if ! $TABIX --test-bgzf "$REVEL_TSV" 2>/dev/null; then
    echo "[RUN] REVEL file is not BGZF — recompressing..."
    gunzip -c "$REVEL_TSV" | $BGZIP > "${REVEL_TSV}.tmp"
    mv "${REVEL_TSV}.tmp" "$REVEL_TSV"
    echo "[OK] REVEL recompressed with bgzip."
  fi
  # Create index if missing
  if [[ ! -f "${REVEL_TSV}.tbi" ]]; then
    echo "[RUN] REVEL .tbi missing — indexing..."
    $TABIX -s 1 -b 2 -e 2 --skip-lines 1 "$REVEL_TSV"
    echo "[OK] REVEL index created."
  fi
fi

# -------------------------
# STEP 6: SpliceAI
# -------------------------
log_step "[STEP 6] SpliceAI (SNV + INDEL VCFs)"

# SpliceAI is gated behind an Illumina BaseSpace login and cannot be
# downloaded automatically. Files must be manually downloaded and placed
# in SPLICEAI_STAGING_DIR before running this script. See comments at the end of this script
#
# Expected files:
#   spliceai_scores.raw.snv.hg38.vcf.gz
#   spliceai_scores.raw.indel.hg38.vcf.gz

mkdir_if_missing "${REF_ROOT}/SpliceAI"

SPLICEAI_SNV="${REF_ROOT}/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="${REF_ROOT}/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"

SPLICEAI_SNV_SRC="${SPLICEAI_STAGING_DIR}/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL_SRC="${SPLICEAI_STAGING_DIR}/spliceai_scores.raw.indel.hg38.vcf.gz"

SPLICEAI_MISSING=0

for src_file in "$SPLICEAI_SNV_SRC" "$SPLICEAI_INDEL_SRC"; do
  if [[ ! -f "$src_file" ]]; then
    echo "[WARN] SpliceAI staging file not found: $src_file"
    SPLICEAI_MISSING=1
  fi
done

if (( SPLICEAI_MISSING == 1 )); then
  echo ""
  echo "[WARN] =================================================================="
  echo "[WARN] SpliceAI files are missing from the staging directory:"
  echo "[WARN]   ${SPLICEAI_STAGING_DIR}"
  echo "[WARN]"
  echo "[WARN] SpliceAI requires a free Illumina BaseSpace account."
  echo "[WARN] Download both files manually from:"
  echo "[WARN]   SNV   : https://basespace.illumina.com/s/otSPW8hnhaZR"
  echo "[WARN]   INDEL : https://basespace.illumina.com/s/jI4aHqRKoGqE"
  echo "[WARN]"
  echo "[WARN] Then place them in: ${SPLICEAI_STAGING_DIR}"
  echo "[WARN] and re-run this script."
  echo "[WARN] =================================================================="
else
  declare -A SPLICEAI_FILES=(
    ["$SPLICEAI_SNV_SRC"]="$SPLICEAI_SNV"
    ["$SPLICEAI_INDEL_SRC"]="$SPLICEAI_INDEL"
  )

  for src in "${!SPLICEAI_FILES[@]}"; do
    dest="${SPLICEAI_FILES[$src]}"
    if [[ -f "$dest" ]]; then
      echo "[SKIP] Already exists: $dest"
    else
      echo "[CP] $src -> $dest"
      cp "$src" "$dest"
      echo "[OK] Copied: $dest"
    fi
    # Index if .tbi is missing
    if [[ ! -f "${dest}.tbi" ]]; then
      echo "[RUN] Indexing $(basename "$dest") with tabix..."
      $TABIX -p vcf "$dest"
      echo "[OK] Index created: ${dest}.tbi"
    else
      echo "[SKIP] Index already exists: ${dest}.tbi"
    fi
  done
fi

# -------------------------
# STEP 7: gnomAD constraints (LOEUF)
# -------------------------
log_step "[STEP 7] gnomAD constraint metrics (LOEUF)"

mkdir_if_missing "${REF_ROOT}/gnomAD_constraints"

LOEUF_GRCH38="${REF_ROOT}/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz"
LOEUF_GRCH38_TBI="${REF_ROOT}/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz.tbi"

echo "[INFO] Downloading gnomAD LOEUF from Google Drive."
echo "       NOTE: Google Drive may silently return an HTML error page instead"
echo "             of the file if the link has expired or the quota is exceeded."
echo "             Verify the downloaded file is not an HTML page if indexing fails."

download_if_missing \
  "https://drive.google.com/uc?export=download&id=1pJMBwZQQLB2K4DelrWCEILzSadmLVybq" \
  "$LOEUF_GRCH38"

download_if_missing \
  "https://drive.google.com/uc?export=download&id=1ZFg9d2VW_PC3VcsaSS84LaKeudw41PT8" \
  "$LOEUF_GRCH38_TBI"

# Sanity-check: the file should be binary (gzip), not HTML
if file "$LOEUF_GRCH38" | grep -q "HTML"; then
  echo "[ERROR] LOEUF download appears to be an HTML page, not a gzip file."
  echo "        The Google Drive link may have expired or hit a download quota."
  echo "        Remove ${LOEUF_GRCH38} and obtain the file manually, then re-run."
  rm -f "$LOEUF_GRCH38"
  exit 1
fi

# -------------------------
# STEP 8: VEP cache
# -------------------------
log_step "[STEP 8] VEP cache (homo_sapiens merged, GRCh38)"

# The merged cache includes SIFT and PolyPhen precomputed scores in addition
# to transcript and regulation data. It is larger (~25 GB) but required for
# full annotation. Use this on servers where disk space is available.
VEP_CACHE_DIR="${REF_ROOT}/homo_sapiens"

VEP_CACHE_TAR="${REF_ROOT}/homo_sapiens_merged_vep_114_GRCh38.tar.gz"
VEP_CACHE_URL="https://ftp.ensembl.org/pub/release-114/variation/indexed_vep_cache/homo_sapiens_merged_vep_114_GRCh38.tar.gz"

if [[ -d "$VEP_CACHE_DIR" ]] && [[ -n "$(ls -A "$VEP_CACHE_DIR" 2>/dev/null)" ]]; then
  echo "[SKIP] VEP cache directory already present and non-empty: $VEP_CACHE_DIR"
else
  # Download with curl (supports resume with -C -, avoids FTP timeout issues)
  if [[ ! -f "$VEP_CACHE_TAR" ]]; then
    echo "[DL] Downloading VEP cache (~15 GB) — this will take a while..."
    echo "     If interrupted, re-run the script: curl will resume from where it stopped."
    curl -L --progress-bar -C - -o "$VEP_CACHE_TAR" "$VEP_CACHE_URL"
    echo "[OK] Download complete: $VEP_CACHE_TAR"
  else
    echo "[SKIP] VEP cache tarball already downloaded: $VEP_CACHE_TAR"
  fi

  echo "[RUN] Extracting VEP cache into: ${REF_ROOT}"
  tar -xzf "$VEP_CACHE_TAR" -C "${REF_ROOT}"
  echo "[OK] VEP cache installed: $VEP_CACHE_DIR"

  # Remove tarball to free space
  rm -f "$VEP_CACHE_TAR"
  echo "[OK] Tarball removed to free disk space."
fi

# -------------------------
# STEP 9: VEP plugins
# -------------------------
log_step "[STEP 9] VEP Plugins"

PLUGINS_DIR="${REF_ROOT}/Plugins"

if [[ -d "${PLUGINS_DIR}/.git" ]]; then
  echo "[SKIP] VEP plugins repo already cloned: $PLUGINS_DIR"
else
  if [[ -d "$PLUGINS_DIR" ]]; then
    echo "[WARN] $PLUGINS_DIR exists but is not a git repo."
    echo "       Remove the directory manually if you want to re-clone, then re-run."
  else
    echo "[RUN] Cloning Ensembl VEP plugins repository..."
    git clone https://github.com/Ensembl/VEP_plugins.git "$PLUGINS_DIR"
    echo "[OK] VEP plugins cloned: $PLUGINS_DIR"
  fi
fi

# -------------------------
# STEP 10: Cancer Hotspots
# -------------------------
log_step "[STEP 10] Cancer Hotspots (hg38)"

mkdir_if_missing "${REF_ROOT}/CancerHotSpots"

HOTSPOT_VCF="${REF_ROOT}/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz"

download_if_missing \
  "https://raw.githubusercontent.com/charlottekyng/cancer_hotspots/master/resources/hg38.hotspots_changv2_gao_nc.vcf.gz" \
  "$HOTSPOT_VCF"

HOTSPOT_TBI="${HOTSPOT_VCF}.tbi"
if [[ -f "$HOTSPOT_TBI" ]]; then
  echo "[SKIP] CancerHotSpots tabix index already exists: $HOTSPOT_TBI"
else
  echo "[RUN] Indexing CancerHotSpots with tabix"
  $TABIX -p vcf "$HOTSPOT_VCF"
  echo "[OK] CancerHotSpots indexed: $HOTSPOT_TBI"
fi

# -------------------------
# STEP 11: Write config file
# -------------------------
log_step "[STEP 11] Writing pipeline configuration file"

CONFIG_FILE="pipeline_config.sh"

# Resource versions / metadata
INSTALL_DATE="$(date +%Y-%m-%d)"
REFERENCE_VERSION="UCSC_hg38"
CHAIN_VERSION="UCSC_hg19ToHg38"
CIVIC_VERSION="civic_nightly_${CIVIC_DOWNLOAD_DATE}"
REVEL_VERSION="v1.3"
SPLICEAI_VERSION="SpliceAIv1.3"
LOEUF_VERSION="custom_grch38"
CANCER_HOTSPOTS_VERSION="charlottekyng/cancer_hotspots_master"
VEP_CACHE_VERSION="GRCh38"
VEP_PLUGINS_VERSION="git_clone_$(date +%Y-%m-%d)"

CLINVAR_VERSION="$(zgrep -m1 '^##fileDate=' "$CLINVAR_VCF" 2>/dev/null | cut -d= -f2 || true)"
if [[ -z "$CLINVAR_VERSION" ]]; then
  CLINVAR_VERSION="unknown_check_${INSTALL_DATE}"
fi

cat > "$CONFIG_FILE" <<EOF
# --------------------------------------------------
# Clinical Genomics Pipeline configuration
# Generated automatically by install_references.sh
# --------------------------------------------------
#
# Installed by : $(whoami)
# Install date : ${INSTALL_DATE}
# Install host : $(hostname)
# --------------------------------------------------

INSTALL_DATE="${INSTALL_DATE}"
REF_ROOT="$(realpath "$REF_ROOT")"

REF_GENOME="\${REF_ROOT}/reference/hg38.fa"
REF_GENOME_VERSION="${REFERENCE_VERSION}"

CHAIN_FILE="\${REF_ROOT}/liftover/hg19ToHg38.over.chain"
CHAIN_FILE_VERSION="${CHAIN_VERSION}"

CLINVAR_VCF="\${REF_ROOT}/ClinVar/clinvar.vcf.gz"
CLINVAR_VERSION="${CLINVAR_VERSION}"

CIVIC_FILE="\${REF_ROOT}/CIVIC/nightly-VariantSummaries.tsv"
CIVIC_VCF="\${REF_ROOT}/CIVIC/civic_grch38.vcf.gz"
CIVIC_VCF_INDEX="\${REF_ROOT}/CIVIC/civic_grch38.vcf.gz.tbi"
CIVIC_VERSION="${CIVIC_VERSION}"

REVEL_FILE="\${REF_ROOT}/REVEL/revel_all_chromosomes.zip"
REVEL_VERSION="${REVEL_VERSION}"

SPLICEAI_SNV="\${REF_ROOT}/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="\${REF_ROOT}/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"
SPLICEAI_VERSION="${SPLICEAI_VERSION}"

LOEUF_FILE="\${REF_ROOT}/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz"
LOEUF_FILE_INDEX="\${REF_ROOT}/gnomAD_constraints/loeuf_dataset_grch38.tsv.gz.tbi"
LOEUF_VERSION="${LOEUF_VERSION}"

VEP_CACHE_DIR="\${REF_ROOT}/homo_sapiens"
VEP_CACHE_VERSION="${VEP_CACHE_VERSION}"

VEP_PLUGINS_DIR="\${REF_ROOT}/Plugins"
VEP_PLUGINS_VERSION="${VEP_PLUGINS_VERSION}"

CANCER_HOTSPOTS_VCF="\${REF_ROOT}/CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz"
CANCER_HOTSPOTS_VERSION="${CANCER_HOTSPOTS_VERSION}"
EOF

echo "[OK] Configuration file created: $CONFIG_FILE"

# -------------------------
# Summary
# -------------------------
echo ""
echo "================================================================================"
echo " DONE"
echo "================================================================================"
echo " Installed by   : $(whoami)"
echo " Install root   : $(realpath "$REF_ROOT")"
echo " Config file    : $(realpath "$CONFIG_FILE")"
echo ""
echo " If a step failed, re-run this script — completed steps will be skipped."
echo "================================================================================"


###############################################################################
# SpliceAI Download Script (HPC)
# This is what I did in case this helps
# Following commands do:
#   1. Installs BaseSpace CLI
#   2. Authenticates (interactive)
#   3. Lists project and dataset
#   4. Downloads dataset using tmux (safe for SSH disconnects)
#
# Usage:
#   chmod +x download_spliceai.sh
#   ./download_spliceai.sh
###############################################################################

# set -e

# ###############################################################################
# # CONFIGURATION
# ###############################################################################

# PROJECT_ID="66029966"
# DATASET_ID="ds.20a701bc58ab45b59de2576db79ac8d0"
# OUTPUT_DIR="${HOME}/spliceai_staging"
# TMUX_SESSION="spliceai_download"

# ###############################################################################
# # STEP 1: Install BaseSpace CLI (if not already installed)
# ###############################################################################

# echo "Installing BaseSpace CLI..."

# # Create local bin directory
# mkdir -p "${HOME}/bin"

# # Download CLI only if not already present
# if [ ! -f "${HOME}/bin/bs" ]; then
#     wget https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs \
#         -O "${HOME}/bin/bs"
#     chmod +x "${HOME}/bin/bs"
# fi

# # Add to PATH if needed
# if [[ ":$PATH:" != *":${HOME}/bin:"* ]]; then
#     export PATH="${HOME}/bin:$PATH"
#     echo 'export PATH="${HOME}/bin:$PATH"' >> "${HOME}/.bashrc"
# fi

# ###############################################################################
# # STEP 2: Authenticate (interactive step)
# ###############################################################################

# echo "Authenticating with BaseSpace..."
# echo "Follow the URL printed below in your browser, then press Enter."

# # This will prompt a URL → open it in your local browser
# bs auth || true

# ###############################################################################
# # STEP 3: Verify project and dataset (optional but useful)
# ###############################################################################

# echo "Checking project..."
# bs list projects | grep "${PROJECT_ID}" || echo "Project not found"

# echo "Checking datasets..."
# bs list datasets --project-id "${PROJECT_ID}" || true

# ###############################################################################
# # STEP 4: Download using tmux (prevents job interruption)
# ###############################################################################

# echo "Starting download in tmux session: ${TMUX_SESSION}"

# # Create output directory
# mkdir -p "${OUTPUT_DIR}"

# # If tmux session already exists → attach
# if tmux has-session -t "${TMUX_SESSION}" 2>/dev/null; then
#     echo "Session already exists. Attaching..."
#     tmux attach -t "${TMUX_SESSION}"
# else
#     # Start new tmux session in background
#     tmux new -d -s "${TMUX_SESSION}" "
#         echo 'Downloading SpliceAI dataset...';
#         bs dataset download \
#             --id ${DATASET_ID} \
#             -o ${OUTPUT_DIR};
#         echo 'Download complete.';
#         bash
#     "

#     echo "Download started in tmux."
#     echo "To monitor progress:"
#     echo "  tmux attach -t ${TMUX_SESSION}"
# fi
# Dont forget to tabix them! tabix -p vcf

# ###############################################################################
# # END
# ###############################################################################

# echo "Script finished."
