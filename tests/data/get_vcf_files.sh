#!/bin/bash
# get_files2test_others_VC.sh
# Subsamples ~1% of PASS variants from standard human chromosomes
# (chr1-22, chrX, chrY) from VCF files in the repo root.
# Outputs bgzip-compressed subsamples to tests/data/.
# Deletes original uncompressed VCFs from the repo root after processing.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="${SCRIPT_DIR}"
OUTPUT_DIR="${SCRIPT_DIR}/data"
HTSLIB_IMAGE="quay.io/biocontainers/htslib:1.21--h566b1c6_1"

mkdir -p "${OUTPUT_DIR}"

# bgzip/tabix via Docker (aliases not available in non-interactive shells)
bgzip() { docker run --rm -v "$(pwd):/data" -w /data "${HTSLIB_IMAGE}" bgzip "$@"; }
tabix() { docker run --rm -v "$(pwd):/data" -w /data "${HTSLIB_IMAGE}" tabix "$@"; }

BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/analysis/SNVs/vcfs/FFG"

# ---------------------------------------------------------------------------
# 1. Download  (comment out if files are already present)
# ---------------------------------------------------------------------------
# cd "${SOURCE_DIR}"
# MuTect2 VCFs
#  curl -O "${BASE}/FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.muTect2.vcf.gz"
# curl -O "${BASE}/FFG_GZ_T_24h-F_VS_WGS_IL_N_1.bwa.muTect2.vcf.gz"
# curl -O "${BASE}/FFG_IL_T_24h_VS_WGS_IL_N_1.bwa.muTect2.vcf.gz"
# SomaticSniper VCFs
#  curl -O "${BASE}/FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.somaticSniper.vcf.gz"
# curl -O "${BASE}/FFG_GZ_T_24h-F_VS_WGS_IL_N_1.bwa.somaticSniper.vcf.gz"
# curl -O "${BASE}/FFG_IL_T_24h_VS_WGS_IL_N_1.bwa.somaticSniper.vcf.gz"
# Strelka VCFs
 curl -O "${BASE}/FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.strelka.vcf.gz"
# curl -O "${BASE}/FFG_GZ_T_24h-F_VS_WGS_IL_N_1.bwa.strelka.vcf.gz"
# curl -O "${BASE}/FFG_IL_T_24h_VS_WGS_IL_N_1.bwa.strelka.vcf.gz"

# ---------------------------------------------------------------------------
# 2. Decompress any remaining .vcf.gz files in SOURCE_DIR
# ---------------------------------------------------------------------------
for gz in "${SOURCE_DIR}"/*.vcf.gz; do
    [ -f "$gz" ] || continue
    echo "Decompressing: $(basename "$gz")"
    gunzip "$gz"
done

# ---------------------------------------------------------------------------
# 3. Subsample, compress output, delete originals
# ---------------------------------------------------------------------------

# Exclude already-subsampled files so we don't re-process them
for vcf in "${SOURCE_DIR}"/*.vcf; do
    [ -f "$vcf" ] || continue
    [[ "$vcf" == *.subsample.vcf ]] && continue

    base="$(basename "${vcf%.vcf}")"
    out_vcf="${OUTPUT_DIR}/${base}.subsample.vcf"
    out_gz="${out_vcf}.gz"

    echo "════════════════════════════════════════"
    echo "Processing: $(basename "$vcf")"

    # Header
    grep "^#" "$vcf" > "$out_vcf"

    # All variants on standard human chromosomes only (chr1-22, chrX, chrY)
    PASS_variants=$(grep -v "^#" "$vcf" \
        | awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$/')

    total=$(echo "$PASS_variants" | grep -c "." || true)

    if [[ "$total" -eq 0 ]]; then
        echo "  No variants on standard chromosomes — skipping."
        rm -f "$out_vcf"
        rm -f "$vcf"
        continue
    fi

    # 1% random subsample
    echo "$PASS_variants" \
        | awk -v seed="$RANDOM" 'BEGIN { srand(seed) } rand() < 0.01 { print }' \
        | sort -k1,1V -k2,2n \
        >> "$out_vcf"

    sampled=$(grep -vc "^#" "$out_vcf" || true)
    chroms=$(grep -v "^#" "$out_vcf" | cut -f1 | sort -u | tr '\n' ' ')
    echo "  Total variants on standard chrs : $total"
    echo "  After 1% subsample              : $sampled"
    echo "  Chromosomes present             : $chroms"

    # bgzip-compress the subsample and index it
    echo "  Compressing → ${out_gz}"
    (cd "${OUTPUT_DIR}" && bgzip -f "$(basename "$out_vcf")")
    (cd "${OUTPUT_DIR}" && tabix -p vcf "$(basename "$out_gz")")
    echo "  Indexed     → ${out_gz}.tbi"

    # Delete original full uncompressed VCF
    rm -f "$vcf"
    echo "  Deleted original: $(basename "$vcf")"
done

# ---------------------------------------------------------------------------
# 4. Move and compress any subsample VCFs already in SOURCE_DIR
# ---------------------------------------------------------------------------
for vcf in "${SOURCE_DIR}"/*.subsample.vcf; do
    [ -f "$vcf" ] || continue
    base="$(basename "$vcf")"
    out_gz="${OUTPUT_DIR}/${base}.gz"

    echo "════════════════════════════════════════"
    echo "Moving existing subsample: ${base}"
    mv "$vcf" "${OUTPUT_DIR}/${base}"

    echo "  Compressing → ${out_gz}"
    (cd "${OUTPUT_DIR}" && bgzip -f "${base}")
    (cd "${OUTPUT_DIR}" && tabix -p vcf "${base}.gz")
    echo "  Indexed     → ${out_gz}.tbi"
done

# ---------------------------------------------------------------------------
# 5. Summary
# ---------------------------------------------------------------------------
echo ""
echo "════════════════════════════════════════"
echo "Done. Output files in ${OUTPUT_DIR}:"
for gz in "${OUTPUT_DIR}"/*.subsample.vcf.gz; do
    [ -f "$gz" ] || continue
    count=$(gzip -dc "$gz" | grep -vc "^#" || true)
    echo "  $(basename "$gz") : ${count} variants"
done
