# SMART Pipeline — Variant Caller Compatibility Verification 2

## Purpose

This verification set checks that the SMART annotation pipeline can process
VCF files produced by variant callers other than DRAGEN. Three tumour/normal
WGS VCFs from the SEQC2 benchmark sample (FFG_GZ_T_24h-B vs WGS_IL_N_1)
plus one tumour-only amplicon VCF from Pisces are used, all in GRCh38/hg38
coordinates.

## Variant callers tested

| VCF file | Caller | Sample type | FILTER field | PASS variants | Notes |
|---|---|---|---|---|---|
| `FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.muTect2.subsample.vcf.gz` | GATK MuTect2 | Tumour/Normal WGS | `t_lod_fstar` / `PASS` | 1,356 | TUMOR column before NORMAL |
| `FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.strelka.subsample.vcf.gz` | Illumina Strelka | Tumour/Normal WGS | `LowEVS` / `PASS` | 1,455 | NORMAL column before TUMOR |
| `FFG_GZ_T_24h-B_VS_WGS_IL_N_1.bwa.somaticSniper.subsample.vcf.gz` | SomaticSniper | Tumour/Normal WGS | `.` (no PASS) | 0 | Requires `--no-pass` flag |
| `To_test_PISCIS.vcf.gz` | Pisces 5.2.10 | Tumour-only amplicon | `q30` / `LowDP` / `PASS` | 4 | Single sample, `VF` instead of `AF` |

## What to check in the output

For each caller verify:

- **Pipeline completes** without errors
- **Tier 1 MAF produced** and passes `tests/verify_maf.py` structural checks
- **VEP annotation** assigns gene, consequence, and HGVS protein change
- **OncoKB annotation** runs correctly regardless of FORMAT structure
- **PASS filter** — MuTect2, Strelka, and Pisces annotate only PASS variants;
  SomaticSniper requires `--no-pass` because it reports `.` instead of `PASS`
- **Sample name** extracted correctly from each filename

## How to run

### MuTect2 and Strelka (PASS filter on, run together)

```bash
docker run --rm \
  -v /path/to/SMART/tests/verification2:/data \
  -v /path/to/SMART/tests/verification2/output:/output \
  -v /Volumes/ExternalSSD/refs:/refs:ro \
  monkiky/smart:1.1.0 \
  "$ONCOKB_TOKEN" \
  --ref-dir /refs \
  --config /data/Config.yaml \
  --input-dir /data \
  --output-dir /output \
  --no-liftover \
  --keep-tmp 2>&1 | tee /path/to/SMART/tests/verification2/output/smart.log
```

### SomaticSniper (PASS filter off)

```bash
docker run --rm \
  -v /path/to/SMART/tests/verification2:/data \
  -v /path/to/SMART/tests/verification2/output:/output \
  -v /Volumes/ExternalSSD/refs:/refs:ro \
  monkiky/smart:1.1.0 \
  "$ONCOKB_TOKEN" \
  --ref-dir /refs \
  --config /data/Config.yaml \
  --input-dir /data \
  --output-dir /output \
  --no-liftover \
  --no-pass \
  --keep-tmp 2>&1 | tee /path/to/SMART/tests/verification2/output/smart_nopass.log
```

### Verify MAF output

```bash
for maf in /path/to/SMART/tests/verification2/output/output/*.maf; do
  echo "--- $maf ---"
  python3 /path/to/SMART/tests/verify_maf.py "$maf"
done
```
