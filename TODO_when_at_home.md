# TODO — Steps requiring the External SSD (reference files)

---

## ✅ 1. CancerHotspots counts file — DONE

Files copied to SSD:
- `/Volumes/ExternalSSD/refs/CancerHotSpots/cancerhotspots_counts.json`
- `/Volumes/ExternalSSD/refs/oncokb_gene_roles.tsv`
- `/Volumes/ExternalSSD/refs/svig_uk_canonical_variants.tsv`

`entrypoint.sh` updated with all new `--` arguments.
`utils/get_ref_files.sh` updated to download/copy these files automatically.

---

## 🔶 2. GENIE lookup — O4 implementation

**Dataset:** AACR Project GENIE v19.0-public (271,837 samples / 227,696 patients)
**Source:** https://www.synapse.org/genie (Synapse ID: syn7222066)
**Requires:** Free Synapse account + accept GENIE data use agreement

**Step A — Download:**
```bash
pip install synapseclient
synapse login
synapse get <syn_id_for_v19_data_mutations.txt>
# file: data_mutations.txt (~1-2 GB)
```

**Step B — Build lookup:**
```bash
mkdir -p /Volumes/ExternalSSD/refs/GENIE
python3 utils/build_genie_lookup.py \
    --input  data_mutations.txt \
    --output /Volumes/ExternalSSD/refs/GENIE/genie_lookup.tsv.gz
```
Expected: ~10–20 MB compressed, ~5–10 min runtime.

**Step C — Already wired** into `entrypoint.sh` (points to `$REF_DIR/GENIE/genie_lookup.tsv.gz`).
File just needs to exist at that path.

**Step D — Test:**
```bash
bash tests/verification1/run_verification1.sh
# Check GENIE_count column in Tier 2 output:
#   BRAF V600E → ~20,000+  |  KRAS G12D → ~15,000+  |  NRAS Q61R → ~2,000+
```

---

## 🔶 3. snpEff hg19 database — Azurify dependency

**Status:** Download started in background. Check:
```bash
ls /tmp/snpEff/data/hg19/
```

If not present, re-run:
```bash
export PATH="/opt/homebrew/opt/openjdk/bin:$PATH"
java -jar /tmp/snpEff/snpEff.jar download hg19
```

Then test Azurify:
```bash
python3 utils/run_azurify.py \
    --smart-output tests/verification1/output/output \
    --azurify-dir  /tmp/Azurify_repo \
    --snpeff-jar   /tmp/snpEff/snpEff.jar
```

---

## 🔶 4. Run verification1 end-to-end

After all above, run the full pipeline to validate all new SVIG-UK columns:
```bash
export ONCOKB_TOKEN=your_token_here
bash tests/run_all_verifications.sh --only verification1
```

Check these columns in `Final_result_tier2.tsv`:
- `SVIG_UK_classification` (also in Tier 3)
- `SVIG_UK_score`
- `SVIG_UK_codes`
- `CancerHotspots_position_count` / `CancerHotspots_aa_change_count`
- `Gene_role_in_cancer`
- `O1_canonical`
- `O5_same_position`
- `O9_candidate`
- `GENIE_count` (if Step 2 complete)

---

## 🔶 5. Docker rebuild + push

```bash
docker buildx build \
  --platform linux/amd64,linux/arm64 \
  --tag monkiky/smart:latest \
  --push .
```

*Updated: May 2026.*
