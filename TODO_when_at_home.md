# TODO — Steps requiring the External SSD (reference files)

These tasks were written while away from the machine with the reference data.
Complete them when connected to `/Volumes/ExternalSSD/refs/`.

---

## 1. CancerHotspots counts file — O7 implementation

**Context:** Added `add_cancerhotspots_counts()` to `post_analysis.py` to enable SVIG-UK
O7 scoring. The function looks for a JSON file with per-position, per-AA-change tumour counts
from the CancerHotspots database.

**What to do:**

```bash
# Copy the pre-downloaded file to the reference directory
cp Files4ThisProject/cancerhotspots_counts.json \
   /Volumes/ExternalSSD/refs/CancerHotSpots/cancerhotspots_counts.json

# Add to get_ref_files.sh (already has a placeholder comment — see utils/get_ref_files.sh)
# The line to add downloads fresh data from the API:
#   curl -s https://www.cancerhotspots.org/api/hotspots/single \
#        > "${REF_ROOT}/CancerHotSpots/cancerhotspots_counts.json"
```

**Wire into entrypoint.sh:**
The `post_analysis.py` step in `entrypoint.sh` needs `--cancerhotspots-counts` added:

```bash
# In the post_analysis.py call block, add:
--cancerhotspots-counts "$REF_DIR/CancerHotSpots/cancerhotspots_counts.json"
```

**Test:**

```bash
# Run verification1 and check two new columns appear in Tier 2 output:
#   CancerHotspots_position_count
#   CancerHotspots_aa_change_count
#
# Expected values for verification1 variants:
#   NRAS Q61R  → position_count=422, aa_change_count=204  → O7 Strong [+4]
#   IDH1 R132H → position_count=~3000+, aa_change_count=~2800+  → O7 Strong [+4]
#   BRAF V600E → position_count=~10000+, aa_change_count=~9000+  → O7 Strong [+4]

python3 tests/verification5/test_post_analysis.py  # unit tests should pass without SSD
bash tests/verification1/run_verification1.sh       # needs SSD + Docker
```

---

## 2. GENIE lookup — O4 implementation

**Context:** `add_genie_counts()` added to `post_analysis.py` and `utils/build_genie_lookup.py`
written. The function looks for a pre-processed lookup file `genie_lookup.tsv.gz`.

**Dataset:** AACR Project GENIE v19.0-public
**Size:** 271,837 samples / 227,696 patients
**Mutations file:** `data_mutations.txt`
**Source:** https://www.synapse.org/genie  (Synapse ID: syn7222066)
**Requires:** Free Synapse account + accept GENIE data use agreement

**Step A — Download GENIE v19.0-public:**
```bash
# 1. Create a free account at https://www.synapse.org
# 2. Accept the GENIE data use agreement at https://www.synapse.org/genie
# 3. Install the Synapse client and download:
pip install synapseclient

synapse login   # enter your Synapse credentials

# Navigate to the v19.0-public release folder and download data_mutations.txt
# Find the exact syn ID by browsing syn7222066 → "19.0-public" folder
synapse get <syn_id_for_v19_data_mutations.txt>

# The file will be named: data_mutations.txt (~1-2 GB uncompressed)
```

**Step B — Build the lookup file (run once, takes ~5–10 min):**
```bash
mkdir -p /Volumes/ExternalSSD/refs/GENIE

python3 utils/build_genie_lookup.py \
    --input  data_mutations.txt \
    --output /Volumes/ExternalSSD/refs/GENIE/genie_lookup.tsv.gz
```
Expected output: ~10–20 MB compressed (gene:protein_change → deduplicated patient count).

**Step C — Wire into entrypoint.sh:**
In the `post_analysis.py` call block, add:
```bash
--genie-counts "$REF_DIR/GENIE/genie_lookup.tsv.gz"
```

**Step D — Test:**
```bash
# Run verification1 and check GENIE_count column in Tier 2 output.
# Expected counts for known drivers in GENIE v19 (227k patients):
#   BRAF V600E  → ~20,000+
#   KRAS G12D   → ~15,000+
#   NRAS Q61R   → ~2,000+
#   IDH1 R132H  → ~3,000+

bash tests/verification1/run_verification1.sh
```

---

## 3. snpEff hg19 database — Azurify dependency

**Context:** Azurify (ML classifier) requires snpEff with the hg19 genome database.
The download was interrupted by a network outage.

**What to do:**

```bash
export PATH="/opt/homebrew/opt/openjdk/bin:$PATH"
java -jar /tmp/snpEff/snpEff.jar download hg19
# Takes ~10-15 min, requires ~3 GB disk space

# Verify:
ls /tmp/snpEff/data/hg19/
```

Then test Azurify integration:

```bash
python3 utils/run_azurify.py \
    --smart-output tests/verification1/output/output \
    --azurify-dir  /tmp/Azurify_repo \
    --snpeff-jar   /tmp/snpEff/snpEff.jar
```

---

## 3. Docker image rebuild (after code changes)

After completing items 1 and 2, rebuild and push the Docker image to pick up
the `post_analysis.py` changes:

```bash
docker buildx build \
  --platform linux/amd64,linux/arm64 \
  --tag monkiky/smart:latest \
  --push .
```

---

*Written: May 2026. All code changes are committed on the `main` branch.*
