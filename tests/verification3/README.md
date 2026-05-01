# SMART Pipeline — Verification 3: Transcript Prioritisation Impact

## Purpose

This verification demonstrates that the preferred-transcript whitelist meaningfully
changes the VEP annotation output, and that the pipeline consistently selects the
correct transcript for clinically actionable variants.

The pipeline is run twice on the same VCF using two different transcript files.
`verify.py` then confirms:

1. **DIFFER** — targeted variants produce different VEP annotations between runs
   (the transcript whitelist is actually driving annotation selection)
2. **PASS** — each run's annotation is confirmed correct by Ensembl REST API
   using that run's expected transcript

## Test Cases

| Variant | Transcript A (preferred) | Result A | Transcript B (alternative) | Result B |
|---|---|---|---|---|
| EGFR L858R (chr7:55191822 T>G) | NM_005228.5 | p.Leu858Arg (**LEVEL_1** in OncoKB) | NM_001346897.2 | p.Leu813Arg (Unknown in OncoKB) |
| TP53 R175H (chr17:7675088 C>T) | NM_000546.6 | p.Arg175His (hotspot, R175H) | NM_001126115.2 | p.Arg43His (different residue entirely) |
| CDKN2A/CDKN2B DEL (chr9:21990000-22005000) | NM_058195.4 | Gene = **CDKN2A** (LEVEL_4) | NM_004936.4 | Gene = **CDKN2B** (no treatment level) |
| FGFR1 N546K (chr8:38417331 G>T) | NM_023110.3 (MANE Select) | p.Asn546Lys (Unknown in OncoKB) | NM_001174067.2 (TSO500 non-MANE) | p.Asn577Lys (**Likely Oncogenic, LEVEL_4** — false positive) |

### Why this matters

- **EGFR L858R**: OncoKB annotation depends entirely on the protein change sent.
  NM_005228 gives the canonical p.Leu858Arg (FDA-approved targeted therapy at LEVEL_1).
  A wrong transcript gives p.Leu813Arg — an unrecognised variant with no therapeutic
  implication, even though the genomic position is identical.

- **TP53 R175H**: Multiple TP53 isoforms exist with different N-terminal lengths.
  NM_001126115 encodes a shorter isoform; the same codon change lands at position 43
  instead of 175, making the hotspot annotation entirely incorrect.

- **CDKN2A/CDKN2B deletion**: The deletion at chr9:21990000-22005000 spans a region
  that overlaps both CDKN2A (p14ARF isoform, NM_058195) and CDKN2B (NM_004936).
  Depending on which transcript is preferred, the pipeline reports a different gene
  as the affected target — a clinically significant distinction.

- **FGFR1 N546K**: A real-world TSO500 panel scenario. The TSO500 gene list includes
  NM_001174067.2 for FGFR1, which is not a MANE transcript; the MANE Select is
  NM_023110.3. The same genomic position (chr8:38417331 G>T) produces p.Asn546Lys
  via the MANE transcript (not currently classified by OncoKB) but p.Asn577Lys via
  the TSO500 non-MANE isoform, which OncoKB returns as Likely Oncogenic with LEVEL_4.
  This is a **false-positive** therapeutic signal arising solely from incorrect
  transcript selection — the most clinically dangerous failure mode.

## Transcript files

| File | Contents |
|---|---|
| `transcripts_A.txt` | NM_005228.5 (EGFR MANE), NM_000546.6 (TP53 MANE), NM_058195.4 (CDKN2A p14ARF), NM_023110.3 (FGFR1 MANE Select) |
| `transcripts_B.txt` | NM_001346897.2 (EGFR alt), NM_001126115.2 (TP53 short isoform), NM_004936.4 (CDKN2B), NM_001174067.2 (FGFR1 TSO500 non-MANE) |

## How to run

```bash
export ONCOKB_TOKEN="<your-token>"
bash tests/verification3/run_verification3.sh
```

This will:
1. Run the pipeline with `transcripts_A.txt` → output in `output_A/`
2. Run the pipeline with `transcripts_B.txt` → output in `output_B/`
3. Run `verify.py` comparing both MAF outputs

## Design rationale — why expected values are hardcoded

`verify.py` contains an `EXPECTED` dictionary with the annotation values predicted for
each variant under each transcript set. These are not guesses — they were looked up via
the Ensembl REST API and OncoKB API before being written into the code, and represent
real external ground truth frozen into the test.

The verification works in three layers:

1. **DIFFER** — confirms the two runs produce different output. "Different" alone is not
   enough: both runs could be wrong, just differently.

2. **VEP API check** — queries Ensembl REST independently to confirm each run's
   annotation is correct for its transcript. This catches pipeline bugs regardless of
   what the expected values say.

3. **OncoKB impact** — reads the OncoKB fields already written into the MAF by the
   pipeline. Without hardcoded expectations, you could only say *"the two runs differ"*
   — not *"run A is clinically correct and run B is a false positive"*. The FGFR1 case
   makes this essential: both runs produce non-empty, plausible-looking results; only
   the expected values let you assert that the LEVEL_4 result from run B is wrong.

If OncoKB updates their database and a result changes, the test will fail and prompt
review of the clinical interpretation — which is the intended behaviour.

## verify.py output

```
DIFFER  — annotations differ between runs (expected for targeted variants)
SAME    — annotations are identical between runs (unexpected — investigate)
PASS    — pipeline annotation matches Ensembl REST for the expected transcript
MISMATCH — pipeline annotation does not match Ensembl REST
ERROR   — variant missing from MAF or API error
```

Results are written to `results.tsv`.

Exit code `0` = all checks passed. Exit code `1` = any SAME, MISMATCH, or ERROR.
