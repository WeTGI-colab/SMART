# Which annotations describe your transcript, and which cannot

SMART is used for clinical analysis, so an annotation that silently describes a different
transcript from the one being reported is a safety problem, not an inconvenience. This page
states, source by source, whether a column describes the transcript you supplied in
`--transcripts-file` — and where it does not, why not.

Nothing here is asserted from intent. Every row of the table below was read out of the code
that performs the annotation, and the checks described at the end run on every verification.

## The short answer

| | |
|---|---|
| Governed by your transcript | everything VEP derives from a transcript, the dbNSFP predictors, LOEUF, SpliceAI in the default `compute` mode, and OncoKB's protein-change query |
| Cannot be transcript-specific | ClinVar, CIViC, COSMIC, Cancer Hotspots, CADD, gnomAD, GENIE — all keyed on genomic position and allele |
| Transcript-specific in the source, not resolved to yours | REVEL, MutationTaster, and OncoKB when it falls back to a genomic query |

## How your transcript governs a run

VEP is run with `--everything`, so it annotates **every** transcript overlapping the variant
and writes one CSQ block per transcript. `vcf2table.py` then selects a single block:

1. a transcript pinned by `oncokb2.0.py`, if any;
2. otherwise `MANE_SELECT` / `MANE` / `MANE_PLUS_CLINICAL` when it is an `NM_` accession;
3. otherwise the RefSeq accession in `Feature`;
4. otherwise the `NM_` accession in `HGVSc`;

and keeps the block whose accession is in your list. Every column VEP derives per transcript —
`Consequence`, `IMPACT`, `HGVSc`, `HGVSp`, `EXON`, `INTRON`, `Protein_position`, `Amino_acids`,
`Codons`, `SIFT`, `PolyPhen`, `DOMAINS`, `BIOTYPE`, `STRAND`, `CANONICAL` — comes from that one
block, so it describes your transcript by construction.

`HGVSc` and `HGVSp` carry the transcript inside the value (`ENST00000269305.9:c.524G>A`), which
makes that claim self-evidencing rather than a promise: the value names the transcript it
describes, and a test asserts the two agree on every row of real output.

### The one failure path, and how it was made visible

The next three paragraphs describe the output **as it was before** the two provenance columns were
added, because that is what they exist to fix. The fix follows immediately after.

If no block matches your list, `vcf2table.py` falls back to the first CSQ block:

```python
if not nm_selected and first_vep is not None:
    # No match found — use the first CSQ record
```

That is an arbitrary transcript, chosen by VEP's output order, and nothing in the output said
so. It happens when the variant is in a gene your list does not cover, or when VEP did not
annotate the listed transcript at all.

**Every other guarantee on this page is downstream of this one.** A predictor faithfully
resolved to the selected transcript is still the wrong answer if the selected transcript was
not the one you asked for.

It was worse than silent. `first_nm` was initialised and never assigned, so a fallback row left
`NM_Transcript` **empty** — which reads as "no transcript available" rather than "a transcript
you did not ask for". The two readings call for opposite responses.

### The fix: two columns that state it, per row

That was the state of the output before this was addressed. It is now stated per row. Every row of
the Tier 2 and Tier 3 sheets carries:

| column | value |
|---|---|
| `TRANSCRIPT_REQUESTED` | the accession from your list that was used, or **empty** when none of your transcripts was annotated for this variant |
| `TRANSCRIPT_MATCH` | `match` · `fallback:<acc>` |

- **`match`** — the selected transcript is one you supplied. The guarantee holds. This covers
  both routes to it: the whitelist walk, and `oncokb2.0.py`'s pin, which sets
  `ONCOKB_PREFERRED_TRANSCRIPT` only when the candidate is in the supplied list and therefore
  can never pin a transcript you did not ask for.
- **`fallback:<acc>`** — none of your transcripts was annotated for this variant, so VEP's
  first block was used. **Every transcript-specific column on that row describes `<acc>`**, not
  a transcript you requested. Read that row as you would an annotation from an unknown
  transcript.

Two values and no third. A first version of these columns reported the pin as its own outcome,
`pinned:<acc>`, which was a mistake worth recording: on the verification cohort as it then stood —
37 rows — it marked 32 of them as something other than a match when the pinned transcript was in
every case the supplied one, and a genuine fallback would have been indistinguishable in that
crowd.

So the question "did this row honour my transcript?" is answered in the row, not in a separate
report. The provenance report still checks it in aggregate and names any affected rows. Rather than
quoting a figure that moves with the cohort, read it off the output:

```bash
python3 scripts/transcript_provenance.py \
    --maf  tests/verification1/output/output/Final_result_tier1.maf \
    --transcripts tests/verification1/verification1_transcripts.txt
```

On the committed `verification1` output that is **45 of 45 rows on a supplied transcript** — every
row, no fallback. `tests/verification5/test_docs_claims.py` asserts that this page's figure and the
script's answer still agree.

## Source by source

| Source | Transcript-specific? | How it is keyed, and why |
|---|---|---|
| VEP core fields | **Yes** | one CSQ block per transcript; the selected block is yours |
| `HGVSc`, `HGVSp` | **Yes, provably** | the transcript is inside the value |
| LOEUF | **Yes, and it matters** | plugin run with `match_by=transcript`. LOEUF is widely described as a gene-level metric; in the dataset it is per transcript, and the values differ. Of the 14,247 genes with more than one transcript in `loeuf_dataset_grch38.tsv.gz`, **13,236 (93%) have a different LOEUF on different transcripts** — KAT6B is 0.057 on `ENST00000650232` (CDS 2,974 bp, 14 coding exons) and 1.891 on `ENST00000649119` (CDS 23 bp, 1 exon), because a short transcript has few possible loss-of-function sites, so the expected count is tiny and the confidence interval's upper bound is large. Looking LOEUF up by gene would therefore give an arbitrary transcript's value and be wrong 93% of the time, sometimes by the distance between "among the most constrained genes" and "no constraint signal at all" |
| dbNSFP predictors | **Yes** (default) | plugin run with `transcript_match=1`. Without it the cell held one value per transcript and marked none — see below |
| SpliceAI, `compute` mode | **Yes** | the model is run against an annotation built from your `--transcripts-file`, and `SpliceAI_pred_SYMBOL` records which transcript each score was computed on (`MET--NM_000245.4`) |
| OncoKB, protein-change query | **Yes** | queried with the selected transcript's `HGVSp_Short` |
| ClinVar | **No** — cannot be | ClinVar classifies a *variant*, not a variant-on-a-transcript. Matched on CHROM+POS+REF+ALT |
| CIViC | **No** — cannot be | evidence is curated per variant; matched on position and allele |
| COSMIC | **No** — cannot be | counts of observations of a genomic change |
| Cancer Hotspots | **No** — cannot be | a recurrently mutated residue or position, matched on position and allele |
| CADD | **No** — by design | a genome-wide deleteriousness model scoring a genomic substitution, with no transcript input |
| gnomAD | **No** — cannot be | an allele frequency in a population; a property of the allele |
| GENIE | **No** — cannot be | counts of patients carrying a genomic change |
| SpliceAI, `lookup` mode | **No** | Illumina's precomputed tables carry no transcript at all: one score per gene, computed on GENCODE v24 canonical of December 2015. This is why `compute` is the default ([#2](https://github.com/WeTGI-colab/SMART/issues/2)) |
| REVEL | **Not guaranteed** | see below |
| MutationTaster | **Not guaranteed** | see below |
| OncoKB, genomic-change fallback | **No** | see below |

"Cannot be" is not a euphemism. A ClinVar classification, a population frequency and a count
of observed cases are statements about a position or an allele; there is no transcript-specific
version of them to look up. Reporting them alongside transcript-specific columns is correct;
what would be wrong is implying they were resolved to your transcript.

## The three that are transcript-specific and still not resolved

### REVEL

REVEL holds one score per transcript. The plugin is run with `no_match=1`, which returns a
score even when the transcript does not match, and the output records no transcript to check
it against — so a REVEL value is not guaranteed to describe your transcript.

That flag was added deliberately, and removing it has a cost. REVEL's file keys on older
Ensembl identifiers: it holds 0.931 for BRAF V600E on `ENST00000288602` while the pipeline
reports `ENST00000646891`. Without `no_match=1` the column was empty for 9 of 71 missense rows
in a 93-variant cohort — every BRAF among them — and the SVIG-UK code O5, worth up to +4, was
silently unavailable for those variants. So the choice is between a score that may come from a
neighbouring transcript and no score at all, and the pipeline currently takes the score.

Transcript-level REVEL scores are usually close, because REVEL is trained on amino-acid
substitutions and a given genomic change often produces the same substitution on related
transcripts. That is a reason to accept it, not a reason to stop stating it.

### MutationTaster

dbNSFP keys its MutationTaster fields on a different transcript set from the rest of the file,
so `transcript_match=1` has no effect on them — the plugin's own documentation says so. Those
two columns still hold every transcript's value, `&`-separated. The provenance report names
them on every run, and a test fails if any column other than MutationTaster is left in that
state.

### OncoKB's genomic-change fallback

When the transcript you report is not the one OncoKB curates on, SMART does not give up: it
re-queries by genomic change. The answer then comes from OncoKB's own canonical transcript, not
yours, and the MAF records that this happened:

```
ONCOKB_QUERY_TYPE  MUTATION_BY_GENOMIC_CHANGE
ONCOKB_MAF_QUERY   Genomic_Change (selected transcript differs from OncoKB...)
```

This is deliberate and, for a therapeutic annotation, right. `chr7:55191822 T>G` is EGFR L858R
whichever transcript names it; answering "Unknown" because your run called it `p.L813R` would
be a false negative on a licensed first-line therapy. The column is therefore not
transcript-specific **and should not be** — but the output says when the fallback was used, so
a reader can see that the verdict came from the genomic change rather than from the transcript
being reported. See [#25](https://github.com/WeTGI-colab/SMART/issues/25).

## What is checked, on every run

`scripts/transcript_provenance.py` runs as part of `tests/run_all_verifications.sh` and prints:

- how many rows were annotated on a transcript you actually supplied, and names any that were
  not;
- the count of columns in each category above;
- for dbNSFP, whether the transcript it matched is the one the row reports — `Ensembl_transcriptid`
  against `Feature`, which is a direct check that the predictors describe your transcript;
- any column still holding an unresolved per-transcript list, with the reason.

`tests/verification5/test_transcript_provenance.py` (24 tests) covers the classification, and
asserts on real output that `HGVSc` names the reported transcript and that dbNSFP matched it on
every row. The per-transcript columns are detected structurally, not from a hardcoded list, so
a new dbNSFP column is covered the day it appears.

## Changing the defaults

```
--dbnsfp-transcript-match 0|1   resolve dbNSFP to the reported transcript (1, default)
--spliceai compute|lookup       compute against your transcripts (default) or look up the
                                precomputed tables, which carry no transcript
```

Both default to the transcript-specific behaviour. `0` and `lookup` exist for reproducing older
results, and each is recorded in the run header so an output can be traced to how it was made.
