# Annotation flags: keeping rows that a filter would drop

`SPLICEAI_MAX_DS` and `ANNOTATION_FLAGS` are two columns of the annotation output. They exist so
that a variant worth reporting is not removed by a filter applied before anyone reads the evidence
that makes it worth reporting.

Requested in, and specified by, two issues that turned out to describe the same problem:
[#14](https://github.com/WeTGI-colab/SMART/issues/14) (synonymous and intronic variants with strong
SpliceAI evidence) and [#15](https://github.com/WeTGI-colab/SMART/issues/15) (the two TERT promoter
hotspot positions). They are solved with one mechanism rather than two special cases, so a further
reason to keep a row needs no further column.

## The problem they solve

SMART annotates SpliceAI and reaches the output with it: `SpliceAI_pred_DS_AG`, `_DS_AL`, `_DS_DG`
and `_DS_DL` are columns 1,093–1,096 of the 1,165 in the Tier 1 MAF of the committed
`tests/verification1` run, and the same four are in the Tier 3 sheet. The two recurrent TERT
promoter variants (C228T and C250T) are annotated correctly too, as `upstream_gene_variant`, which
SMART maps to `Variant_Classification = 5'Flank`.

Neither is missing. What happens is that a filter on `Variant_Classification` or on VEP's `IMPACT`
removes the row — as `Silent`, `Intron` or `5'Flank`, or as `LOW`/`MODIFIER` — before anyone reaches
column 1,093. The data is there; the row is gone.

(The exact column positions move with the cohort: the OncoKB evidence columns are padded to the
widest row in a run. The point is that they are a thousand columns to the right of the ones a
filter reads, not that they are at any particular index.)

The volume is low and the stakes are not. Across three cohorts (10,124 rows) roughly one variant in
five is synonymous, intronic or UTR; of those, 11 carry a SpliceAI delta of 0.5 or above and about
25 reach 0.2. Those are intronic changes three to five bases into the splice region with deltas
between 0.90 and 0.99, in genes where loss of function matters. The TERT promoter variants are among
the most frequently reported non-coding variants in oncology.

## The two columns

| Column | Contents |
|---|---|
| `SPLICEAI_MAX_DS` | the highest of the four SpliceAI delta scores, as a number, so splice evidence can be filtered and sorted on one value. Empty where SpliceAI returned no prediction. Computed for every row, whatever its impact. |
| `ANNOTATION_FLAGS` | reasons to keep the row, `;`-separated. Empty where none apply — the column is always present, so a filter can rely on it. |

Flags currently emitted:

| Flag | Meaning |
|---|---|
| `splice_predicted(<score>)` | a `LOW` or `MODIFIER` variant whose `SPLICEAI_MAX_DS` reaches 0.2 ([#14](https://github.com/WeTGI-colab/SMART/issues/14)) |
| `svig_canonical` | a position on the SVIG-UK canonical variants list; today the two TERT promoter hotspots ([#15](https://github.com/WeTGI-colab/SMART/issues/15)) |

## How to use them

The filter changes from

```
drop rows where Variant_Classification in (Silent, Intron, 5'Flank)
```

to

```
drop rows where Variant_Classification in (Silent, Intron, 5'Flank)
                and ANNOTATION_FLAGS is empty
```

or, for an impact-based filter, keep `HIGH` and `MODERATE` **or** any row with a non-empty
`ANNOTATION_FLAGS`. One line on the consuming side, and one column carries every reason, so a
further reason to keep a row does not add a further column.

A laboratory that wants a stricter splice cut than 0.2 filters on `SPLICEAI_MAX_DS` itself; the
score is in the sheet for every row regardless of impact.

## Why the threshold is 0.2

(The threshold, the scope and the flag name were settled on the issues — see the discussion on
[#14](https://github.com/WeTGI-colab/SMART/issues/14#issuecomment-6059166063) and
[#15](https://github.com/WeTGI-colab/SMART/issues/15#issuecomment-6059416021).)

0.2 is the delta at which the classifier starts awarding O6. Rescuing from 0.5 — the value SMART's
`SpliceAI_cutoff` passes — would mean that everything between 0.2 and 0.5 is filtered out before the
classifier sees it, even though the classifier would have awarded evidence to it. Flagging from the
same value means every variant that could earn O6 reaches the classifier.

## Why only LOW and MODIFIER

Those are the impacts an impact-based filter removes, and they are where splice-disrupting variants
are lost. `HIGH` and `MODERATE` rows survive such a filter anyway, so flagging them would flag much
of the sheet without keeping a single extra row. Their scores are still in `SPLICEAI_MAX_DS`.

## Why `Variant_Classification` is not rewritten

Relabelling a rescued row as `Splice_Region` or `SpliceAI_Rescued` would be the shorter change, and
it is the wrong one (the point is argued on
[#14](https://github.com/WeTGI-colab/SMART/issues/14), where it is also reported from experience
that relabelling lost the original classification and needed a custom vocabulary list):

- it would no longer be true — the variant is still synonymous;
- `Variant_Classification` is a defined MAF vocabulary, and a value outside it breaks third-party
  tools that read the MAF, and needs a custom vocabulary list in tools such as maftools;
- it destroys the record of what the consequence actually was, which is the thing a scientist needs
  in order to judge the prediction.

The flag adds the reason without removing the fact.

## Why the canonical positions are read from a reference file

`svig_canonical` is set for every entry of
[`reference_lists/svig_uk_canonical_variants.tsv`](../reference_lists/svig_uk_canonical_variants.tsv)
whose `match_rule` begins with `genomic:` — the same file, and the same entries, that the classifier
uses for O1. Writing `chr5:1295113` into the annotation code would create a second definition that
can drift from the first, and the set extends by itself if further non-coding positions are added to
the list.

Matching is genomic, not by HGVS. A promoter variant is described relative to a transcript, so
`c.-124C>T` depends on which TERT transcript was annotated; the genomic position does not. Entries
matched by protein change or by a bespoke rule are not flagged — a position is what an annotation
flag can test.

## Provenance

Both columns are computed by SMART, so nothing in any VCF describes them and their descriptions come
from `Config.yaml` — see [`field_provenance.md`](field_provenance.md) for how that division works.
Both are tier 3, so they reach the sheet clinical scientists read.

## Verified against a real MAF

The flags run inside the per-sample transform, after the annotation columns have been expanded and
before `drop_columns`. On the real per-sample MAF from the committed `tests/verification1` run,
`FINAL_Table/verification1.oncoKB.maf` — 45 rows, 1,151 columns — both columns are added,
`SPLICEAI_MAX_DS` parses the real cells (VEP writes them as a numeric column, so empty cells arrive
as `NaN` rather than as `.`), and both reach the Tier 2 and Tier 3 sheets.

`add_annotation_flags()` on that MAF flags **4 of the 45 rows**, one per route into the mechanism:

| Row | Impact / class | `SPLICEAI_MAX_DS` | Flag |
|---|---|---|---|
| TERT chr5:1295113 G>A | MODIFIER / `5'Flank` | — | `svig_canonical` |
| TERT chr5:1295135 G>A | MODIFIER / `5'Flank` | — | `svig_canonical` |
| MET chr7:116771992 | LOW / `Intron` | 0.58 | `splice_predicted(0.58)` |
| TP53 chr17:7675994 | LOW / `Splice_Region` | 0.58 | `splice_predicted(0.58)` |

Each is exactly the row a consequence- or impact-based filter drops: two promoter hotspots VEP calls
`5'Flank`, an intronic change and a protein-synonymous splice-region change. None of them would
survive the filter this mechanism exists to correct, and all four now carry a reason not to drop
them.

The two `LOW`/`MODIFIER` rows that are **not** flagged are as important, because they are where a
looser rule would produce false positives:

- a third TERT `5'Flank` row, chr5:1295113 **G>T**. Same position as the first hotspot, different
  allele, and only `G>A` is on the canonical list — matching on position alone would flag it.
- a CDK4 5'UTR row for which SpliceAI returned no prediction at all. An empty `SPLICEAI_MAX_DS` is
  not a low score, and is not evidence.

The sample used to contain no flaggable row at all, which is why earlier versions of this page said
so. Its only `Silent`/`LOW` row then was a GNAQ variant that was synonymous only because the
fixture's allele was strand-flipped; correcting that ([#19](https://github.com/WeTGI-colab/SMART/issues/19))
turned it into the intended `p.Q209L` missense, and the cohort has since gained the MET, TP53 and
TERT records above. **There is no `Silent` row in the MAF today.** The unit tests still construct
eligible rows synthetically, for the boundaries that a 45-variant cohort of coding hotspots does not
happen to sit on; the integration tests assert the four rows above.

## Tests

`tests/verification5/test_annotation_flags.py` (63 tests): the 0.2 boundary from both sides, the
`LOW`/`MODIFIER` scope, both TERT positions and the near misses that must not match, the two flags
on one row, the empty forms SpliceAI writes, and that `Variant_Classification` is unchanged and
unflagged rows are written empty rather than as the word `None`. The last two groups cover the
dtypes real data arrives in (`float64` with `NaN`, and an integer `Start_Position` matched against
a list that holds it as text) and the real MAF described above.
