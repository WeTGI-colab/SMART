# Testing SMART

SMART is an annotation pipeline, so almost every defect is of one kind: a value in the output
that does not mean what the column says it means. The tests are arranged around three questions
that find that kind of defect.

**Does the logic do what the guideline says?** Pure-Python unit tests, no Docker, no network.
Fast enough to run on every change.

**Does the pipeline still run end to end?** Containerised runs over small fixtures, checking that
VCFs from different callers are handled, that the transcript whitelist changes what it should,
and that parallel execution matches a sequential one.

**Does each annotated value match what its source actually says?** Every column is traced back to
the OncoKB, VEP, CIViC, COSMIC or CADD API, or to the file on disk the pipeline read it from.
This is the layer that catches the quiet failures, and most of this year's real findings came
from it.

## The layers

| | Where | Needs | What it answers |
|---|---|---|---|
| Unit tests | `verification5/` | python3 + pytest | the logic, 1,183 tests in about fifteen seconds |
| Caller compatibility | `verification2/` | Docker, refs | four callers' VCF conventions (manual) |
| Transcript choice | `verification3/` | Docker, refs, token | the whitelist changes the annotation and the OncoKB call |
| Parallel jobs | `verification4/` | Docker, refs, token | `--jobs N` matches a sequential run |
| API verification | `verification1/` | Docker, refs, token | **47 columns** against the sources' own APIs |
| Source verification | `verify_local_sources.py` | refs | a further **146 columns**, against the files on disk |
| Variant presence | `verify_maf.py` | — | expected variants are in an output |

### What the API verification cohort covers

`verification1` holds 35 records on 22 panel transcripts: 30 SNVs, one in-frame deletion and
four CNAs. Thirteen were added deliberately to cover what the suite's own field coverage report
said was untested, rather than for being well known:

- **Breadth of OncoKB evidence.** KRAS Q61H (15 treatments across levels 1/2/3A/4/R1),
  NRAS G12D (14, plus Dx2 and Px1), KIT D816V (LEVEL_R2 and 17 citations), PDGFRA D842V
  (LEVEL_1 with LEVEL_R1), IDH2 R140Q (two diagnostic and three prognostic implications) and
  TP53 Y220C (nine prognostic implications).
- **`ONCOGENIC` values nothing else produces.** ABL1 T315I is `Resistance`, MLH1 V384D is
  `Inconclusive`, STK11 F354L is `Likely Neutral`. Before these, only `Oncogenic` and
  `Likely Oncogenic` were ever exercised.
- **`ONCOKB_VUS`.** CDH1 R63* is the one variant found with `vus=True`.
- **Intronic and splice variants.** MET `c.3028+3A>G` (intron 14/20), TP53 `c.376-1G>A`
  (splice acceptor, intron 4/10) and TP53 `c.782+1G>A` (splice donor, intron 7/10). These fill
  the `INTRON` field, and they are the first records in this cohort that exercise SpliceAI,
  `ANNOTATION_FLAGS` and O6 at all — the cohort was otherwise coding hotspots, which is not
  where splice prediction matters.

**Every added record was confirmed three ways before going in**: its REF against the GRCh38
FASTA, and OncoKB's `byGenomicChange` endpoint returning the alteration the variant ID claims.
That third check is not ceremony. Ensembl's variant recoder offered `chr4:54733167 A>T` for
KIT D816V, which OncoKB resolves to **D820V** — a different variant — and the same happened for
IDH2 (R149Q) and KRAS (S65I). The VEP HGVS endpoint returned all three TP53 variants
**strand-flipped**, because it reports `allele_string` on the transcript strand and TP53 is on
the minus strand: the identical defect as [#19](https://github.com/WeTGI-colab/SMART/issues/19),
reappearing. Taking the first coordinate a service returns is how #19 happened.

Three OncoKB fields remain uncovered and **no variant can cover them**: `LEVEL_3B`,
`ONCOKB_diagnosticSummary` and `ONCOKB_prognosticSummary` all require a tumour-type query,
and the pipeline queries pan-cancer ([#26](https://github.com/WeTGI-colab/SMART/issues/26)).

## Running them

```bash
export ONCOKB_TOKEN="<your-token>"
bash tests/run_all_verifications.sh --refs /path/to/refs

# one at a time
bash tests/run_all_verifications.sh --only verification5    # no Docker, no token, no refs
bash tests/run_all_verifications.sh --only verification1
bash tests/run_all_verifications.sh --only localsources     # needs --refs
```

The unit tests on their own:

```bash
python3 -m pytest tests/ -q
```

Each suite has its own README with the detail; this file is the map.

`tests/verify_maf.py` is a standalone helper rather than a suite: it checks that specific known
variants are present in a MAF with the correct ID values, as a quick sanity check after a run.

`tests/e2e_test.py` and `tests/test_multi_transcript.py` run under pytest with the rest.

## A worked example: BRAF V600E

One variant, and what each layer asserts about it. This is a useful example precisely because
two of the answers are surprising.

The row, as the pipeline produces it:

```
Chromosome              chr7
Start_Position          140753336
Reference_Allele        A
Tumor_Seq_Allele2       T
Feature                 ENST00000646891        <- the transcript we report
HGVSp_Short             p.V600E
Variant_Classification  Missense_Mutation
```

**The unit tests** (`verification5`) never see this row. They assert that, given a variant with
these properties, each evidence code fires as the guideline text says — 211 tests over the
classifier alone, with their expectations taken from the guideline rather than from observed
output, so that a wrong implementation cannot make the test agree with it.

**Verification 1** asks OncoKB's API about this variant and compares field by field. The output
says `ONCOKB_ONCOGENIC = Oncogenic`; the API must say the same, and 506 of the 813 committed
checks are OncoKB fields like this one. It does the same against the VEP, CIViC, COSMIC and CADD
APIs. It cannot check gnomAD, ClinVar, SpliceAI, REVEL, LOEUF or Cancer Hotspots, because those
are annotated from local files and have no API — its own header says so.

**Source verification** covers those. For this variant:

```
gnomADe_AF = 1.369e-06     the VEP REST API returns gnomade = 1.369e-06            PASS
ClinVar_CLNSIG = Conflicting_classifications_of_pathogenicity
                           clinvar.vcf.gz holds exactly that at 7:140753336 A>T    PASS
SpliceAI DS_* = 0.00       the precomputed table agrees                            PASS
```

and then the two surprises.

### Surprise 1 — the most famous hotspot in oncology is not flagged as a hotspot

```
CancerHotspots_HOTSPOT         (empty)
CancerHotspots_HOTSPOT3D       1 / BRAF / p.V600E
```

That is **correct**, and the verifier passes it. The reason is allele representation. In the
Cancer Hotspots file the linear hotspot record for V600E is a two-base MNV one position to the
left, while the 3D record is the single-base substitution our caller emits:

```
chr7  140753335  CA  TT   HOTSPOT;HOTSPOT_GENE=BRAF;HOTSPOT_HGVSp=p.V600E     <- not our spelling
chr7  140753336  A   T    HOTSPOT3D;HOTSPOT3D_GENE=BRAF;HOTSPOT3D_HGVSp=p.V600E
```

So an exact-match annotation finds the 3D record and misses the linear one. Anyone filtering on
`CancerHotspots_HOTSPOT` alone would not see BRAF V600E. The classifier reads neither column, so
no classification is affected — but it is exactly the sort of thing that is invisible until
something compares the output with the source.

### Surprise 2 — REVEL was silently absent, and it changed a classification

```
REVEL = (empty)            new_tabbed_revel_grch38.tsv.gz holds 0.931 for this variant
```

A genuine defect, and this is how it was found. REVEL's file records one transcript per score
using older Ensembl identifiers: `0.931` is filed under `ENST00000288602`, while the pipeline
reports `ENST00000646891`. The plugin matched by transcript, found nothing, and wrote an empty
cell — indistinguishable from "no score exists".

It mattered because `rule_o5` reads that column and returns `None` when it is empty, with no
fallback. So O5, worth +4 or +2, could not be awarded for any BRAF missense variant. 9 of the 71
missense rows of a 93-variant cohort were affected.

Fixed by running the plugin with `no_match=1`, since REVEL scores the amino-acid substitution
rather than the transcript structure and its file states the substitution in `aaref`/`aaalt`. The
verifier now checks that proviso instead of assuming it: a cross-transcript score whose
substitution agrees passes **with a note naming the transcript it came from**, and one whose
substitution differs does not.

## What the tests do not cover

Worth stating plainly, because it is where the gaps have been.

Of the 427 described columns, **47** are verified against the sources' own APIs and **146** against
the files on disk -- 193 distinct columns in all. The rest come from the input VCF, are computed by
SMART, or have no external source to compare against.

An earlier version of this section said "219 ... against an API". 219 was a count of source
attributions in Config.yaml, not of anything verified, and it overstated API coverage by 4.7x.
The two figures above are counted from `verification1/results.tsv` and
`verification1/local_source_results.tsv`, and a test now checks them against those files.
Within the 146:

- **GENIE** is not compared but **recomputed** — its columns are aggregates. The O4 match columns
  are deliberately left as invariants only: reproducing the matching rules would be
  reimplementing the classifier, and a verifier running the pipeline's own logic would agree with
  it by construction and prove nothing.
- **Indels** are skipped by several modules, because matching them against a source needs the
  same left-alignment the pipeline applies.
- **dbNSFP** returns multi-transcript values with no indication of which is ours, and its
  positional conservation scores are unreachable for non-missense variants. Both are
  [issue #20](https://github.com/WeTGI-colab/SMART/issues/20).

Coverage of the Python itself is about 28% of statements, concentrated entirely in the pure
logic; the uncovered lines are the pipeline orchestrators and HTTP calls, which the
container-based suites exercise instead. `verification5/README` has the measurement.

## Reading the results files

**Never read a MAF or a results TSV with `comment="#"`.** pandas honours a comment character
anywhere in a line, not only at its start, and real annotation values contain `#` — ClinVar's
`CLNVI` holds `UniProtKB:P01111#VAR_006847`. The row is truncated there and every later column
comes back empty. On verification1's output that silently discards **895 of 1,035 columns on 16
of the 22 rows**. Count the leading metadata lines and pass `skiprows` instead:

```python
skip = 0
with open(path) as fh:
    for line in fh:
        if line.startswith("#"):
            skip += 1
        else:
            break
df = pd.read_csv(path, sep="\t", skiprows=skip)
```

This is [issue #4](https://github.com/WeTGI-colab/SMART/issues/4). It is the mistake that
produced 350 phantom mismatches on the first run of the source verifier.

## Columns that never carry a value

`ONCOKB_VUS` was empty on every row of every run for years. `oncokb2.0.py` read the OncoKB
field as `VUS`; its name is `vus`. Nothing caught it, and nothing could have, because three
different situations are the same empty cell by the time they reach the MAF:

- the source has no value for this variant,
- the column does not apply to this kind of variant,
- the column is broken and will never carry a value.

`dict.get` returns `None` for an absent key exactly as it does for a null value, and the
pipeline writes `.` for both. The bug was found only because a variant was added that has
`vus=True` — with the previous cohort the column would have looked legitimately empty forever.

`scripts/audit_empty_columns.py` turns that into something a test can hold. It asks of every
column: *did any variant, in any run, give it a value?* Those that never did must be listed in
`verification1/empty_columns_allowed.tsv` with a reason and one of two kinds:

| kind | meaning |
|---|---|
| `structural` | no variant could fill it as the pipeline runs. Changing this means changing the pipeline — e.g. everything needing a tumour-type query ([#26](https://github.com/WeTGI-colab/SMART/issues/26)) |
| `cohort` | a variant would fill it and we have none. The reason names the kind of variant needed |

On the current output: 1,148 columns, 106 never carry a value — 76 structural, 30 `cohort`.
The `cohort` half is printed on every run as a to-do list, and it is the answer to "what should
we add next": a non-coding hotspot, a variant in a miRNA, one inside a regulatory motif, an
upstream variant for `DISTANCE`, a ClinVar "included" variant, one with Dx1 evidence, a
transcript VEP flags.

The audit fails in **both** directions, and the second is the one that matters:

- a column that is empty and unlisted — possibly the next `ONCOKB_VUS`;
- a listed column that now carries data — the entry is stale and must be removed.

So adding a variant forces the blank it fills to be re-examined, instead of letting it
accumulate unexamined. `tests/verification5/test_audit_empty_columns.py` tests the logic on
synthetic tables, needing no Docker, refs or token, and runs the audit over the real output
when it is present.

A related guard covers the cause rather than the symptom:
`tests/verification5/test_oncokb_response_keys.py` tokenises `oncokb2.0.py`, extracts every
literal key it reads off an OncoKB response, and requires each to exist in a frozen real
response. All 26 other keys were audited and were correct; `VUS` was the only one wrong.

## Fixture defects, and what now prevents them

`tests/verification1` had two, both fixed and both recorded in
[issue #19](https://github.com/WeTGI-colab/SMART/issues/19) and in `verification1/README.md`.
They are worth reading before comparing against any figure published from this fixture, because
the figures changed.

Its two committed copies held different variants — the `.vcf.gz` the pipeline reads had 22, the
plain `.vcf` had 18 — and three records carried a REF allele absent from GRCh38, two of them
minus-strand genes whose alleles had been written on the gene's strand. VEP does not reject such
a record, so the suite had been verifying variants nobody submitted.

Two guards now exist, and each addresses one of them. `scripts/check_reference.py` runs before
VEP and by default refuses a run whose REF alleles disagree with the reference, so that class of
defect cannot reach an output again. And `tests/verification5/test_check_reference.py` asserts
that the fixture's two copies stay identical and that its corrected alleles stay corrected,
because regenerating one copy from the other is precisely how four variants stopped being
tested without anything reporting it.
