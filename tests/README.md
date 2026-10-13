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
| Unit tests | `verification5/` | python3 + pytest | the logic, 632 tests in about a second |
| Caller compatibility | `verification2/` | Docker, refs | four callers' VCF conventions (manual) |
| Transcript choice | `verification3/` | Docker, refs, token | the whitelist changes the annotation and the OncoKB call |
| Parallel jobs | `verification4/` | Docker, refs, token | `--jobs N` matches a sequential run |
| API verification | `verification1/` | Docker, refs, token | 219 columns against the sources' own APIs |
| Source verification | `verify_local_sources.py` | refs | the other 151, against the files on disk |
| Variant presence | `verify_maf.py` | — | expected variants are in an output |

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

Of the 425 described columns, 219 are verified against an API, 151 against the files on disk, and
55 come from the input VCF or are computed by SMART and have no external source to compare
against. Within the 151:

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

## Known defects in the fixtures

`tests/verification1/verification1.vcf.gz` holds 22 variants and is what the pipeline reads;
`verification1.vcf` holds 18, is two months older, and reading it to see what is tested gives
the wrong answer. Three of its records carry a REF allele that does not match GRCh38, two of
them minus-strand genes whose alleles were written on the gene's strand. Both problems, and the
correct alleles, are in [issue #19](https://github.com/WeTGI-colab/SMART/issues/19) and in
`verification1/README.md`.

`scripts/check_reference.py` now refuses a run whose REF alleles disagree with the reference, so
this class of defect cannot reach an output again.
