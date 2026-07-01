# Reference lists

Ready-to-use gene / transcript lists distributed with SMART.

## OncoKB canonical transcripts
| File | Description |
|------|-------------|
| `oncokb_transcripts.txt` | OncoKB's canonical RefSeq **NM** transcript per curated gene, one per line. Pass directly to the pipeline as `--transcripts-file` (Tier 1 of the transcript whitelist). |
| `oncokb_transcripts_summary.tsv` | The same list with extra context — gene symbol, NM accession, Ensembl transcript — for review. |

Both are produced by [`utils/get_oncokb_transcripts.py`](../utils/get_oncokb_transcripts.py); re-run it to regenerate against the current OncoKB release.

## TSO500 panel lists
| File | Description |
|------|-------------|
| `TSO500_genes_list.txt` | Genes covered by the TruSight Oncology 500 panel. |
| `TSO500_transcripts_list.txt` | Preferred RefSeq NM transcript per TSO500 gene — a ready-made `--transcripts-file` whitelist. |
| `TSO500_transcript_MANE.txt` | TSO500 transcripts with their MANE Select mapping. |
| `TSO500_transcript_list_explained.md` | Notes explaining how the TSO500 transcript list was compiled. |
