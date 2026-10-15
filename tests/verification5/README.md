# Verification 5 — unit tests

Pure-Python unit tests over the pipeline's logic. **No Docker, no OncoKB token and no reference
files are needed** — the suite runs in any standard Python environment in about a second, which
is what makes it the layer to run on every change.

```bash
python3 -m pytest tests/ -q                      # from the repository root
bash tests/verification5/run_verification5.sh    # or through the runner
bash tests/verification5/run_verification5.sh --coverage
```

See [`tests/README.md`](../README.md) for how this fits with the other layers.

## What is tested

| Script | Functions tested |
|--------|-----------------|
| `vcf2table.py` | `calculate_end` (all variant classes), `convert_hgvsp_short` (3→1 letter AA conversion), `load_nm_transcripts` (version stripping, edge cases) |
| `oncokb2.0.py` | `get_cancer_from_filename`, `clean_value`, `stringify_list`, `_extract_one_letter_protein`, `classify_variant` (all SVTYPE + MantaID combinations), `map_svtype_to_cna`, `load_preferred_transcripts`, `get_csq_index`, `parse_csq_format` |
| `post_analysis.py` | `load_field_config`, `resolve_field_meta` (exact + wildcard pattern), `standardize_format_column`, `normalize_format_fields`, `add_vaf_column`, `clean_column_values`, `calculate_end_position`, `drop_columns`, `add_annotation_flags` |
| `svig_uk_classifier.py` | every evidence code against the guideline text, the scoring thresholds and the combination restrictions (211 tests) |
| `vcf_metadata.py` | harvesting field provenance from the run's own VCF headers, precedence against `Config.yaml`, and the five consistency checks |
| `check_reference.py` | reading a FASTA through its `.fai`, which records are checkable, and detection of strand flips, off-by-one REFs and a wrong assembly |
| `verify_local_sources.py` | the comparison logic: contig naming, list separators, rounding tolerance, flags against numbers, and the GENIE recomputation |
| `verify_reference_files.py` | the reference-integrity check: the fused-header defect found in the lab's CIViC VCF, the bootstrapping hazard (a manifest generated from a corrupt directory must not bless it), truncation and error-page stubs, plain gzip where bgzip is required, stale and mismatched indexes, `.fai` consistency, required versus optional sources, and the exit statuses `entrypoint.sh` switches on |

Per file: `test_svig_uk_classifier` 211, `test_oncokb` 70, `test_annotation_flags` 63, `test_verify_local_sources` 60, `test_post_analysis` 48, `test_vcf2table` 43, `test_check_reference` 42, `test_vcf_metadata` 28, `test_vep_only_mode` 27, `test_multi_transcript` 20, `test_oncokb_maf_query` 8, `test_build_dbnsfp` 7, `test_genie_annotation` 5.

**1,183 tests, all passing.** Real variants from the verification1 dataset (NRAS Q61R, IDH1 R132H, PIK3CA E545K/H1047R, BRAF V600E, TP53 R175H/R248W, EGFR exon 19 deletion) are used as parametrised ground-truth cases.

## Code coverage

Coverage is measured over the pure-Python lines of each script. Lines excluded from coverage are those that require Docker, a live OncoKB API connection, or a real VCF file to execute (the main pipeline loops — `annotate_vcf`, `merge_maf_files`, `main` — and the HTTP query functions).

| Script | Statements | Covered | **Coverage** |
|--------|-----------|---------|--------------|
| `vcf2table.py` | 216 | 40 | **19%** |
| `oncokb2.0.py` | 401 | 122 | **30%** |
| `post_analysis.py` | 342 | 108 | **32%** |
| **Total** | **959** | **270** | **28%** |

The uncovered lines are concentrated in three pipeline orchestrators (`main`, `annotate_vcf`, `merge_maf_files`) and the OncoKB HTTP query functions. These are exercised end-to-end by verifications 1–4. The unit tests cover **100% of the pure, stateless logic** that can be tested without external infrastructure.

To reproduce the coverage report:

```bash
cd tests/verification5
python3 -m pytest . --cov=../../scripts --cov-report=term-missing
# HTML report (open coverage_html/index.html in a browser):
python3 -m pytest . --cov=../../scripts --cov-report=html:../../coverage_html
```

---

---

The 49 dbNSFP columns, the GENIE aggregates and every other annotated value are checked by a
different layer — against the files they were annotated from, not against this suite's
expectations. See [`tests/verify_local_sources.py`](../verify_local_sources.py) and
[`tests/README.md`](../README.md).
