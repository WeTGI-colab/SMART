# Verification 4 — Parallel job processing (`--jobs`)

## Purpose

Verifies that running SMART with `--jobs N` correctly processes multiple samples
in parallel on a single machine, producing the same output as a sequential run.

## What is `--jobs`?

The `--jobs N` flag tells SMART to process up to **N samples at the same time**
on the same machine (same CPU and RAM). It is process-level parallelism — not
distributed computing.

```
--jobs 1   →  sequential (default, one sample at a time)
--jobs 3   →  three samples run concurrently, sharing the machine's cores
```

This is useful when you have multiple VCF files and spare CPU capacity.
A typical VEP + OncoKB run takes 20–40 minutes per sample; running 3 in parallel
compresses wall-clock time to roughly the duration of the slowest sample.

### Resource guidance

Each parallel job runs a full VEP process (~4–8 GB RAM) alongside the OncoKB
API calls. Recommended limits by available RAM:

| RAM   | Safe `--jobs` value |
|-------|---------------------|
| 16 GB | `--jobs 2`          |
| 32 GB | `--jobs 3`          |
| 64 GB | `--jobs 6`          |

### Why not HPC nodes?

`--jobs` is intentionally designed for a **single machine with internet access**.
SMART requires internet connectivity for the OncoKB REST API, which is unavailable
on most HPC compute nodes. All parallelism therefore happens on the local machine
(e.g. your workstation or a cloud VM with internet access).

## How parallelism works internally

`entrypoint.sh` uses bash background jobs and `wait -n` (bash ≥ 4.3) as a
semaphore:

1. Each sample is processed by a `process_sample()` function run in a subshell
   (`&`), so all samples are fully isolated from each other.
2. When N jobs are already running, the main process blocks with `wait -n` until
   one finishes before launching the next.
3. Each sample writes its output to uniquely named files — no shared state during
   processing.
4. Per-sample stdout/stderr goes to `logs/<sample>.log` so output does not
   interleave in the terminal.
5. After all samples finish, per-sample variant count rows are merged into the
   shared `variant_counts.txt` in launch order.
6. `post_analysis.py` runs only after every sample has succeeded.

## What this test checks

`run_verification4.sh` copies `verification1.vcf.gz` twice under different sample
names (`par_sample_A`, `par_sample_B`) and runs the pipeline with `--jobs 2`.

`verify.py` then checks — without any external API calls:

| Check | What is verified |
|-------|-----------------|
| `log_exists:<sample>` | Per-sample log file was created in `logs/` |
| `log_complete:<sample>` | Log contains the `--- <sample> complete ---` marker |
| `counts_file_exists` | `variant_counts.txt` was produced |
| `counts_row_count` | Exactly one data row per sample (count-merge worked) |
| `counts_has_sample:<sample>` | Each sample name appears in the counts file |
| `maf_exists:<sample>` | Output MAF exists and is non-empty |
| `maf_identical:A_vs_B` | Both MAFs are content-identical after blanking `Tumor_Sample_Barcode` (same input VCF → same annotations) |

## Running

```bash
export ONCOKB_TOKEN="your-token"
bash tests/verification4/run_verification4.sh
```

Or via the full suite:

```bash
bash tests/run_all_verifications.sh --only verification4
```

## Files

| File | Description |
|------|-------------|
| `run_verification4.sh` | Sets up input VCFs and runs the Docker pipeline |
| `verify.py` | Checks parallelism plumbing and output correctness |
| `vcfs/` | Runtime-created input directory (not committed) |
| `output/` | Runtime-created output directory (not committed) |
| `results.tsv` | Latest verification results |
