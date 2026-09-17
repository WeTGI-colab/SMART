# Changelog

All notable changes to SMART are recorded here. The version in [`VERSION`](VERSION) is
bundled in the Docker image, printed at the start of every run and written as the first
line of every output MAF (`#SMART_VERSION x.y.z`), so any result file can be traced back
to the entry below.

Docker images are published as `monkiky/smart:<version>`.

## [1.0.8] — 2026-09-17

### Fixed
- **Malformed `HGVSp_Short` for frameshifts and multi-residue changes.** The short
  protein notation only handled substitutions: frameshifts came out as, e.g.,
  `p.T15GlufsTer` (mixed one- and three-letter codes, stop position lost) instead of
  `p.T15Efs*55`, and in-frame deletions, delins, duplications and synonymous changes
  stayed in three-letter form. `HGVSp_Short` is the column MafAnnotator uses to query
  OncoKB. On the internal validation set the fix changed `HGVSp_Short` for about a third
  of rows; OncoKB tolerated most of the old frameshift strings, but exact known variants
  such as NPM1 `p.W288Cfs*12` were not recognised as such (`VARIANT_IN_ONCOKB` False → True).
  **Re-running samples processed with earlier versions is recommended when frameshift or
  in-frame variants matter.**
- **A single dropped connection to the OncoKB API aborted the whole sample.** API calls
  are now retried with exponential backoff (5 attempts) on dropped connections, timeouts,
  rate limiting and gateway errors. If the API stays unreachable the sample fails with a
  clear error instead of producing variants without OncoKB annotation.
- README: the command to read the version bundled in an image pointed at the wrong path.

### Added
- This changelog.

### OncoKB data updates
OncoKB content changes independently of SMART releases, so the same SMART version can
annotate differently over time. Every output row records `ONCOKB_DATA_VERSION`; compare it
before comparing results between runs.

- **v7.2 → v7.5** (observed between internal validation runs in 2026). With
  `HGVSp_Short` unchanged, 26 oncogenicity calls and 11 highest-sensitive-level values
  changed on the validation set, mostly splice-region and intronic variants moving from
  `Unknown` to `Likely Oncogenic`.

### Known issues
- **OncoKB protein-change queries assume OncoKB's transcript numbering.** When the
  transcript chosen for a gene numbers residues differently from OncoKB's canonical
  transcript, the query can return `Unknown` and drop treatment levels. Example: NF1 on
  NM_001042492 carries 21 more residues than OncoKB's NM_000267, so a frameshift at
  p.1406 on the former is p.1385 on the latter; querying with the first numbering returns
  `Unknown`, the second returns `Likely Oncogenic` with LEVEL_1. CUX1 on NM_001913 versus
  OncoKB's NM_181552 is affected the same way. Planned fix: query OncoKB by genomic
  change, which does not depend on the transcript.

## [1.0.7] — 2026-08-05

### Fixed
- Embedded newlines and tabs inside OncoKB text fields could split a row and break the
  MAF. Cell values are now sanitised before writing.

## [1.0.6] — 2026-06-30

### Changed
- Post-analysis is memory-scalable: output tables are written with a two-pass streaming
  approach instead of loading every sample into memory.

## [1.0.5] — 2026-06-19

### Added
- Optional CADD deleteriousness annotation, enabled when CADD files are present in the
  reference directory.

## [1.0.4] — 2026-06-19

Docker image only; there is no matching version commit in git. Built from the 18 June
follow-up changes to 1.0.3: COSMIC columns surfaced in the Tier 2 output and field
descriptions in `Config.yaml` grounded in their source documentation.

## [1.0.3] — 2026-06-18

### Added
- Optional COSMIC Cancer Mutation Census annotation, enabled when a COSMIC VCF built with
  `utils/cosmic_cmc_to_vcf.py` is present in the reference directory. COSMIC is
  licence-gated and never shipped with SMART.

## [1.0.2] — 2026-05-30

### Changed
- `--output-dir` is mandatory: results must be written to a mounted volume.

### Fixed
- Post-analysis working directory when `--output-dir` differs from the container's
  working directory.

## [1.0.1] — 2026-05-29

### Fixed
- Documentation of Tier 1 column counts and of example variant counts.

## [1.0.0] — 2026-05-16

First release, matching the Zenodo release tag.
