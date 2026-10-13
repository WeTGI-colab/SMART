# Changelog

All notable changes to SMART are recorded here. The version in [`VERSION`](VERSION) is
bundled in the Docker image, printed at the start of every run and written as the first
line of every output MAF (`#SMART_VERSION x.y.z`), so any result file can be traced back
to the entry below.

Docker images are published as `monkiky/smart:<version>`.

## [2.0.0] — unreleased (in development)

SMART stays an annotation pipeline (Tiers 1–3); classification is a separate output.

### Added
- **Every annotation source is now verified against its own data.** `tests/verify_local_sources.py`
  gained modules for dbNSFP, GENIE and gnomAD, closing the set: on the 93-variant public cohort,
  11,253 checks over ClinVar, Cancer Hotspots, SpliceAI, REVEL, LOEUF, dbNSFP and GENIE give
  10,597 PASS, 600 SKIP and 56 MISMATCH, and gnomAD adds 416 checks with none. GENIE is verified
  by recomputing its same-nucleotide counts from the counts table rather than by comparing a
  value, since its columns are aggregates; the O4 match columns are checked for invariants only,
  because reproducing the matching rules would be reimplementing the classifier and would agree
  with it by construction. Normalisation had to absorb four list separators -- VEP rewrites `,`,
  `|` and dbNSFP's `;` all to `&` -- a relative tolerance, since the VEP cache writes three
  significant figures where the REST API writes four, and the difference between a flag's
  presence written as 1 and a genuine score of 1.0, which is decided where both sides are in view
  rather than in the normalisation.
- **Acknowledgement of AACR Project GENIE in the manuscript**, as its terms of access require
  alongside citation of the consortium publication and a statement of the dataset version.

- **gnomAD verified against the VEP API; 26 more columns covered.** gnomAD has no file in the
  pipeline -- `gnomADe_AF`, `gnomADg_AF` and the population columns come from the VEP cache via
  `--everything` -- so it was left out of the local-source verifier for want of anything to look
  up, and out of `verification1` on the grounds that the REST API's gnomAD release might differ
  from the cache's and so could not prove a mismatch wrong. The API returns them and they agree:
  416 checks across 16 variants, no disagreement, `gnomade = 4.857e-05` against our
  `gnomADe_AF = 4.857e-05` and so on for every population. The exclusion was unnecessary.
  Comparison needs a relative tolerance because the cache writes three significant figures and
  the API four, and the retry policy treats a 4xx other than 429 as a permanent refusal -- a 400
  here means the variant itself is malformed, which is how the strand-flipped allele of #19
  presents.

- **The REF alleles of a VCF are verified against the reference genome before annotation**
  (`scripts/check_reference.py`, `--ref-check fail|warn`). Nothing checked this, and no other
  step catches it: VEP does not reject a VCF whose REF disagrees with the FASTA, it takes the
  reference base from the genome and keeps the submitted ALT, so the consequence it reports
  belongs to a variant nobody submitted. Found in SMART's own test data —
  `tests/verification1/verification1.vcf` holds three records whose REF does not match
  GRCh38. Two are minus-strand genes whose alleles were written on the gene's strand:
  `GNAQ_Q209L` (chr9:77794572, VCF `A`, genome `T`) was reported as `p.Gln209=` `Silent`,
  because the flipped ALT equalled the true reference and the variant annotated itself out of
  existence, and `DICER1_E1705K` (chr14:95094139, VCF `G`, genome `C`) was reported as
  `p.E1705*` `Nonsense_Mutation`, a different substitution from the one its name states. The
  third is an off-by-one REF on an EGFR deletion. None produced a warning; the only step that
  noticed was the SpliceAI lookup, and only by omission, since it matches on CHROM+POS+REF+ALT
  and those records silently received no splice annotation. The check reads the reference
  through its `.fai` index, skips records that state no literal sequence (symbolic ALTs,
  breakends, placeholder REFs), distinguishes a wrong assembly from a handful of strand flips,
  and runs on the file VEP is about to read — after liftover, so a GRCh37 input does not fail
  by design. The default stops the run, since one wrong allele yields a wrong consequence and
  therefore a wrong classification. Documented in `docs/reference_check.md`.

- **Annotation flags, so a filter cannot drop a variant before its evidence is read**
  ([#14](https://github.com/WeTGI-colab/SMART/issues/14),
  [#15](https://github.com/WeTGI-colab/SMART/issues/15)). SpliceAI already reached the output
  (`SpliceAI_pred_DS_*`, columns 113–116 of 807) and the two TERT promoter hotspots were already
  annotated correctly as `5'Flank`, but a filter on `Variant_Classification` or on `IMPACT`
  removes those rows as `Silent`, `Intron`, `5'Flank`, `LOW` or `MODIFIER` before anyone reaches
  the evidence: the data was there and the row was gone. Across three cohorts (10,124 rows) that
  is 11 variants with a SpliceAI delta of 0.5 or above and about 25 at 0.2, several of them
  intronic changes three to five bases into the splice region with deltas of 0.90–0.99. Two new
  columns: `SPLICEAI_MAX_DS`, the highest of the four deltas as one sortable number, computed for
  every row; and `ANNOTATION_FLAGS`, `;`-separated reasons to keep a row — `splice_predicted(<score>)`
  for a `LOW`/`MODIFIER` row reaching 0.2, the delta from which the classifier awards O6, so every
  variant that could earn O6 now reaches the classifier; and `svig_canonical` for a position on the
  SVIG-UK canonical variants list, read from the same reference file the classifier uses for O1 and
  matched genomically, since a promoter variant's `c.` notation depends on the transcript annotated
  while its position does not. The consuming filter gains one condition and one column carries every
  reason. `Variant_Classification` is deliberately unchanged: a splice-predicted synonymous variant
  is still synonymous, a value outside the MAF vocabulary breaks third-party readers, and relabelling
  destroys the record of the consequence. Verified on a real per-sample MAF (22 rows, 1,035
  columns) that both columns are added, parse the numeric cells VEP actually writes and reach
  the Tier 2 and Tier 3 sheets. Documented in `docs/annotation_flags.md`.

- **Field provenance taken from the VCFs** (`scripts/vcf_metadata.py`, [#16](https://github.com/WeTGI-colab/SMART/issues/16)). The description, source
  and version of every column used to come from `Config.yaml`, maintained by hand, which had
  drifted: `hotspot` was attributed to Cancer Hotspots when the variant caller writes it,
  `FractionInformativeReads` to Manta when DRAGEN writes it, six further DRAGEN fields had no
  entry, and 475 of 807 Tier 1 columns were undocumented. A VCF documents its own fields, for the
  run that produced them, so the headers are now harvested — the caller's own `##INFO`/`##FORMAT`
  entries and `##source`, VEP's version and cache, each custom annotation's resource file, each
  plugin's version, and OncoKB's fields. The VCF supplies description, source and version for
  what it covers; `Config.yaml` keeps the tier, the MAF standard columns, what SMART computes and
  the OncoKB MafAnnotator columns, whose version is the OncoKB release observed in the run.
  Checks report fields lost between the caller and the table, an input `INFO` name colliding with
  an annotation field (an error: the caller's value would silently replace it), a run straddling
  two OncoKB releases, and samples within a run disagreeing on a resource. Documented in
  `docs/field_provenance.md`; without `--vcf-metadata` the pipeline behaves as before.

- **Gene-specific splice provisions (classifier v2.3).** SVIG-UK §7.11 states that synonymous
  changes at the last base of exons 4, 6 and 9 of TP53 "should be classed as driver variants as
  they affect splicing of these exons". Those positions are not part of the canonical +/-1,2
  dinucleotide and in silico predictors can miss them, so such a variant was previously assessed
  as a synonymous change. A matching variant now enters the O2 exon-skipping branch, so its
  strength comes from the PVS1 tree (reading frame, NMD, fraction of protein lost) rather than a
  fixed weight, and a review flag records that the guideline's "driver" wording is not itself a
  code or a strength. New reference list `reference_lists/svig_uk_gene_specific_splice.tsv`
  (`--gene-specific-splice`), with the positions verified against the UCSC ncbiRefSeq exon
  structure.
- **B4 exon-edge caution.** B4 "should be used with caution for synonymous variants ... located at
  the first base or last three bases of the exon" [Suppl. p.16], so it is no longer applied there;
  B3 is unaffected and a review flag is raised. Applies to every gene, not only TP53.
- **Functional-assay evidence for O10/B6 (classifier v2.3).** O10 was derived only from a proxy
  (OncoKB mutation effect, capped at supporting), while the guideline permits up to very strong
  and gives the operative criterion that a variant is functionally abnormal when protein activity
  is <25% of the wild-type level. Assay measurements now drive O10 where they cover the variant,
  at the strength recorded per assay, with the OncoKB proxy as the fallback and a review flag
  noting which was used. New reference lists `reference_lists/svig_uk_functional_assays.tsv` and
  `svig_uk_functional_assay_scores.tsv` (`--functional-assays`), built by
  `utils/build_functional_assays.py` from MaveDB. The first assay is the median TP53
  transcriptional activity across eight promoter-specific assays (Fortuno et al. 2021,
  PMID 34273903, the ClinGen TP53 VCEP model; 2,314 variants, CC0). The boundaries the guideline
  does not fix — "intermediate" and the B6 side — are left unset, so neither is applied until a
  laboratory records them, and the assay strength ships at the conservative end of the permitted
  range for the same reason.

- **Tier 4 — SVIG-UK oncogenicity classification** (`Final_result_tier4.tsv`), written after
  post-analysis by `scripts/svig_uk_classifier.py`. Two-row header documenting the
  classification system (description | source | version), one column per SVIG-UK code with
  the strength and points applied, the evidence for each code and the codes removed by
  combination rules. New pipeline options `--tumour-type`, `--o4-on-target` and
  `--no-classification`.
- SVIG-UK classifier following ACGS/SVIG-UK guidelines v1.1 (26/05/2026), O1 canonical
  variants list and COSMIC resource builder.
- **Tier 4 review support.** `SVIG_UK_REVIEW_STATUS` states that the classification is an
  automated suggestion to be confirmed by a clinical scientist, and `SVIG_UK_REVIEW_FLAGS`
  lists reasons to look closely (missing resources, minimum-two-codes downgrade, conflicting
  evidence, score within 1 point of a class boundary, unknown gene role, CHIP correction, codes
  not automated). New structured export `Final_result_tier4.json` (`--tier4-json`) for review tools.
- **Gene-specific population thresholds for hereditary cancer genes** (SVIG-UK O3 and B1 require
  ClinGen VCEP cut-offs where they exist): PM2 -> O3, BA1 -> B1 stand-alone, BS1 -> B1 strong for
  TP53, RUNX1, PTEN, CDH1, ATM, PALB2, BRCA1, BRCA2, APC, MLH1, MSH2, MSH6, PMS2, VHL and DICER1
  (`reference_lists/svig_uk_vcep_population_thresholds.tsv`, from the ClinGen CSpec registry).
  GrpMax filtering allele frequencies are estimated from gnomAD allele frequencies.
- **GENIE case counts as annotation** (Tiers 1–3), independent of classification:
  `GENIE_SAME_NUCLEOTIDE_*`, `GENIE_O4_MATCH_*` (counted as for O4), split into all /
  haematological / solid tumours, with the most frequent OncoTree codes, patients seen in more
  than one tumour type and the cohort size (`scripts/genie_annotation.py`, enabled when
  `GENIE/genie_o4_counts.tsv.gz` is present; the builder now also writes cohort sizes).

### Fixed
- **O5 could not be awarded for any gene whose reported transcript is absent from REVEL's
  file, BRAF among them.** `rule_o5` reads the table's `REVEL` column and returns None when it
  is empty, with no fallback -- the classifier's own REVEL lookup is used only for the
  *reference* variant's score, not the assessed one. Because the plugin matched by transcript
  and REVEL's file records older Ensembl identifiers, the column was empty for 9 of the 71
  missense rows (13%) of the 93-variant cohort: all 8 BRAF variants on `ENST00000646891`, and
  one FLT3. So O5, worth +4 moderate or +2 supporting, was silently unavailable for those
  variants. Fixed by running the plugin with `no_match=1` (see the REVEL entry above);
  classifications computed before that fix did not have O5 available for the affected genes.

- **REVEL silently returned nothing for variants annotated on a newer transcript, BRAF V600E
  among them.** The plugin matches by transcript by default, and REVEL's source file records one
  transcript per score using older Ensembl identifiers: it holds 0.931 for BRAF V600E on
  `ENST00000288602` while the pipeline reports `ENST00000646891`, so the lookup found nothing and
  the most frequently reported variant in the panel had no REVEL score, with nothing in the
  output to say whether that meant "no score" or "not assessed". Found by the local-source
  verifier. REVEL scores the amino-acid substitution rather than the transcript structure, so a
  score computed on another transcript is the right answer provided the substitution is the same,
  which its file states in its `aaref`/`aaalt` columns; the plugin now runs with `no_match=1` and
  the column description says what the score is and is not. The verifier tests the proviso rather
  than assuming it: a cross-transcript score whose substitution agrees passes with a note naming
  the transcript it came from, and one whose substitution differs does not. Of 13 missense rows in
  `verification1` only BRAF was affected, and it was V600E.

- **`VAF` was described as a population frequency when it is the tumour allele fraction.**
  `Config.yaml` described `VAF` as "Frequency of existing variant in 1000 Genomes" with
  `source: VCF`, but SMART computes it from the caller's allelic depth as ALT/(REF+ALT). Both
  `AF` (VEP's population frequency) and `VAF` are tier 3, so the sheet carried two columns
  described as population frequencies, one of which was the patient's allele fraction. The
  description and source are corrected; the computed values never changed. Found while adding
  the annotation flags, and the same class of drift as
  [#16](https://github.com/WeTGI-colab/SMART/issues/16).
- **Duplicate `VAF` key in `create_config/ConfigORIGINAL.yaml`.** The field was defined twice;
  in YAML the later definition silently wins, so the first was dead. Consolidated into one.

- **OncoKB did not recognise known variants when the selected transcript is not OncoKB's.**
  OncoKB curates protein changes on one RefSeq transcript per gene; with another transcript
  the residue numbering can differ (NF1 NM_001042492 p.R1769* is NM_000267 p.R1748*;
  CUX1 NM_001913 vs NM_181552), so those variants came back as Unknown (1.0.8 known issue).
  Variants whose selected transcript differs from OncoKB's curated transcript are now
  queried by genomic change (GRCh38), both in the per-transcript `ONCOKB_*` annotation
  (`ONCOKB_QUERY_TYPE` = `MUTATION_BY_GENOMIC_CHANGE`) and in MafAnnotator (new column
  `ONCOKB_MAF_QUERY`). If OncoKB's genomic-change service fails, SMART falls back to the
  protein-change query and records it; no variant is dropped. The transcript list reported
  in SMART's own columns is unchanged.

- GENIE O4 undercounted truncating and in-frame entries next to repeats: GENIE indels are
  left-aligned, which on the plus strand placed repeat-shifted events one residue upstream
  (NPM1 p.W288Cfs*12 counted 89 instead of 839 patients at residue >=288). Residues are now taken
  after shifting indels 3' in transcript direction, as HGVS does.

### Changed
- Canonical variants list and combination-rule citations updated to guideline v1.1:
  RUNX1 c.601G>A corrected to c.602G>A, H3-3B c.111G>C removed, NPM1 types A, B and D
  listed; combinations renumbered (#10 O9/O2, #12 O9/O6, #9 O8/O7, #15 O10/O2) and
  O5/O7 (#5), O5/O8 (#8), O5/O10 (#16) cited.
- Optional dbNSFP 5.0a annotation (VEP dbNSFP plugin, enabled when
  `dbNSFP/dbNSFP5.0a_grch38.smart.gz` is present): 46 columns from 30 predictors and
  conservation scores that SMART did not already annotate, one version per tool, no rank
  scores; SIFT, PolyPhen, REVEL and CADD stay from their current sources. New builder
  `utils/build_dbnsfp.py`; column descriptions, sources and tiers in `Config.yaml`
  (AlphaMissense in Tier 3, the rest in Tier 2). dbNSFP is licence-restricted and not shipped.
- O2 follows the full PVS1 decision tree using RefSeq exon structure (UCSC `ncbiRefSeq`)
  and the GRCh38 FASTA: NMD predicted from the stop position (final exon / last 50 bp of
  the penultimate exon escape), >10% protein rule for NMD-escaping variants, exon-skipping
  frame and NMD for canonical splice variants, cryptic-site and first/final-intron
  handling. Frameshifts without a VEP protein change are translated from the genome
  (checked against all 177 VEP-annotated frameshifts in the validation set: 177 identical).
  On the WGLS set (969 variants) sensitivity rose from 0.664 to 0.724 with no new false
  positives (in development; not a frozen evaluation). With the O3 CHIP correction, GENIE
  O4 and O5 added below, sensitivity is 0.789 with specificity 0.996 (2 false positives,
  unchanged).
- O3 CHIP correction for DNMT3A, ASXL1, TET2 and TP53 as described in the guideline:
  gnomAD v4.1 frequency recalculated from carriers with VAF 0.30–0.65 when neither O3 nor
  B1 applies and total AF >0.001%. New builder `utils/build_gnomad_chip_ab.py` extracts the
  allele-balance histograms from the public gnomAD VCFs.
- O4 from AACR Project GENIE with the guideline's thresholds: unique patients, on-target
  lineage (OncoTree), GRCh37 entries lifted to GRCh38 and matched by genomic change so that
  counts do not depend on GENIE's transcript numbering. New builder `utils/build_genie_o4.py`
  (output is GENIE-derived and must not be redistributed). COSMIC remains the fallback.
  The tumour type can be an OncoTree code, and `--o4-on-target code|tissue|lineage|all`
  chooses which GENIE entries are on-target (default `lineage`).
- O5 missense route: a different (likely) oncogenic amino acid change at the same residue
  (SVIG-UK canonical list, optional CGI catalogue references) compared by REVEL, with the
  O5/O7, O5/O8 and O5/O10 combination rules. New builder `utils/build_o5_references.py`.

### Fixed
- Splice O2 was applied to repeat-shifted indels that VEP flags as splice donor/acceptor
  although the 3'-shifted HGVS change lies outside the ±1/2 dinucleotide
  (e.g. a 4-bp deletion at +5..+8 inside a repeat that starts at the donor GT).

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
