# Changelog

All notable changes to SMART are recorded here. The version in [`VERSION`](VERSION) is
bundled in the Docker image, printed at the start of every run and written as the first
line of every output MAF (`#SMART_VERSION x.y.z`), so any result file can be traced back
to the entry below.

Docker images are published as `monkiky/smart:<version>`.

A section headed **unreleased** is work in progress: no git tag and no Docker image exist for
it yet, and it is written as the work lands rather than at release, so that anyone following
the repository between releases can see what is changing and why. An output file stamped with
an unreleased version came from a development run, and the entry below is where to find out
what that run did differently from the last published image.

## [2.1.0] — unreleased (in development since 2026-10-13)

SMART stays an annotation pipeline (Tiers 1–3); classification is a separate output.

### Added
- **SpliceAI is computed against the reported transcripts by default, or looked up on request**
  (`--spliceai compute|lookup`, `--spliceai-mask 0|1`; [#2](https://github.com/WeTGI-colab/SMART/issues/2),
  [#22](https://github.com/WeTGI-colab/SMART/issues/22)). Illumina's precomputed tables carry no
  transcript — one score per gene, computed on GENCODE v24 canonical of December 2015 — and were
  generated within a ±50 bp window, measured as exactly 50 over 199,971 records. Neither is
  configurable: both are properties of the files. For 16 of the 519 TSO500 transcripts the
  gene-level score therefore describes a different transcript from the one reported, MET, HRAS,
  AKT1, FGFR1 and CDKN2A among them; and where a transcript grew between GENCODE releases, part of
  it has no record at all — TET2's reported transcript has 37 kb, 28% of its length, uncovered,
  and 93 of one cohort's 168 TET2 variants sit there.

  Computing removes both, because the model takes a gene annotation as input. One annotation is
  built per run from the same `--transcripts-file` used for prioritisation
  (`scripts/build_spliceai_annotation.py`), rows named `GENE--REFSEQ`; since SpliceAI reports its
  annotation's NAME as the SYMBOL field, **the output records which transcript each score was
  computed on**, which the tables never allowed. The model runs on the file VEP annotates, after
  liftover; VEP runs without the plugin; and the scores are merged into the existing
  `SpliceAI_pred_*` columns, so the tier split, the `ANNOTATION_FLAGS` threshold, O6 and the
  cryptic-site check are unchanged whichever mode produced them.

  Validated against the Broad's own service, which is an independent implementation: for
  `16:2317763 T>A` a local run gives `0.01|0.21|0.04|0.00|1|106|128|497`, identical on all four
  deltas and all four positions, where our precomputed table gives 0.08 with every delta position
  within ±2 bp. Ten tests assert it against a frozen fixture of six public variants
  (`tests/verification5/test_spliceai_lookup_agreement.py`).

  Reference data for the two modes differs by four orders of magnitude: **91 GB** of precomputed
  tables against **10 MB** for `gencode.v43.annotation.txt.gz`, the model riding in the image
  (`tensorflow-cpu==2.15.*`, `spliceai==1.3.1`; TensorFlow pinned because the released package
  uses `tf.keras`, which Keras 3 removes). The cost moves to CPU, about 0.6 s per variant.

  Masked is the default, which is what the model's developers recommend for variant interpretation
  and what the laboratory uses; the tables we hold are the raw set, so this was previously raw by
  omission rather than by choice. Rows the model does not score are left empty rather than
  carrying a precomputed value, so a column never mixes the two sources silently, and
  multi-nucleotide variants — which the released model declines without saying so — are counted and
  warned about per sample.

  **Known limitation.** Supplying scores where there were none activates
  [#21](https://github.com/WeTGI-colab/SMART/issues/21): three TET2 variants in the development
  cohort move from Likely Oncogenic to VUS because the O2 cryptic-site check fires on an acceptor
  gain at +1 bp, which is the model re-detecting the junction the variant destroys rather than an
  alternative site. Masking does not prevent it. The rule is wrong independently of this change,
  but it is only reachable once the data exists, so #21 should be settled before cohorts are
  re-run on the new default.

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
  (`SpliceAI_pred_DS_*`, columns 1093–1096 of the 1,165-column Tier 1 MAF) and the two TERT
  promoter hotspots were already
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
  entry, and 475 of the Tier 1 columns were undocumented — of 807 at the time, 1,165 now. A VCF
  documents its own fields, for the
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
- **Two of the classifier's four reference lists were missing from the Docker image.** The
  `Dockerfile` copied `svig_uk_canonical_variants.tsv` and
  `svig_uk_vcep_population_thresholds.tsv` and not the rest. The classifier resolves each
  default path relative to `scripts/`, so inside a container the others were simply not found
  — and their absence does not fail a run:

  - O10/B6 fell back to the OncoKB mutation-effect proxy, capped at **supporting** — 1 point
    where the shipped assay table gives 2;
  - the TP53 exon 4/6/9 synonymous provision the guideline requires **never applied**.

  Neither emitted a note, so the same variant could be classified differently inside and
  outside a container with nothing in the output to explain it. These are the guideline's own
  content, already in the repository, under no licence restriction. All five are now copied,
  and a test asserts the Dockerfile copies every list the classifier names a default for —
  the general form, rather than the two that happened to be missing.

- **The reference verifier required the wrong VEP cache**, so it certified reference sets that
  cannot be annotated against and rejected ones that can. It required `homo_sapiens`; the
  pipeline runs `vep --merged`, which looks for `homo_sapiens_merged`. Both directions were
  wrong, and the worse one is that a directory holding only the plain cache passed
  verification and then failed inside VEP — which is precisely what the verifier exists to
  prevent. A machine with both caches sees neither failure, which is how it survived.

  `entrypoint.sh` now checks the merged cache **by name** in `REQUIRED_PATHS` rather than
  leaving VEP to discover it, so a missing cache fails before the run starts instead of deep
  in the annotation step.

- **Three documented tumour types could never be matched from a filename.**
  `get_cancer_from_filename` split the basename on `[_.-]` and required a single part to equal
  a `CANCER_MAP` key. Three of the seventeen keys contain hyphens, so after splitting on `-`
  no part could ever equal them: **HNSC, SNUC and THCA were unreachable**, and README names
  two of them as supported, spelled as the token a user would put in a filename. All three
  fell through to `UNKNOWN`, which makes OncoKB answer pan-cancer rather than refuse — so the
  failure was silent and the output looked complete.

  Spans are now rejoined and tried **longest first**, and that ordering is the substance of
  the fix: `thyroid` is THY and `thyroid-carcinoma` is THCA, two different OncoTree codes.
  Shortest-first would resolve `sample_thyroid-carcinoma.vcf` to THY and send OncoKB a tumour
  type the sample does not have, changing every tumour-specific level it returns.

  **A test asserted the defect.** `test_oncokb.py` expected `thyroid-carcinoma.vcf` to give
  `THY`, with the comment "split on `-` → `thyroid` matched first" — the mechanism of the bug,
  written down as the expectation. Someone saw the behaviour, understood why, and recorded it
  as correct. Corrected, with a note saying so; 17 further tests cover the rest, including one
  asserting that **every** key in the map is reachable, which is the general form of the bug.

- **`--tumour-column` was unreachable from the container.** `vcf2table` refuses an ambiguous
  tumour/normal VCF and names that flag in the refusal — and `entrypoint.sh` never parsed it,
  invoking vcf2table with only `--vcf --transcripts --output --debug`. So the message pointed
  at a remedy a Docker user had no way to apply, and this changelog offered it as the escape
  hatch. Both `--tumour-column` and `--normal-column` are now parsed, passed through with the
  `:+` form so that unset means absent rather than empty, documented in `--help`, and tested
  in both halves: that the entrypoint accepts the flag the refusal names.

- **CHANGELOG internal contradictions, found by auditing this file against itself.**

  There was **no `## [2.0.0]` section**, yet `smart:2.0.0` was named here as a published image
  and the image exists. The preamble promises that any result file can be traced back to an
  entry, so a MAF stamped `#SMART_VERSION 2.0.0` could not be. The section had been *renamed*
  to 2.1.0 in `476bfa8`, taking its content with it. A 2.0.0 section now says exactly that,
  and warns that 2.1.0's entry has since grown to describe more than that image contains.

  **Two `### Fixed` blocks under one version**, with `### Changed` between them — the same
  rename artefact. The second held one bullet; merged.

  **The two CIViC entries contradicted each other within this release.** One said the rebuilt
  file holds 781 records with none symbolic; the other still asserted that 480 of 1,261 carry
  a symbolic ALT, "however current the snapshot is", and called #24 open. The second was
  measured on the corrupt file and is now marked superseded in place rather than deleted,
  because the figure was published and someone may be holding it.

  Also corrected: `#24` described as open after it was closed; "columns 113–116 of 807" (they
  are 1093–1096 of 1,165); "475 of 807 Tier 1 columns" now dated, since 807 was true when
  written; and a "632 unit tests" figure inside the entry that was itself about fixing stale
  test counts, now dated rather than silently restated.

- **README numbers and links.** The tier column counts were a May snapshot: Tier 3 given as 77
  (78), Tier 2 as 670 (700), Tier 1 as 1,028 (1,165, stated in four places and inconsistently
  in a fifth). `ONCOKB_TX_0` … `ONCOKB_TX_33` where the indices run to 35. "7 ancestry groups"
  where gnomAD exomes carry 9 and genomes 10. `git clone` pointed at a personal account rather
  than the organisation that holds the repository, and the "Full field reference" link went to
  that account's GitHub Pages site. Added, because they were conspicuously missing: the
  **SVIG-UK score thresholds** (≥10 Oncogenic, 6–9 Likely Oncogenic, 0–5 VUS, −1 to −6 Likely
  Benign, ≤−7 Benign), quoted from `CLASS_THRESHOLDS`; and a note that **breakends get no
  fusion annotation**, since `MantaBND` goes to the mutation endpoint and
  `/annotate/structuralVariants` is never called — with the honest qualification that no
  SVIG-UK code depends on it and no published somatic SV framework exists to implement.

- **Eight of the nine `Config.yaml` files described the tumour allele fraction as a
  population frequency.** `VAF` read "Frequency of existing variant in 1000 Genomes", and the
  description is printed in the Tier 2/3 metadata row -- so a reporting scientist read a 0.33
  somatic allele fraction annotated as a 1000 Genomes frequency, in the sheet the README calls
  "optimised for clinical scientists". Only the root copy had the correction.

  That was one symptom of a larger failure: there are **nine** Config files, they are
  hand-maintained, and corrections had been landing in the root and nowhere else. Measured:
  **76 (column, file) pairs had drifted across 23 columns** -- LOEUF's version, citation and
  the gene-level claim corrected earlier in this release; the gnomAD AF versions; `hotspot`
  and `FractionInformativeReads` attributed to the wrong tool. Every one was a fix made once
  and propagated nowhere, including three the same release had just announced as fixed. The
  copies legitimately differ in *which* columns they list -- each cohort keeps the config it
  was run with -- but not in what a column means, and that distinction had been lost.

  Propagated by splicing the root's own text blocks, so the formatting matches by
  construction rather than being re-wrapped. **A test now asserts that every shared column is
  described identically to the root**, with three further tests naming the specific claims
  that drifted, because a generic comparison would not say which. Hand-patching nine files is
  how this happened; the test is the part that stops it.

  It paid for itself immediately: it failed on the first run because **19 entries in the root
  still claimed `gnomAD_v4.0`** -- the per-ancestry frequencies, which the earlier fix in this
  release had missed while correcting only two of twenty-one.

- **531 Config entries claimed a VEP version the image has never contained.** `Dockerfile`
  pins `ARG VEP_VERSION=114.0` and a real run's annotated VCF header reads `##VEP="v114.0"`,
  while 59 entries in each of the nine Config files, plus `README.md`, `docs/fields.md`,
  `docs/SMART_pipeline_metadata.md` and the generated example page, said 114.2. The Tier 2/3
  metadata row printed 114.2 directly beneath a run banner saying 114.0. Corrected, and the
  test now reads the pin out of the Dockerfile and compares, so the two cannot diverge again.

- **Three README claims broke a fresh install**, each invisible on a machine that already
  worked:

  The 91 GB SpliceAI tables were documented as needed "ONLY for `--spliceai lookup`", and the
  default `compute` mode's reference data as "10 MB". Both paths are in `entrypoint.sh`'s
  `REQUIRED_PATHS`, checked with no reference to the mode, so a user who downloaded the 10 MB
  file got `ERROR: Required resource not found`. The VEP cache was documented as
  `homo_sapiens/`; the pipeline runs `vep --merged`, which needs `homo_sapiens_merged/` -- and
  the developer's machine has both, so this only ever failed on someone else's install. And
  `--output-dir` was absent from the "Required" block although the run aborts without it,
  while the token and transcript list were listed as required although `--vep-only` needs
  neither. All three are now checked by a test that reads both files and compares.

- **Every headline number in the test documentation was wrong**, and one of them was not a
  count but a claim. "**219 columns against the sources' own APIs**" appeared in both READMEs
  and in `verify_local_sources.py`: 219 is a count of source attributions in `Config.yaml`,
  not of anything verified. The API verifier compares **47** distinct columns, and the
  file-based verifier a further **146**. It overstated API-verified coverage by **4.7x**, in
  the section a reader consults to judge how much to trust the output.

  Also corrected: "632 tests" and, in a second README for the same suite, "722" (1,033);
  the classifier's "146 unit tests" (211); the cohort as "22 records" (45). The numbers are
  now checked against the repository -- counts within a tenth, so adding three tests does not
  fail the docs, while 632-against-1,033 does; the column and cohort figures exactly, because
  they come from committed files. Past-tense notes about what the fixture *had* are left
  alone: the issue-19 write-up needs to say 22 to explain what changed.

- **No check in `local_source_results.tsv` could be traced to a variant.** All 5,986 rows
  carried `Variant_ID = "."`, because the expression was
  `row.get("Variant_ID") or row.get("dbSNP_RS") or "."` and **neither column exists in the
  MAF**, so the fallback was the only branch that ever ran. 4,996 passing checks were no
  comfort: the file exists to be read when something disagrees, and there was nothing in it
  to say which variant had disagreed. Now the MAF's `ID`, with a coordinate fallback that
  cannot be blank. 41 distinct labels where there was one.

- **`--only <suite>` aborted on the runner's last line.** `V1_DIR` was set inside the
  verification1 branch but used by the report section at the end, so a single-suite run
  reached it unset and died under `set -u` -- after the summary had printed PASS. The run
  looked successful and then failed, which is the worst of both. The paths are set
  unconditionally now, and a test asserts the assignment is unindented, which is the testable
  form of "outside any conditional".

- **A multi-allelic VCF record loses its allele fraction and every allele count, silently**
  ([#33](https://github.com/WeTGI-colab/SMART/issues/33)). Measured: `AD=55,23` gives
  VAF 0.295; `AD=30,23,25` gives **NaN** for `VAF`, `t_ref_count`, `t_alt_count` and
  `t_depth`. Empty is the right answer while the allele is unresolved -- a count taken from
  the wrong `AD` element is a plausible number attributed to the wrong variant -- but it
  happened in silence, and `VAF` is a Tier 3 column a geneticist reads.

  Nothing in the pipeline splits multi-allelic records: no `bcftools norm -m-`, nor any
  equivalent. So the row keeps `Tumor_Seq_Allele2 = "A,G"` and there is no single allele for
  it to be *about* -- which allele its HGVS, consequence and OncoKB answer describe is
  unspecified. Fixing the counts without resolving the allele would be worse than a blank.
  `vcf2table` now detects them and names them on stderr with the remedy. The code comment
  that said "see the open issue" pointed at an issue that did not exist; it does now.

  Two things found writing it: `sys` was never imported in `vcf2table.py`, so the warning
  would have raised `NameError` on exactly the input it exists to report; and the cohort
  contains no multi-allelic record at all, so the path had never run in a test. A verified
  candidate is recorded in the issue -- GNAQ `chr9:77794572 T>A,G`, REF confirmed against the
  GRCh38 FASTA, both ALTs curated by OncoKB (Q209L Oncogenic, Q209P Likely Oncogenic).

- **SMART's own OncoKB client never sent the consequence, so it under-annotated every
  variant whose interpretation depends on it.** `oncokb2.0.py` queried `byProteinChange`
  with the gene and the protein string and nothing else -- `consequence` did not appear
  anywhere in the file. OncoKB then answers only about that literal string. Measured against
  the API for TP53 `T125=`:

      alteration=T125=                                   -> hotspot False, Unknown, no Px level
      alteration=T125=&consequence=splice_region_variant -> hotspot True,  Likely Oncogenic, LEVEL_Px1
      alteration=T125=&consequence=synonymous_variant    -> hotspot False, Unknown, no Px level

  In the delivered output that cost `ONCOKB_HOTSPOT` (False where OncoKB says True),
  `ONCOKB_PROG_LVL` (empty where OncoKB says **LEVEL_Px1**, a clinical prognostic level) and
  a variant summary reading "is not annotated by OncoKB" -- all three written by this client,
  **in the same row where MafAnnotator's own columns had it right**, because MafAnnotator
  does send the consequence. One row, two clients, two answers.

  Verified not to be an endpoint difference: the batch endpoint MafAnnotator uses returns
  exactly what the single endpoint returns when both are given the consequence.

  The client now takes `Consequence` from the CSQ and sends its leading term. One term,
  because OncoKB does not understand a compound: the `%2B`-joined form returns Unknown, and
  `splice_region_variant&synonymous_variant` only appeared to work because `&` ends a query
  parameter, so OncoKB never saw the second half. VEP orders the terms of one annotation by
  severity, so the leading term is the one carrying the interpretation -- and for this
  variant the trailing one, `synonymous_variant`, is precisely the answer that was wrong.
  The consequence is part of the cache key, because the same protein string with a different
  consequence is a different question. The fallback path needed its own lookup: it returns a
  gene and protein change but not the consequence, so that is found again **by `Feature`** --
  taking the first CSQ's would sometimes describe a different transcript, and a wrong
  consequence is a wrong question that OncoKB answers. 24 tests, no network.

- **The API verifier asked OncoKB a weaker question than the pipeline, then reported the
  pipeline as wrong.** `verification1/verify.py` sent the alteration alone. Its two standing
  MISMATCHes on TP53 T125T were its own: the MAF held `Likely Oncogenic`, which is correct,
  and the verifier's call returned `Unknown`. It now sends `Variant_Classification`, the
  header `AnnotatorCore` reads, so the two queries are one question. **A check that fails
  correct output costs more than one that passes broken output**: it teaches the reader to
  ignore the two red lines that are always there.

  Fixing it is also what exposed the client defect above -- with the verifier asking properly,
  three further fields disagreed, and those three turned out to be ours.

- **The empty-column allowlist is now per table.** `post_analysis` pads the OncoKB indexed
  columns out to a fixed count, so `ONCOKB_TX_*_abstracts` is wholly populated in
  `FINAL_Table` and partly empty in the Tier 1 MAF. One list describing both was reported
  stale for the first table while being required by the second, and "remove the entry" was
  then wrong advice. `--allowlist` is repeatable; the shared list carries what is true of any
  output and a second carries what is only true later. A test asserts the second is not a
  subset of the first, so a redundant file gets folded back rather than lingering.

  The 118 columns the Tier 1 audit reported were triaged rather than waved through: 86 are
  OncoKB treatment sub-fields that fill only for an FDA-approved indication, an excluded
  cancer type or an abstract citation (22 siblings do carry data, so the path is exercised);
  22 are CIViC evidence and assertion fields this cohort does not reach, while its CIViC
  *variant* fields are populated, which locates the gap precisely; 4 are the normal-sample
  columns of a tumour-only fixture.

- **The OncoKB annotator is pinned.** `Dockerfile` cloned it with no commit, so which
  annotator an image carried was a function of its build date rather than of SMART's version.
  Measured across the published images: `smart:1.0.7` ships annotator `4bfe720` (2026-06-25)
  while `1.0.8`, `2.0.0` and `2.1.0` all ship `77521cd` (2026-08-25) -- a difference nothing
  recorded. Two images carrying the same SMART version could annotate differently, and a
  collaborator's bug report could not be reproduced, because the code they ran cannot be
  identified from the version they name. That is how [#32](https://github.com/WeTGI-colab/SMART/issues/32)
  came to describe behaviour that no current image has. Pinned at `77521cd` with the commit
  asserted after checkout, and the `|| true` that covered the optional requirements file moved
  off that assertion -- on one chain it would have left an image built from an unknown commit
  reporting success.

- **The LOEUF constraint table is verified against the gnomAD release it claims to be**
  (`scripts/verify_loeuf_against_gnomad.py`), and three things it claimed were wrong.

  Its only integrity check was `file ... | grep HTML` -- that the download was not an error
  page. It came from an **unversioned private Google Drive link**, so nothing established
  which gnomAD release the numbers were, and nothing would have noticed if that link had
  later served a different one. LOEUF and missense Z feed O8 and the Tier 3 output, so a
  substitution changes classifications with no visible symptom. A plain checksum is not
  enough here either: the distributed file is not a gnomAD artefact byte-for-byte, but the
  released TSV with its first column marked `#` and gene spans appended as `chrom/start/end`
  for tabix, so a checksum pins a file without saying whether its numbers are gnomAD's. The
  check is therefore that **every data column of every row equals a named, citable release**,
  with the index columns excluded as having no counterpart. Measured: all 55 columns of all
  **211,363 rows identical to gnomAD v4.1**. It catches a single altered `mis.z_score` among
  300 rows, distinguishes a wrong release from an unreachable one (exit 1 against exit 3,
  because the first must stop a clinical run and the second is the caller's judgement), and
  rejects a truncated reference download -- which would otherwise read as "the local file has
  extra genes", the opposite of the truth. 17 tests, no network and no reference files.

  **It was recorded as the wrong release.** `Config.yaml` said `gnomAD_v4.0`. Every one of
  the 211,363 rows matches v4.1; against v4.0 only 230 match and 210,882 differ. v4.0 also
  has a different column set, so the two are not interchangeable in any case.

  **The citation was for the wrong dataset, and the first correction was wrong too.** It
  credited Karczewski et al., Nature 2020, "the mutational constraint spectrum quantified from
  variation in 141,456 humans" -- the **v2.1.1** paper. I replaced that with Chen et al.,
  Nature 625:92-100 (2024), which is also wrong: that paper is the **Gnocchi non-coding
  1 kb-window** constraint map built on gnomAD v3.1.2 *genomes*, a different product. The
  per-gene constraint table **has no accompanying paper**; it is documented in the release
  notes and the Hail Table README, with methods descending from Samocha 2014, Lek 2016 and
  Karczewski 2020. Config.yaml now says that, and says explicitly that Chen et al. is not the
  source, so the same mistake is harder to make again.

  **And it was still described as gene-level** in `tests/verify_local_sources.py` and
  `docs/SMART_pipeline_metadata.md`, after that claim had already been corrected in
  `docs/transcript_specificity.md`. It is per transcript: 93% of multi-transcript genes carry
  different values. The verifier's comment now gives the plugin's refusal of upstream and
  downstream variants its real reason -- `match_by=transcript` needs a transcript to attribute
  the value to, and an upstream variant has none.

  The table is **verified rather than rebuilt**, deliberately. The released TSV has no
  start/end coordinates, and for its RefSeq half not even a chromosome, so rebuilding the
  index needs a gene-span source; the transcript model this pipeline ships resolves 17,904 of
  the release's 18,204 gene symbols, so a rebuild from what is on hand would drop 300 genes
  the present file has. Verifying gives the same assurance without the regression. The 160
  release rows the indexed table does not hold were checked rather than assumed: they are 61
  readthrough fusions and uncharacterised `LOC` loci, none of them cancer genes and none in
  the panel.

- **The reference manifest's contig list described genes, not chromosomes.** It was read
  from column 1 for every kind of file. That is CHROM in a VCF, but a tabix-indexed TSV puts
  the sequence name wherever its builder chose: LOEUF's is **column 56**, and column 1 is a
  gene symbol. So the fact meant to record which chromosomes a reference covers held
  **18,143 gene names**, which is how the autosome-only gap stayed invisible in the one
  artefact built to make reference properties checkable -- and why the claim that the manifest
  already recorded it, made earlier in this work, was wrong.

  The column now comes from the index itself, read before the scan. LOEUF records
  `chr1 ... chr22`, in genomic rather than lexical order so a complete list does not look
  wrong at a glance. A list longer than a plausible number of sequences is now reported as a
  misidentified column instead of being stored, so the same mistake cannot recur silently.
  The gene-keyed SVIG-UK tables, which are read whole and never queried by position, declare
  that they have no sequence column rather than being scanned for one -- without that, the
  new check called two sound files damaged.

  Side effect: the manifest shrank from **378 KB to 26 KB**, the difference being lists of
  gene names recorded as contigs.

- **The chrX gap in O8 is upstream, and upstream has since closed it**
  ([#31](https://github.com/WeTGI-colab/SMART/issues/31)). The per-gene note added earlier
  distinguished "looked up and below threshold" from "never looked up"; the open question was
  rebuilding the table to include the sex chromosomes.

  **gnomAD v4.1 -- the release in use -- publishes no chrX or chrY constraint.** Its official
  file has no sex-chromosome row, and ATRX, KDM6A, BCOR, STAG2, DDX3X, PHF6, GATA1, BTK, AR,
  DMD, F8, MECP2 and RPS4X are absent from it entirely. So the gap was never ours: our file
  is a faithful copy of a release that excludes them. gnomAD's stated reason was a missing
  model input rather than a methodological objection -- methylation data for chrX and chrY
  were unavailable when v4.0 was built, and the mutational model conditions expected counts
  on CpG methylation.

  **gnomAD v4.1.1 (2026-03-25) added them, and v4.1.2 (2026-09-23) is current.** Verified
  against the released file: 113 columns against our 55, 221,898 rows, **3,581 chrX rows and
  113 chrY rows**, and native v4-model values for every gene listed above -- ATRX
  `mis.z_score` 4.05, STAG2 6.96, SMC1A 9.29, DDX3X 5.95, BTK 5.85, KDM5C 7.33. 728 chrX MANE
  Select genes, 145 of them reaching the guideline's threshold of 3.09.

  **But it is not a file swap, and that is the finding that matters.** The upgrade
  recalibrates the autosomes too. Measured on the 17,460 genes with a MANE Select row in both
  releases, `mis.z_score` >= 3.09 is reached by **1,488 genes in v4.1 and 2,217 in v4.1.2 --
  750 genes gain O8 and 21 lose it**, none of them on a sex chromosome. gnomAD also moved
  their own LOEUF threshold to < 0.45 and publish no updated missense-Z threshold at all, so
  SVIG-UK's 3.09 has never been recalibrated against either release and that +49% is the size
  of the unexamined drift. Two X-linked genes come back flagged by gnomAD themselves and must
  not be read naively: DMD `mis.z_score` -5.43 (`outlier_mis`, `outlier_syn`) and MECP2 -2.97
  (`outlier_syn`).

  So the table is **not** upgraded in this release. Doing it is a scored-rule revalidation,
  not a download, and #31 stays open carrying that measurement.

  Two things recorded while establishing this. Backfilling chrX from **v2.1.1** -- the obvious
  move before v4.1.1 existed -- would have been wrong twice over: on the autosomal genes
  present in both, v2.1.1 reaches 3.09 for only 749 of the 1,470 genes v4.1 reaches it for, a
  49% shortfall that would have made an X-linked gene about half as likely to meet O8 as an
  autosomal one; and v2.1.1's own sex-chromosome values were computed with the
  autosome-calibrated model and no documented ploidy correction, whereas v4.1.1 is the first
  release to fit separate plateau models for chrX non-PAR and chrY non-PAR. And the `NA`
  values in the `chromosome` column are not damage: that column is populated on the Ensembl
  rows and left `NA` on the RefSeq ones, in the new release as in ours.

  gnomAD **v2.1.1** did cover both sex chromosomes, so backfilling from it is the obvious
  move, and it is the wrong one. On the autosomal genes present in both releases, v2.1.1
  reaches the guideline's missense-Z threshold of 3.09 for only **749 of the 1,470 genes
  v4.1 reaches it for -- it would miss 49%**. Dropped into the same column and read against
  the same fixed threshold, that would make an X-linked gene roughly half as likely to meet
  O8 as an autosomal one, with nothing in the output to show why: a sex-dependent bias in
  classification, invisibly. So the gap is reported and not filled. The reasoning is recorded
  in the verifier, and a test asserts the expected contig list is exactly the 22 autosomes,
  so a future file carrying chrX fails verification and forces the note and the Config
  description to be revisited rather than quietly bypassed.

- **`End_Position` was empty on every row of every delivered output.** 44 of 44 in the Tier 1
  MAF the laboratory reads, copy-number calls included -- a segment with no end cannot be
  interpreted. Three defects composed: `compute_end` computed a value only when `ID == "."`,
  so a caller that names its records (Manta does, and so does any VCF with an rsID) fell
  through to `None`; for the Manta rows it read `DUPSVLEN`, which no caller writes and which
  is not a column of that table, while `END` sat correct in the same row; and upstream
  `vcf2table` only computes its own `END` when *no* record in the file has INFO `END`, so one
  Manta record skipped the calculation for every variant.

  The data was never missing: `END` was present and right for all four copy-number rows (MET
  116798037, CDKN2A 22009812, CDK4 57793240, ERBB2 39730426) and discarded. Rewritten with an
  explicit order -- the caller's `END` first, because a symbolic ALT has no sequence to
  measure; then `SVLEN`; then `start + len(REF) - 1`, which is the VCF definition and what
  `oncokb_maf_query.py` already computed for the OncoKB query, so the two were inconsistent
  with one right; and nothing for a breakend, which is a junction rather than an interval.

  **The audit built to catch this class of defect was pointed one step too early.** It ran on
  the intermediate `FINAL_Table` MAF, written before `post_analysis`, so a column that
  `post_analysis` itself emptied passed. It now runs on the Tier 1 MAF as well, where it
  immediately reports 115 further unexamined blanks.

  **And the fix above was itself one step too late**, which a collaborator's report on an
  older release ([#32](https://github.com/WeTGI-colab/SMART/issues/32)) exposed. It was made
  in `post_analysis`, which runs *after* MafAnnotator, so the table OncoKB is queried from
  still carried `vcf2table`'s `df["End_Position"] = df["END"]` -- INFO `END`, which only
  structural variants write. Measured on the 44-variant cohort: **4 of 44 rows** reached the
  annotator with an end.

  That is not only a reporting defect. MafAnnotator resolves the query type **per row** and
  falls back to `Genomic_Change` when a row has no usable protein change -- and that fallback
  tests only that the genomic *headers* exist, not that they hold values. So for every variant
  with no protein change (promoter, canonical splice, splice region, regulatory) it built a
  query from a blank end, OncoKB answered nothing, and the row was reported
  `ONCOGENIC=Unknown`. `End_Position` is now computed in `vcf2table`, per row where INFO `END`
  is absent rather than only when the whole column is: a VCF holding one Manta record has the
  column, empty for every SNV and indel, which is how gating on the column's absence skipped
  the calculation for an entire file. 40 of 44 rows went from blank to filled, the remaining
  four being breakends, which have no span. `calculate_end`'s default branch now uses the REF
  span too, rather than returning the start position, so a row VEP gave no `VARIANT_CLASS`
  gets the VCF-defined end; the three implementations of this calculation now agree where
  previously only one was right.

- **Every frameshift was queried against OncoKB as a missense substitution at the same codon.**
  `oncokb2.0.py` built the protein change with
  `re.match(r"([A-Za-z*]{3})(\d+)([A-Za-z*]{3})")`, and `re.match` anchors at the start
  without requiring the end: `Trp288CysfsTer12` matched `Trp`+`288`+`Cys` and the `fs*12` was
  discarded. Measured against the API, `W288Cfs*12`, `W288Cfs` and `W288fs` all return
  Oncogenic with **LEVEL_Dx1**, while `W288C` returns "has not specifically been reviewed". So
  NPM1 lost the diagnostic level that defines an AML entity, and BRCA2 `S1982Rfs*22` -- a
  truncating variant in a tumour suppressor -- was reported as a VUS with "no available
  functional data about the BRCA2 S1982R mutation".

  The danger was not an unparseable query. Where the truncation produced something OncoKB
  cannot read it normalises correctly: `Arg248*`, `R248*` and `Arg248Ter` all return the same
  answer, so nonsense variants were never affected. The frameshift case was dangerous because
  the truncated string is a **valid different variant** that OncoKB answers about, so the
  wrong answer came back looking right.

  A second defect in the same expression: anything that did not match that pattern was
  returned unconverted, so `Met1673_Leu1676del` reached OncoKB with three-letter codes. Fixed
  with the whole-string replacement `vcf2table.convert_hgvsp_short` already used -- the two
  implementations had diverged and only one was right. A first attempt added lookarounds that
  were wrong in both directions, refusing the `Cys` in `CysfsTer12` because `f` follows it and
  the `Pro` in `delinsPro` because `s` precedes it. 17 tests, needing no token or network.

- **A tumour/normal VCF had its tumour values taken from the normal.** Every value SMART
  derives from a sample column describes the tumour, and it took whichever column came first
  without checking there was only one. For Strelka and SomaticSniper, which write NORMAL
  first, a variant at 55% in the tumour was reported at **VAF 0.02 with GT 0/0** -- a plausible
  subclonal call, not a visible error. The repository already documented the hazard:
  `tests/verification2/README.md` records that Mutect2 3.8 writes TUMOR first and Strelka and
  SomaticSniper write NORMAL first, in a suite that is manual and not automated.

  The tumour is now identified in order of authority (`vcf2table.select_sample_columns`):
  one sample column is the tumour (the
  routine DRAGEN input, unchanged); otherwise `##tumor_sample=` from the header (Mutect2 4.x);
  otherwise the column name against TUMOR/TUMOUR and NORMAL, which is what the three paired
  callers in the repository use; otherwise the run **refuses**, naming the columns, with
  `--tumour-column` to resolve it. Refusing is the safer failure because a wrong sample is a
  plausible low allele fraction.

  The matched normal is now carried rather than discarded, in the MAF format's own fields --
  `t_depth`, `t_ref_count`, `t_alt_count`, `n_depth`, `n_ref_count`, `n_alt_count`,
  `Matched_Norm_Sample_Barcode` -- none of which SMART emitted before, even for a tumour-only
  VCF. Plus `VAF_NORMAL`, which is the reason to carry it rather than only refuse: a fraction
  near 0.5 or 1.0 in the normal says the variant is **germline**, an interpretation the output
  could not previously support. The normal's FORMAT fields are expanded under an `N_` prefix,
  so which sample a cell describes cannot depend on which source claimed a column name first.

- **The CIViC reference file was structurally corrupt, and nothing could have noticed.** Its
  whole meta-header, the `#CHROM` line and the first data record were one physical line with
  no newlines, so the `##INFO` declarations were unparseable and `tabix 1:8013449` returned
  nothing.

  The cause was not a bad download. `civic_formating.py` drops symbolic ALTs deliberately
  (line 364, with a comment explaining that VEP's `exact` matching cannot match one) and the
  file on disk carried **480** of them, so it cannot have come from the committed builder -- it
  was written by a pre-repository version. What kept it alive is `get_ref_files.sh`, whose
  guard is **an existence test**: `if [[ -f "$CIVIC_VCF" ]]; then echo "[SKIP]"`. The correct
  builder had never been allowed to run over it.

  Rebuilt from the nightly TSV with liftover: 8 header lines where there was 1, three
  parseable `##INFO` declarations where there were none, 781 records, no symbolic ALTs, and
  `1:8013449` retrievable. The old file is kept outside the repository. **This also reframes
  [#24](https://github.com/WeTGI-colab/SMART/issues/24)**: the "38% of CIViC records carry a
  symbolic ALT" figure was a property of the stale artefact, not of what the builder produces.

- **Reference files are now verified before anything is annotated against them**
  (`scripts/verify_reference_files.py`, `generate` and `verify`). `get_ref_files.sh` performed
  no integrity checking of any kind -- no checksum, no structural check -- so a truncated file,
  a quota stub or a mis-concatenated one was accepted silently and annotated against for
  months, which is exactly what happened to CIViC.

  Per file: a SHA-256, the byte size, and the cheap structural facts that catch the failures
  that actually occur -- for a VCF that the header is more than one line, the record count, the
  contig list, index freshness and a tabix probe taken from the file's own first indexed
  record; for a FASTA that the `.fai` is present and arithmetically consistent. Wired into
  `get_ref_files.sh` after building and into `entrypoint.sh` **before annotating**, so a
  clinical run refuses a damaged reference rather than writing empty columns that look like
  results. It proved itself immediately by blocking a verification run on the corrupt CIViC
  file.

  The bootstrapping hazard is handled explicitly, because it is the obvious way to get this
  wrong: a manifest generated from a corrupt directory would record the corruption as correct,
  so `generate` writes the defects into the entry, refuses to certify the directory and exits
  non-zero, and `verify` fails an entry carrying recorded defects even when every byte matches.

- **The OncoKB annotator's ASCII fold is removed at build time**
  ([#27](https://github.com/WeTGI-colab/SMART/issues/27)). `append_annotation_to_file` ended
  with `rowstr.encode('ascii', 'ignore').decode('ascii')`, deleting every non-ASCII character
  from every row it wrote -- and it rewrites the whole MAF, so the loss was not confined to the
  columns OncoKB adds. 430 characters across 18 kinds were being removed:
  `α-ketoglutarate` in the IDH1 mechanism became `-ketoglutarate`, `PI3Kβ` became `PI3K`, and
  `a DOR of ≥6 months` became `a DOR of 6 months`, which states a different clinical claim.

  Safe to remove because the fold is the last operation before the write and feeds nothing; the
  API calls happen earlier on the parsed values. The patch asserts itself in three steps --
  target present before, absent after, file still parses -- so an upstream change fails the
  build instead of silently restoring the deletion.

  **And it exposed a defect of ours.** While that issue was open,
  `tests/verify_local_sources.py` folded *both* sides to ASCII so as not to report a difference
  the pipeline could not avoid. With the fold gone from the pipeline, keeping it would be worse
  than useless: both sides would still agree, so a real disagreement in exactly those
  characters would pass -- a test comparing two mutilated values and finding them equal, which
  is how `ONCOKB_VUS` survived for years. Left as a named identity function with the reasoning,
  and three tests asserting a difference in `≥` or `β` is caught.

- **O8 was recorded as "not met" for genes it had never been looked up for**
  ([#31](https://github.com/WeTGI-colab/SMART/issues/31)). The "not evaluated" guard tested
  whether the gnomAD constraint *table* had loaded, which it always has; the table covers 22
  autosomes and carries no chrX or chrY rows, so every X-linked tumour suppressor -- ATRX,
  KDM6A, BCOR, STAG2, DDX3X, PHF6, GATA1, BTK, AR among them -- returned nothing from the
  lookup and was scored as lacking missense constraint. A per-gene note now distinguishes
  "looked up and below threshold" from "never looked up". **A note only: no score and no code
  changes, so no classification moves.** The table itself needs rebuilding, which is #31.

- **HTTP calls in the verifiers retry with backoff and jitter.** `verification1/verify.py` made
  none -- a single attempt with a 15-second timeout -- which is why the same cohort reported
  130, then 26, then 52 errors on consecutive runs, and why the report said something
  different about the same data depending on the day. Six attempts, doubling from one second
  with a cap and jitter, honouring `Retry-After`, and not retrying a 400 or a 404, which are
  answers rather than congestion. The full suite now completes with no errors at all, and
  verification1's effective checks rose from 884 to 1,262 -- each error had been taking every
  field of that variant with it.

- **Eight more variants in the API-verification cohort, to 45 records**, each closing a class
  the suite had never contained and documented with its reason in
  `tests/verification1/README.md`: NPM1 W288fs (`Frame_Shift_Ins`, and the **LEVEL_Dx1** the
  empty-column audit had listed as unreachable), BRCA2 S1982fs (`Frame_Shift_Del`, LEVEL_1 with
  26 treatments), ASXL1 G646fs, ERBB2 Y772_A775dup (`In_Frame_Ins`), TP53 T125T (`Silent` at
  the protein level, which OncoKB calls `X125_splice`), BTK C481S (the first record on a **sex
  chromosome**, and `Resistance` with LEVEL_R1) and a constructed BRAF MNV, because
  `entrypoint.sh` warns that the released SpliceAI model silently scores no multi-nucleotide
  variant and that warning path had never run in a test. TERT `C228A` was added last, from
  [#32](https://github.com/WeTGI-colab/SMART/issues/32): it gives the cohort its first case of
  **two different ALTs at one position**, which is where a merge keyed on coordinates alone
  loses a row, and it covers the blank-protein query path that issue is about.

  Each was confirmed three ways before being added -- the REF against the GRCh38 FASTA, the
  consequence from VEP, and OncoKB's `byGenomicChange` returning the alteration the ID claims.
  All three were needed: Ensembl's recoder offered `chr4:54733167 A>T` for KIT D816V, which
  OncoKB resolves to **D820V**, and the VEP HGVS endpoint returned every minus-strand gene's
  alleles on the transcript strand -- five times in this work, the same defect as #19.

- **Two columns state, per row, whether the annotation describes the transcript that was
  asked for**: `TRANSCRIPT_REQUESTED` and `TRANSCRIPT_MATCH` (`match`, or `fallback:<acc>`).
  VEP annotates every transcript overlapping a variant and `vcf2table.py` keeps one block;
  where none of the supplied transcripts was annotated it fell back to VEP's first block --
  an arbitrary transcript, chosen by output order -- and nothing in the output said so.

  It was worse than silent: `first_nm` was initialised and never assigned, so a fallback row
  left `NM_Transcript` **empty**, which reads as "no transcript available" rather than "a
  transcript you did not ask for". Those two readings call for opposite responses. `tier: 3`,
  so both columns appear in the Tier 2 and Tier 3 sheets, and the question is answered in the
  row rather than in a separate report.

  The first implementation had three outcomes, reporting `oncokb2.0.py`'s pin as its own
  value. Real output showed 32 of 37 rows as `pinned:` -- and the pinned accessions were the
  supplied ones, because that code sets `ONCOKB_PREFERRED_TRANSCRIPT` only when
  `candidate in preferred_tx`, the same list. So a pin can never be a transcript the user did
  not ask for, and the column would have told a reader that 32 of 37 rows carried a transcript
  from somewhere else, with a genuine fallback indistinguishable among them. The test is now
  membership of the supplied list rather than the code path, with two values and a test
  asserting a third cannot reappear.

- **`docs/transcript_specificity.md`** states, source by source, which columns describe the
  supplied transcript and why the rest cannot: ClinVar, CIViC, COSMIC, Cancer Hotspots, CADD,
  gnomAD and GENIE are keyed on position and allele, and a classification, a population
  frequency or a count of observed cases has no transcript-specific form to look up. The three
  that are transcript-specific in the source and still unresolved -- REVEL (`no_match=1`),
  MutationTaster (keyed on a different transcript set) and OncoKB's genomic-change fallback --
  are given with the reason and the trade each represents. Every row of that table was read
  out of the code that performs the annotation.

- **dbNSFP predictors now describe the transcript the user supplied**
  (`--dbnsfp-transcript-match 0|1`, default 1;
  [#29](https://github.com/WeTGI-colab/SMART/issues/29)). dbNSFP holds one value per
  transcript, and the plugin was run without `transcript_match=1`, so every cell kept every
  transcript's value, `&`-separated and parallel to `Ensembl_transcriptid`, with nothing
  marking which element was the reported one. On the test cohort that was 550 cells across 25
  columns: 190 where the first element -- what a reader takes when a cell looks like a single
  value -- held a *different* value from the reported transcript's, and 143 where the reported
  transcript had no value at all while the cell showed a stranger's. For TP53 the reported
  `ENST00000269305` sat at index 23 of 28.

  Measured on a real run, 529 cells changed: 367 from a list to the reported transcript's
  single value, 162 from a list to empty because dbNSFP carries no value for that transcript.
  The only columns still holding `&` lists are `MutationTaster_pred` and `MutationTaster_score`
  -- exactly the two the plugin documents as keyed on a different transcript set and therefore
  unaffected by the flag, and nothing else.

  Turning it on was safe because these columns are read by nothing: they appear in
  `Config.yaml`, which assigns their tier and description, and in no script. No classifier rule
  consumes them, so no classification changes and no benchmark needs re-running. The only
  consumer is the person reading the sheet, which is who was being misled. The flag is reported
  in each run's header, and `--dbnsfp-transcript-match 0` restores the old behaviour.

  Two limits stated rather than discovered later: `MutationTaster` stays unresolved, and a row
  whose reported transcript dbNSFP does not carry now shows nothing instead of someone else's
  value -- DICER1's reported transcript is absent from dbNSFP entirely.

- **`verification1` compared a variant's consequences on separator and order.** VEP writes
  several consequences `&`-joined into the MAF while the REST API returns them comma-joined in
  a different order, and the comparison split on `,` alone -- so it tested a one-element set
  against a two-element one. MET `c.3028+3A>G` is `splice_donor_region_variant&intron_variant`
  in the MAF and `intron_variant,splice_donor_region_variant` from the API: the same annotation,
  reported as a mismatch. Unnoticed until a variant with more than one consequence was added to
  the cohort. Both separators are now split and the comparison is set-based; a genuine
  difference, and a subset, still fail.

- **The two TERT promoter hotspots are now in the test cohort, so
  [#15](https://github.com/WeTGI-colab/SMART/issues/15) is verified rather than asserted.**
  `svig_uk_canonical_variants.tsv` carries them as `genomic:` rows and
  `load_canonical_genomic()` returns exactly `{('5','1295113','G','A'),
  ('5','1295135','G','A')}`, which was true and demonstrable all along. What was not
  demonstrable is that the flag fires on them: `verification1` contained no TERT variant, so
  the table posted in that issue's comment -- captioned "here is the real output" -- showed a
  TERT row that the test constructs, either from a literal dict or by moving a real row to
  that position. The same comment's GNAQ row read `Silent / LOW`, which it was only because
  of the #19 strand flip.

  `chr5:1295113 G>A` (c.-124C>T) and `chr5:1295135 G>A` (c.-146C>T) are now in the fixture,
  both REF-verified against GRCh38, on `NM_198253.3` which is in the TSO500 panel. They are
  also the cohort's first 5'Flank records, a consequence class it had never exercised.

- **The API-verification cohort grew from 22 records to 35**, chosen from what the suite's own
  field coverage report said was untested: breadth of OncoKB evidence (KRAS Q61H with 15
  treatments across levels 1/2/3A/4/R1; NRAS G12D; KIT D816V for `LEVEL_R2`; PDGFRA D842V;
  IDH2 R140Q and TP53 Y220C for diagnostic and prognostic implications); the three `ONCOGENIC`
  values nothing previously exercised (ABL1 T315I `Resistance`, MLH1 V384D `Inconclusive`,
  STK11 F354L `Likely Neutral`); CDH1 R63* for `ONCOKB_VUS`; and the first intronic and splice
  records this cohort has ever held -- MET `c.3028+3A>G`, TP53 `c.376-1G>A` and `c.782+1G>A` --
  which fill `INTRON` and are the only records here that exercise SpliceAI, `ANNOTATION_FLAGS`
  and O6 at all. Field checks went from 631 to 993.

  Every record was confirmed three ways before being added: its REF against the GRCh38 FASTA,
  and OncoKB's `byGenomicChange` endpoint returning the alteration its ID claims. Both of the
  other checks were needed. Ensembl's variant recoder offered `chr4:54733167 A>T` for KIT
  D816V, which OncoKB resolves to **D820V**, a different variant -- likewise IDH2 (R149Q) and
  KRAS (S65I). And the VEP HGVS endpoint returned all three TP53 records **strand-flipped**,
  because it gives `allele_string` on the transcript strand and TP53 is minus-strand: the
  identical defect as #19, caught by the same check.

### Fixed
- **Three of `verification1`'s twenty-two records carried a REF allele that does not exist at
  that position in GRCh38, so the suite had been verifying the wrong variants**
  ([#19](https://github.com/WeTGI-colab/SMART/issues/19)). `GNAQ_Q209L` and `DICER1_E1705K` had
  their alleles written on the gene's strand rather than the genome's -- both are minus-strand
  genes -- and `EGFR_del19` was written as the deleted bases without an anchor base, with the
  last one wrong. VEP does not reject such a record: it takes the reference base from the FASTA
  and keeps the submitted ALT, so it annotated a variant nobody submitted. DICER1 came back as
  `p.E1705*` rather than `p.E1705K`, and OncoKB's genomic-change endpoint returned
  `geneExist=false` for all three, which reads identically to "OncoKB does not curate this gene".

  Corrected to `chr9:77794572 T>A`, `chr14:95094139 C>T` and
  `chr7:55174772 GGAATTAAGAGAAGCA>G`, each verified twice: the REF against the GRCh38 FASTA, and
  the resulting variant against OncoKB, which now resolves them to exactly the alterations the
  fixture names claim -- `Q209L`, `E1705K` and `E746_A750del` (the canonical 15 bp in-frame
  deletion, `c.2236_2250del`). The three records now annotate far enough to be checked at all,
  so the suite runs more checks than before and completes with no mismatch.
  `tests/run_all_verifications.sh` no longer passes `--ref-check warn` for this fixture, so a
  regression now stops the run.

- **The fixture's two committed copies held different variants.** `verification1.vcf.gz`, which
  the pipeline reads, had 22 records; the plain `verification1.vcf` had 18, missing
  `EGFR_T790M`, `KRAS_G12C`, `PIK3CA_E545K` and `TP53_R248W` -- the source of the wrong figure
  "18 variants" in earlier documentation. Both are now generated from one file and are identical,
  and a test asserts they stay that way, because regenerating one from the other is exactly how
  four variants stop being tested without anything saying so.

- **The local-source verification ran for the first time, and 212 of its 214 mismatches were
  its own comparison.** It had never executed in a full run, because `verification3`'s failure
  aborted the suite before it. Of the 214: 142 were the SpliceAI module comparing compute-mode
  output against the precomputed tables -- different sources by design, which is the whole
  point of #2 -- now skipped when `SpliceAI_pred_SYMBOL` carries the `GENE--REFSEQ` form that
  only a local run writes; 70 were a per-transcript list of sentinels (`.;.;.;.` across 28
  transcripts) reported as disagreeing with an empty cell, now normalised to empty; 2 were
  #27's non-ASCII loss, now folded on both sides identically so a real difference still fails.
  A further 24, which appeared once the TP53 splice variants were added, were the dbNSFP
  plugin's consequence gate (#28) and are now skipped with the reason named.

  Skipping the SpliceAI rows was the first answer and it was the wrong one: it took the module
  from 137 checks to none, so in the default compute mode nothing looked at the one thing that
  can actually go wrong there. `spliceai_patch_maf.py` merges the model's VCF into the table on
  CHROM+POS+REF+ALT, and a key that matched the wrong row would put one variant's deltas on
  another, with nothing to notice. The module now verifies compute mode on its own terms:
  `SpliceAI_pred_SYMBOL` carries `GENE--REFSEQ`, because SpliceAI reports its annotation's NAME
  as the symbol, so it must agree with the gene and transcript the same row reports -- a
  disagreement is a misaligned merge, not two sources differing. The deltas are then required
  to be probabilities and the positions to fall within the `-D 500` the pipeline passes.
  21 tests, one of them the misaligned merge itself.

  That first implementation then had a defect of its own, found by looking at why 96 checks
  were still being skipped: `str(row.get(col, "") or "")` turns `0.0` into `""`, because zero
  is falsy in Python, and a delta of exactly 0.00 is the most common value SpliceAI reports.
  95 real checks were being reported as "the model returned no score". It is the same mistake
  as treating NaN as a value -- confusing *a value that happens to be empty or zero* with
  *no value* -- which has now appeared twice in this work, so the tests pin it explicitly for
  `0.0`, `0`, `"0.0"`, `"0"` and `"0.00"`, with `None` still a skip.

- **An audit of every column that never carries a value, because that is the only way this
  class of defect is visible.** `ONCOKB_VUS` below is the example: an absent API key, a null
  value and a field that does not apply to a variant all arrive as the same empty cell, so a
  broken column is indistinguishable from an empty one. `scripts/audit_empty_columns.py` asks
  of each of the output's 1,148 columns whether any variant ever gave it a value, and requires
  the 106 that never do to be listed in `tests/verification1/empty_columns_allowed.tsv` with a
  reason and a kind: `structural` (76 — no variant could fill it as the pipeline runs, mostly
  the tumour-type query of #26) or `cohort` (30 — a variant would fill it and we have none).

  The `cohort` half is printed on every run as a list of what to add next. The audit fails both
  ways: an empty column that is unlisted, and a listed column that now carries data, so adding
  a variant forces the blank it filled to be re-examined rather than accumulating unexamined.
  24 tests, the logic on synthetic tables and the audit over the real output. Wired into
  `run_all_verifications.sh`.

- **`verification1` was recorded as PASS unconditionally by the suite runner.** The check
  counted rows whose **second** column matched `MISMATCH|ERROR|COVERAGE_GAP`; the second column
  is `Sample`, which never contains any of those words, so the count was always zero. That is
  why the summary said `verification1: PASS` while `verify.py` printed `RESULT: FAILED`. It now
  finds the `Result` column by name, separates a real mismatch from a transient API error, and
  folds in the empty-column audit.

- **`ONCOKB_VUS` has been `.` on every variant of every run, because the code read the wrong
  key.** `oncokb2.0.py` had `data.get("VUS")`; OncoKB's field is `vus`, lower case. A key that
  is absent returns `None`, `clean_value` turns `None` into `"."`, and `"."` is
  indistinguishable from a legitimately empty annotation -- so the column looked empty rather
  than broken. Proven against the live API: for CDH1 R63* the response is `{"vus": true}`,
  `data.get("VUS")` is `None`, and the parallel `ONCOKB_vus` column -- written by the
  flattening path, which uses the API's own key names -- held `True` all along.

  All 26 other keys the extractor reads off a response were audited against a real one and are
  correct; this was the only mismatch. `tests/verification5/test_oncokb_response_keys.py`
  (7 tests) now tokenises the source, extracts every literal key it reads, and requires each to
  exist in a frozen real response, so the next such typo fails a test instead of emptying a
  column silently. Verified by reintroducing the bug: three of the seven fail.

  `verification1/verify.py` carried the identical typo, which is why the check had passed: it
  compared an empty API value against an empty column and found them equal. Fixing only the
  pipeline turned 32 rows into mismatches, which is how the second half was found. Both sides
  now read `vus`.

- **`verification3` had been failing since SMART 2.0 landed, and its failure silently aborted
  the rest of the suite** ([#25](https://github.com/WeTGI-colab/SMART/issues/25)). Its OncoKB
  expectations were written in the repository's first commit and record 1.x behaviour: a variant
  annotated on a transcript that puts it at a different protein position was queried by that
  protein change, OncoKB had no record of it, and every field came back `Unknown`. 2.0
  deliberately replaced that with a genomic-change fallback (`780ea2a`), which returns OncoKB's
  real evidence whichever transcript was selected — correct, because `chr7:55191822 T>G` is EGFR
  L858R however it is named, and answering `Unknown` because the run called it `p.L813R` would be
  a false negative on a licensed first-line therapy. The expectations now record 2.0's behaviour,
  with the reason written next to each. Because `run_all_verifications.sh` runs under `set -e`,
  the 8 mismatches aborted the run before `verification4` and the local-source checks executed,
  so neither had been running at all. **Note for the manuscript:** transcript choice still
  changes `NM_Transcript`, `HGVSc`, `HGVSp` and anything keyed on protein position, and those
  checks still assert it, but it no longer changes OncoKB's verdict. Section 3.1.2 rests on the
  old behaviour.

- **A transient Ensembl REST failure failed `verification3` and therefore the whole suite.**
  Its exit status counted `ERROR` — a 500 or a timeout from a third-party service, several per
  run — alongside real mismatches. Errors are still reported and still written to `results.tsv`,
  but only a `SAME` or a `MISMATCH` now fails the suite, which is how `verification1` already
  behaved.

- **CIViC: the verifier blamed the wrong thing** ([#24](https://github.com/WeTGI-colab/SMART/issues/24),
  since closed). VEP annotates CIViC with `exact` matching, which requires the alleles to agree,
  and a symbolic ALT (`<inframe_deletion>`, `<N_A>`, `<frameshift_variant>` …) can therefore
  never match. The misleading note the verifier printed ("snapshot may be outdated") is
  corrected: an unmatched CIViC record is a property of the record, not of the snapshot's age.

  **Superseded within this release, and the figures above it were wrong.** This entry said 480
  of CIViC's 1,261 records carried a symbolic ALT — "over a third of it contributes nothing,
  however current the snapshot is". That was measured on a **corrupt file**: the one on disk had
  been written by a pre-repository builder, and the committed `civic_formating.py` drops symbolic
  ALTs on purpose. The rebuilt file holds **781 records and none symbolic**, so there is no 38%
  and nothing to fix in the matching. See the CIViC entry earlier in this release for what was
  actually wrong and how it stayed hidden. The reasoning against switching to `overlap` still
  stands on its own: it would give every SNV inside a large deletion that deletion's therapeutic
  evidence.

  The 7 CIViC coverage gaps the current run reports — EGFR `E746_A750del` among them — are
  variants CIViC does not hold at that position, not records it holds unmatchably.

- **Writing a string into a MAF column that was empty in every row will raise on a future
  pandas, and two such writes were already warning.** `read_csv(dtype=str)` gives an all-empty
  column `float64`, because no string ever reaches it to make it object dtype; filling it from
  another source is then a dtype change, which pandas has deprecated and will remove. This made
  a column being empty for a given cohort a latent failure rather than the ordinary case it is.
  Every column `_apply_oncokb_fixes` writes is now cast up front in one auditable place
  (`_as_object`), covering all fifteen rather than the two the test fixtures happened to expose.

- **`strip_embedded_whitespace` relied on the same deprecated downcasting.**
  `DataFrame.replace(regex=True)` across the whole frame re-infers dtypes and silently downcast
  the all-NaN columns, so this too was due to start failing. The substitution now runs only on
  the values that are actually strings, which is identical in effect -- nothing numeric can hold
  a newline -- and keeps NaN as NaN instead of the string `"nan"`. The suite now passes with
  `-W error::FutureWarning -W error::DeprecationWarning`.

- **The unreliable source of CNA treatment levels could overwrite the reliable one.** The
  structured `ONCOKB_TX_{n}_level` / `_drugs` columns are written upstream by `oncokb2.0.py` from
  the API response; the fallback that parses the stringified `ONCOKB_treatments` blob ran *after*
  them and wrote without testing whether the cell was already filled. Where the blob did parse it
  therefore replaced a correct value with whatever it yielded. It is now a fallback in fact as
  well as intent, filling only cells still empty.

- **`expand_oncokb_json_arrays` reported unparseable values as "all rows empty".** The three
  OncoKB array columns reach `post_analysis` as Python reprs whose unescaped apostrophes `json`
  rejects, not as empty cells -- on `verification1`, 15 of 18 rows for `ONCOKB_treatments`, not
  zero. Nothing was lost (the per-entry columns come from `oncokb2.0.py`, and the source column is
  kept rather than dropped), but the message said the opposite of what had happened and hid how
  widespread it is. It now counts and names the real reason.

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

- Splice O2 was applied to repeat-shifted indels that VEP flags as splice donor/acceptor
  although the 3'-shifted HGVS change lies outside the ±1/2 dinucleotide
  (e.g. a 4-bp deletion at +5..+8 inside a repeat that starts at the donor GT). This bullet
  sat under a second `### Fixed` heading of its own, an artefact of the 2.0.0 section being
  renamed rather than a separate set of changes.

### Changed
- **Test documentation brought up to date.** The README described four verification suites and
  161 unit tests; there are five, a sixth verifier for the annotation sources, and 632 unit
  tests **as of this entry** — the suite has since grown past 1,000, and those figures are now
  checked against the repository rather than restated by hand. It also gave verification1 as 18 variants, a figure taken from the stale plain-text
  copy of a fixture that exists in two differing versions ([#19](https://github.com/WeTGI-colab/SMART/issues/19));
  it is 22. `tests/verification1/README.md` recorded the same 18, and its variant table
  documented `GNAQ_Q209L` as `p.Gln209=` and `DICER1_E1705K` as `p.E1705*` — the annotations
  produced by that fixture's strand-flipped alleles, written up as though intended. Both READMEs
  now state what each suite does not cover: `verify.py` reaches only sources with a REST API,
  which is 219 of the 425 described columns, and `tests/verify_local_sources.py` covers the
  other 151. `run_all_verifications.sh` gained a `localsources` suite. Both READMEs also warn
  against reading a MAF with `comment="#"` ([#4](https://github.com/WeTGI-colab/SMART/issues/4)).

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

## [2.0.0] — superseded, no separate entry

**An image with this stamp exists**, so a MAF whose first line reads `#SMART_VERSION 2.0.0`
needs an explanation, and the preamble above promises one.

There is no 2.0.0 section because there were never two sets of changes. The section was written
as `## [2.0.0] — unreleased (in development)` and renamed to `## [2.1.0]` in commit `476bfa8`,
taking its content with it: Tier 4 SVIG-UK classification, dbNSFP, and the O2/O3/O4/O5 work are
all described under 2.1.0 above. Nothing was dropped and nothing moved to a different release.

So for a file stamped 2.0.0, read the 2.1.0 section, with one caveat: it has grown since, so it
describes **more** than that image contains. The pinned annotator commit is the one fact that
distinguishes them — `smart:2.0.0` and `smart:2.1.0` both ship oncokb-annotator `77521cd`,
while `smart:1.0.7` ships `4bfe720`. Anything else in 2.1.0's entry may postdate the 2.0.0
image. Treat results from it as provisional and re-run rather than reconcile: the 2.1.0 work
includes the `End_Position`, frameshift-query and tumour/normal fixes, each of which changes
delivered values.

The underlying mistake was renaming a version heading rather than opening a new one. Recorded
here instead of being tidied away, because the next person to rename a heading should know what
it costs.

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
