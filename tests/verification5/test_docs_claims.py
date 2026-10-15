"""The claims in docs/ must match the code and the run output they describe.

Every figure corrected in these pages was correct when it was written. That is the whole
problem: a doc claim that nothing checks has no way of failing, so it goes stale silently and
is then quoted — in a data submission, in a manuscript, to a user following a setup — long
after the thing it describes has moved.

What had drifted, and what each group below now holds:

  - `SMART_pipeline_metadata.md` documented SMART **1.0.1** in thirteen places, including the
    image tag and its own `head -1` verification recipe, against a `VERSION` of 2.1.0 and a
    committed MAF whose first line reads `#SMART_VERSION 2.1.0`. The document is for data
    submission, so the version is the one figure in it that must not be wrong.
  - both of its "exact command lines used" omitted `--output-dir`, which `entrypoint.sh`
    treats as fatal (`ERROR: --output-dir is required.`). Neither command ran.
  - it named the VEP cache `homo_sapiens/`. The pipeline passes `--merged`, which reads
    `homo_sapiens_merged/`; the annotated VCF's own `##VEP=` header records the merged path.
    On a machine that happens to have both -- the developer's -- this never shows. The same
    error in README is held by test_readme_matches_entrypoint.py.
  - it listed `--canonical`, `--mane` and `--cache` as key VEP flags. None of the three is
    passed (`--everything` implies the first two, `--offline` the third) and `--merged`, which
    is passed and changes which cache is read, was absent.
  - the tier column counts were 1,028/670/77 against a real 1,165/700/78.
  - `annotation_flags.md` said the verification sample "contains no row that should be
    flagged", citing a Silent GNAQ row that no longer exists -- it was a strand-flip artefact
    (#19) and the corrected fixture makes it a MODERATE missense. `add_annotation_flags()`
    flags **four** rows of that MAF. A page whose worked example claims zero where the code
    produces four is worse than no example.
  - `reference_check.md` said two tests assert the three bad records are "still present".
    They assert the opposite: the corrected alleles are present and the flipped ones are gone.
  - `transcript_specificity.md` quoted 37 of 37 rows on a supplied transcript; the cohort is
    now 45.
  - `field_provenance.md` quoted 754 OncoKB columns and 423 field definitions against 848 and
    427 + 37 patterns.

Numbers that move with the cohort are therefore recomputed here from the committed run rather
than compared against a second hardcoded copy, and the docs state them against that named run.
Everything that needs the run output skips when it is absent: `tests/verification1/output/` is
gitignored, so a fresh clone has none of it.
"""
import importlib.util
import json
import os
import re
import subprocess
import sys

import pytest

from conftest import post

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
DOCS = os.path.join(ROOT, "docs")
V1 = os.path.join(ROOT, "tests", "verification1")
OUT = os.path.join(V1, "output")

TIER1 = os.path.join(OUT, "output", "Final_result_tier1.maf")
TIER2 = os.path.join(OUT, "output", "Final_result_tier2.tsv")
TIER3 = os.path.join(OUT, "output", "Final_result_tier3.tsv")
SAMPLE_MAF = os.path.join(OUT, "FINAL_Table", "verification1.oncoKB.maf")
VCF_META = os.path.join(OUT, "logs", "verification1.vcf_metadata.json")
PROVENANCE = os.path.join(OUT, "logs", "run_provenance.txt")

NO_RUN = "verification1 has not been run (needs Docker, refs and an OncoKB token)"
needs_run = pytest.mark.skipif(not os.path.exists(TIER1), reason=NO_RUN)
needs_sample_maf = pytest.mark.skipif(not os.path.exists(SAMPLE_MAF), reason=NO_RUN)


def doc(name):
    return open(os.path.join(DOCS, name), encoding="utf-8").read()


def version():
    return open(os.path.join(ROOT, "VERSION"), encoding="utf-8").read().strip()


def entrypoint():
    return open(os.path.join(ROOT, "entrypoint.sh"), encoding="utf-8").read()


def columns(path, drop_leading_comment=False):
    """The header row of a committed output table, as a list of column names."""
    with open(path, encoding="utf-8") as fh:
        line = fh.readline()
        if drop_leading_comment and line.startswith("#"):
            line = fh.readline()
    return line.rstrip("\n").split("\t")


def collected(target):
    """How many tests a file collects, for the counts the docs quote."""
    out = subprocess.run([sys.executable, "-m", "pytest", target, "--collect-only", "-q"],
                         cwd=ROOT, capture_output=True, text=True).stdout
    m = re.search(r"(\d+) tests? collected", out)
    return int(m.group(1)) if m else None


def stated(text, pattern):
    """Every integer the text states for `pattern`, commas stripped."""
    return [int(m.replace(",", "")) for m in re.findall(pattern, text)]


class TestTheDocumentedVersion:
    """A submission document that names the wrong pipeline version misattributes the run.

    1.0.1 survived to 2.1.0 because nothing compared the two. The doc states the version in
    prose, in an image tag and in a `head -1` recipe, and the recipe disproved the prose.
    """

    DOC = "SMART_pipeline_metadata.md"

    def test_no_owned_doc_still_claims_1_0_1(self):
        for name in ("SMART_pipeline_metadata.md", "annotation_flags.md",
                     "reference_check.md", "spliceai_local.md",
                     "transcript_specificity.md", "field_provenance.md"):
            assert "1.0.1" not in doc(name), (
                f"docs/{name} still names SMART 1.0.1; VERSION is {version()}")

    def test_the_metadata_doc_states_the_current_version(self):
        text = doc(self.DOC)
        assert f"**Version:** {version()}" in text, (
            f"the submission document does not state version {version()}")

    def test_every_image_tag_in_it_is_the_current_version(self):
        tags = set(re.findall(r"monkiky/smart:([\w.]+)", doc(self.DOC)))
        assert tags, "the document no longer names a Docker image"
        assert tags == {version()}, (
            f"image tags {sorted(tags)} against VERSION {version()}; images are published "
            f"as monkiky/smart:<version>")

    def test_the_head_recipe_matches_what_the_pipeline_writes(self):
        """The doc's own verification recipe. It is the check that caught the staleness."""
        assert f"#SMART_VERSION {version()}" in doc(self.DOC)

    @needs_run
    def test_the_committed_maf_agrees(self):
        with open(TIER1, encoding="utf-8") as fh:
            first = fh.readline().strip()
        assert first == f"#SMART_VERSION {version()}", (
            f"the committed Tier 1 MAF begins {first!r}, not #SMART_VERSION {version()}")
        assert first in doc(self.DOC)


class TestTheDocumentedCommandLines:
    """`--output-dir` is fatal when absent, so a command line without it is not a command.

    Both "exact command lines used" omitted it. A reader following either one gets
    `ERROR: --output-dir is required.` and nothing points back at the document.
    """

    def test_output_dir_is_still_required_by_entrypoint(self):
        """The premise. If this stops being true, revisit the doc, not this test."""
        assert "ERROR: --output-dir is required." in entrypoint()

    def test_every_docker_run_block_passes_output_dir(self):
        text = doc("SMART_pipeline_metadata.md")
        # --entrypoint overrides entrypoint.sh entirely, so those blocks take none of its flags
        blocks = [b for b in re.findall(r"```bash\n(.*?)```", text, re.S)
                  if "docker run" in b and "--entrypoint" not in b]
        assert blocks, "the document no longer shows a docker run command"
        for b in blocks:
            assert "--output-dir" in b, (
                "a documented docker run command omits --output-dir, which entrypoint.sh "
                f"treats as fatal:\n{b}")

    @pytest.mark.parametrize("flag", ["--ref-dir", "--config", "--transcripts-file"])
    def test_the_other_required_flags_are_shown_too(self, flag):
        text = doc("SMART_pipeline_metadata.md")
        blocks = [b for b in re.findall(r"```bash\n(.*?)```", text, re.S)
                  if "docker run" in b and "ONCOKB_TOKEN" in b
                  and "--entrypoint" not in b]
        for b in blocks:
            assert flag in b, f"a full-pipeline command omits the required {flag}:\n{b}"


class TestTheDocumentedVepInvocation:
    """Flags listed as "key" must be flags the pipeline passes.

    `--canonical`, `--mane` and `--cache` were listed and are never passed; `--merged` is
    passed, decides which cache directory is read, and was missing. A reader building a
    reference set from this list downloads the wrong cache.
    """

    def vep_call(self):
        text = entrypoint()
        start = text.index("    vep \\")
        return text[start:text.index("\n\n", start)]

    @pytest.mark.parametrize("flag", ["--merged", "--everything", "--offline"])
    def test_the_flags_the_doc_lists_are_passed(self, flag):
        assert flag in self.vep_call(), f"entrypoint.sh no longer passes {flag}"
        assert f"`{flag}`" in doc("SMART_pipeline_metadata.md"), \
            f"{flag} is passed but the document does not list it"

    @pytest.mark.parametrize("flag", ["--canonical", "--mane", "--cache"])
    def test_the_flags_it_does_not_pass_are_not_listed_as_passed(self, flag):
        call = self.vep_call()
        assert flag not in call, (
            f"entrypoint.sh now passes {flag}; the document says it does not and should be "
            f"corrected rather than this test relaxed")
        text = doc("SMART_pipeline_metadata.md")
        m = re.search(r"\| Key flags \| (.*?) \|\n", text)
        assert m, "the Key flags row is gone from the document"
        assert flag not in m.group(1), f"{flag} is listed as a key flag but is never passed"

    def test_the_cache_is_documented_as_merged(self):
        """`--merged` reads homo_sapiens_merged/. Documenting homo_sapiens/ sends a user to
        the wrong 25 GB download, and only a machine with both caches hides it."""
        assert "--merged" in self.vep_call()
        text = doc("SMART_pipeline_metadata.md")
        assert "homo_sapiens_merged" in text
        for line in text.splitlines():
            if "homo_sapiens" not in line:
                continue
            bare = re.sub(r"homo_sapiens_merged", "", line)
            if "homo_sapiens" in bare:
                assert "not" in bare, (
                    f"this line names the unmerged cache without saying it is the wrong "
                    f"one:\n{line}")

    @needs_run
    def test_the_annotated_vcf_header_agrees(self):
        vcf = os.path.join(OUT, "AnnotatedVcf", "verification1_annotated.vcf")
        if not os.path.exists(vcf):
            pytest.skip("the annotated VCF was not kept (--keep-tmp)")
        line = next(l for l in open(vcf, encoding="utf-8") if l.startswith("##VEP="))
        cache = re.search(r'cache="([^"]+)"', line).group(1)
        assert "homo_sapiens_merged" in cache, (
            f"the run used cache {cache}, which the document does not describe")
        vep = re.search(r'##VEP="v([\d.]+)"', line).group(1)
        assert f"VEP | {vep} |" in doc("SMART_pipeline_metadata.md"), (
            f"the run used VEP {vep}; the document's Software Versions table disagrees")


class TestTheDocumentedTierColumnCounts:
    """1,028/670/77 against a real 1,165/700/78.

    The counts move with the cohort, because post_analysis pads the OncoKB evidence columns
    out to the widest row in the run. So the document states them against a named run and
    this reads that run rather than holding a second copy of the numbers.
    """

    @needs_run
    @pytest.mark.parametrize("path,lead,label", [
        (TIER1, True, "Final_result_tier1.maf"),
        (TIER2, False, "Final_result_tier2.tsv"),
        (TIER3, False, "Final_result_tier3.tsv"),
    ])
    def test_the_count_in_the_doc_is_the_count_in_the_file(self, path, lead, label):
        n = len(columns(path, drop_leading_comment=lead))
        text = doc("SMART_pipeline_metadata.md")
        row = next((l for l in text.splitlines() if label in l and "columns" in l), None)
        assert row, f"the Output Files table no longer describes {label}"
        claimed = stated(row, r"([\d,]+) columns")
        assert claimed == [n], (
            f"{label} has {n:,} columns; the document's row says {claimed}")

    @needs_run
    def test_tier_1_is_the_widest_and_tier_3_the_narrowest(self):
        """The tiers are a progressive narrowing; a doc describing them otherwise is wrong
        about what the files are for, not just about a number."""
        n1 = len(columns(TIER1, drop_leading_comment=True))
        n2, n3 = len(columns(TIER2)), len(columns(TIER3))
        assert n1 > n2 > n3, f"tier widths {n1}, {n2}, {n3} are not decreasing"


class TestTheFlaggedRowsInAnnotationFlags:
    """The page's worked example claimed zero flagged rows; the code produces four.

    It was true once. The sample's only LOW row was a GNAQ variant that VEP called synonymous
    because the fixture's alleles were strand-flipped (#19); correcting that made it a
    MODERATE missense, and the cohort has since gained the MET, TP53 and TERT records that
    are flagged. A worked example that disagrees with the code teaches the wrong filter.
    """

    def real(self):
        import pandas as pd
        return pd.read_csv(SAMPLE_MAF, sep="\t", low_memory=False)

    @needs_sample_maf
    def test_the_documented_flagged_count_is_what_the_code_produces(self):
        out = post.add_annotation_flags(self.real())
        flagged = out["ANNOTATION_FLAGS"].fillna("").astype(str).str.len() > 0
        n = int(flagged.sum())
        text = doc("annotation_flags.md")
        assert stated(text, r"flags \*\*(\d+) of the \d+ rows\*\*") == [n], (
            f"add_annotation_flags sets a flag on {n} row(s) of the committed MAF; "
            f"docs/annotation_flags.md states something else")

    @needs_sample_maf
    def test_the_documented_shape_of_that_maf_is_its_shape(self):
        df = self.real()
        text = doc("annotation_flags.md")
        assert stated(text, r"([\d,]+) rows, [\d,]+ columns") == [len(df)]
        assert stated(text, r"[\d,]+ rows, ([\d,]+) columns") == [df.shape[1]]

    @needs_sample_maf
    def test_every_row_the_doc_names_is_actually_flagged(self):
        """The table of four rows, read back out of the output. A count alone would pass
        with four different rows flagged for four different reasons."""
        out = post.add_annotation_flags(self.real())
        flagged = out[out["ANNOTATION_FLAGS"].fillna("").astype(str).str.len() > 0]
        text = doc("annotation_flags.md")
        for _, r in flagged.iterrows():
            locus = f"{r['Hugo_Symbol']} {r['Chromosome']}:{r['Start_Position']}"
            assert locus in text, f"{locus} is flagged but the document does not list it"

    @needs_sample_maf
    def test_no_silent_row_is_described_in_a_cohort_that_has_none(self):
        df = self.real()
        assert "Silent" not in set(df["Variant_Classification"].astype(str)), (
            "a Silent row is back in the cohort; the document says there is none")
        assert "Silent GNAQ" not in doc("annotation_flags.md"), (
            "the document still cites the Silent GNAQ row that #19 removed")

    @needs_run
    def test_the_documented_spliceai_column_positions_are_real(self):
        cols = columns(TIER1, drop_leading_comment=True)
        idx = [i for i, c in enumerate(cols, 1) if c.startswith("SpliceAI_pred_DS_")]
        assert len(idx) == 4, f"expected four SpliceAI delta columns, found {len(idx)}"
        text = doc("annotation_flags.md")
        assert stated(text, r"columns ([\d,]+)–[\d,]+ of") == [idx[0]]
        assert stated(text, r"columns [\d,]+–([\d,]+) of") == [idx[-1]]
        assert stated(text, r"of the ([\d,]+) in the Tier 1 MAF") == [len(cols)]


class TestTheDocumentedTestCounts:
    """A page that says "(42 tests)" of a file holding 45 is the mildest of these errors and
    the most common. Exact rather than tolerant: the figure is one `--collect-only` away."""

    @pytest.mark.parametrize("page,pattern,target", [
        ("reference_check.md", r"test_check_reference\.py` \((\d+) tests\)",
         "tests/verification5/test_check_reference.py"),
        ("annotation_flags.md", r"test_annotation_flags\.py` \((\d+) tests\)",
         "tests/verification5/test_annotation_flags.py"),
        ("transcript_specificity.md", r"test_transcript_provenance\.py` \((\d+) tests\)",
         "tests/verification5/test_transcript_provenance.py"),
        ("spliceai_local.md", r"test_spliceai_lookup_agreement\.py`?\s*\n?\((\d+) tests\)",
         "tests/verification5/test_spliceai_lookup_agreement.py"),
    ])
    def test_the_stated_count_is_what_the_file_collects(self, page, pattern, target):
        claimed = stated(doc(page), pattern)
        assert claimed, f"docs/{page} no longer states a test count for {target}"
        actual = collected(target)
        assert actual, f"could not collect {target}"
        for n in claimed:
            assert n == actual, (
                f"docs/{page} says {n} tests in {os.path.basename(target)}; it collects "
                f"{actual}")


class TestTheReferenceCheckTestsAreDescribedCorrectly:
    """The page said two tests assert the bad records are "still present". They assert the
    opposite, and the same page's own body says the records were corrected -- so it
    contradicted itself, in the direction that tells a reader the fixture is still broken."""

    def test_the_two_tests_exist_under_the_names_the_doc_gives(self):
        src = open(os.path.join(ROOT, "tests", "verification5", "test_check_reference.py"),
                   encoding="utf-8").read()
        text = doc("reference_check.md")
        for name in ("test_the_fixture_carries_the_corrected_alleles",
                     "test_the_strand_flipped_alleles_are_gone"):
            assert f"def {name}(" in src, f"{name} no longer exists"
            assert name in text, f"docs/reference_check.md does not name {name}"

    def test_the_doc_does_not_say_the_defect_is_still_present(self):
        text = doc("reference_check.md")
        assert "still\npresent" not in text and "still present" not in text, (
            "the page says the three bad records are still present; the tests assert the "
            "correction and the page's own body says they were corrected")

    def test_the_corrected_alleles_the_doc_prints_are_the_ones_asserted(self):
        """Both files name the three records' alleles. They must be the same three."""
        src = open(os.path.join(ROOT, "tests", "verification5", "test_check_reference.py"),
                   encoding="utf-8").read()
        block = src[src.index("CORRECTED = {"):src.index("PRE_FIX = {")]
        text = doc("reference_check.md")
        for chrom, pos, ref, alt in re.findall(
                r'\("(chr[\dXY]+)",\s*(\d+),\s*"(\w+)",\s*"(\w+)"\)', block):
            assert f"{chrom}:{pos} {ref}>{alt}" in text, (
                f"the corrected record {chrom}:{pos} {ref}>{alt} is asserted by the test "
                f"but not stated by the page")


class TestTheTranscriptProvenanceFigure:
    """37 of 37 against a cohort of 45. The script computes it, so this runs the script."""

    @needs_run
    def test_the_documented_figure_is_what_the_script_reports(self):
        out = subprocess.run(
            [sys.executable, os.path.join(ROOT, "scripts", "transcript_provenance.py"),
             "--maf", TIER1,
             "--transcripts", os.path.join(V1, "verification1_transcripts.txt")],
            cwd=ROOT, capture_output=True, text=True)
        m = re.search(r"selected transcript is one the user supplied: (\d+)/(\d+)", out.stdout)
        assert m, f"transcript_provenance.py printed no supplied-transcript line:\n{out.stdout}"
        used, total = m.group(1), m.group(2)
        text = doc("transcript_specificity.md")
        assert f"**{used} of {total} rows on a supplied transcript**" in text, (
            f"the script reports {used} of {total}; docs/transcript_specificity.md states "
            f"something else")

    def test_the_page_does_not_claim_the_fallback_is_unreported(self):
        """Lines 48-62 read as a live contradiction: "nothing in the output says so",
        then "Two columns now state it". It is narrative order, so the page must say which
        half describes the past."""
        text = doc("transcript_specificity.md")
        assert "nothing in the output says\nso" not in text, (
            "the page still asserts in the present tense that nothing in the output reports "
            "a fallback, two paragraphs before describing the columns that report it")
        assert "TRANSCRIPT_REQUESTED" in text and "TRANSCRIPT_MATCH" in text


class TestTheFieldProvenanceFigures:
    """754 OncoKB columns and 423 field definitions, against 848 and 427 + 37 patterns.

    The 475-of-807 figure is genuinely historical -- it is the measurement that motivated
    reading the VCF headers -- so the page keeps it and labels it. The current figures are
    recomputed here through the same resolver the pipeline uses.
    """

    def config(self):
        import yaml
        return yaml.safe_load(open(os.path.join(ROOT, "Config.yaml"), encoding="utf-8"))

    def test_the_config_counts_are_current(self):
        cfg = self.config()
        text = doc("field_provenance.md")
        assert stated(text, r"explicit field definitions \| ([\d,]+)") == [len(cfg["fields"])]
        assert stated(text, r"wildcard patterns \| ([\d,]+)") == [len(cfg["field_patterns"])]

    def test_the_historical_figures_are_labelled_as_history(self):
        text = doc("field_provenance.md")
        i = text.index("475 had no entry")
        assert "history, not the current" in text[:i], (
            "the 475-of-807 measurement is stated without saying it is historical, so it "
            "reads as the current coverage")

    @needs_run
    def test_the_current_resolution_counts_are_what_the_resolver_produces(self):
        spec = importlib.util.spec_from_file_location(
            "vcf_metadata_docs", os.path.join(ROOT, "scripts", "vcf_metadata.py"))
        vm = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(vm)

        cfg = self.config()
        meta = json.load(open(VCF_META, encoding="utf-8"))
        oncokb_version = None
        if os.path.exists(PROVENANCE):
            m = re.search(r"OncoKB data:\s*(.+)", open(PROVENANCE, encoding="utf-8").read())
            oncokb_version = m.group(1).strip() if m else None

        cols = columns(TIER1, drop_leading_comment=True)
        sources = [vm.resolve(c, meta, cfg["fields"], cfg["field_patterns"],
                              oncokb_version).get("source", "") for c in cols]
        oncokb = sum(1 for s in sources if s == "OncoKB")
        unresolved = sum(1 for s in sources if not s)

        text = doc("field_provenance.md")
        assert stated(text, r"\| Tier 1 columns \| ([\d,]+) \|") == [len(cols)]
        assert stated(text, r"the MafAnnotator table\) \| ([\d,]+) \|") == [oncokb]
        assert stated(text, r"no source at all \| ([\d,]+) \|") == [unresolved]
        assert stated(text, r"named annotation source[^|]*\| ([\d,]+) \|") == \
            [len(cols) - oncokb - unresolved]

    def test_the_oncokb_majority_claim_is_stated_against_a_named_run(self):
        """The count grows with the cohort, so an unqualified figure is the next 754."""
        text = doc("field_provenance.md")
        m = re.search(r"resolve\s+to\s+source `OncoKB`", text)
        assert m, "the page no longer states how many columns resolve to OncoKB"
        i = m.start()
        assert "tests/verification1" in text[max(0, i - 300):i], (
            "the OncoKB column count is stated without naming the run it was measured on")


class TestSpliceaiLocalDescribesTheShippedDefault:
    """The page read as an unmerged experiment -- "Exploratory work on branch
    spliceai-local-compute", "Nothing here changes the pipeline", "**Not decided**" -- while
    describing what every run has done by default since `SPLICEAI_MODE="compute"` landed.
    README links it as the authoritative detail page, so a reader was told the default is
    something it is not. The measurements and the comparison are the valuable part and stay.
    """

    def test_compute_is_still_the_default(self):
        """The premise. If the default reverts, this page's framing must revert with it."""
        assert 'SPLICEAI_MODE="compute"' in entrypoint()

    def test_the_page_does_not_call_itself_undecided(self):
        text = doc("spliceai_local.md")
        for phrase in ("Nothing here changes the pipeline",
                       "**Not decided**",
                       "Exploratory work on branch"):
            assert phrase not in text, (
                f"docs/spliceai_local.md still says {phrase!r}, which contradicts "
                f'SPLICEAI_MODE="compute"')

    def test_the_page_says_compute_is_the_default(self):
        text = doc("spliceai_local.md")
        assert "This is what SMART does by default" in text
        assert 'SPLICEAI_MODE="compute"' in text

    def test_the_deleted_branch_is_not_presented_as_live(self):
        """The branch was merged and deleted. A page pointing at it sends a reader nowhere."""
        text = doc("spliceai_local.md")
        i = text.find("spliceai-local-compute")
        assert i != -1, "the branch is no longer mentioned at all; nothing to check"
        assert "merged and deleted" in text[i:i + 400], (
            "the page names the branch without saying it has been merged and deleted")

    def test_the_image_really_ships_what_the_page_says_it_does(self):
        dockerfile = open(os.path.join(ROOT, "Dockerfile"), encoding="utf-8").read()
        text = doc("spliceai_local.md")
        for pin in ("spliceai==1.3.1", "tensorflow-cpu==2.15.*"):
            assert pin in dockerfile, f"the Dockerfile no longer installs {pin}"
            assert pin in text, f"the page does not state the shipped pin {pin}"

    def test_the_open_questions_are_kept_rather_than_dropped(self):
        """Rewriting the framing must not quietly promote the caveats to solved. #21 is open
        and the 91 GB is still required on disk, both of which a reader needs."""
        text = doc("spliceai_local.md")
        assert "## Open questions" in text
        assert "issues/21" in text
        assert "REQUIRED_PATHS" in text

    def test_the_tables_are_still_required_on_disk(self):
        """The claim behind that caveat: the 91 GB saving is not yet realised, because the
        existence check runs before anything and ignores SPLICEAI_MODE."""
        text = entrypoint()
        block = text[text.index("REQUIRED_PATHS=("):]
        block = block[:block.index("done")]
        assert "SPLICEAI_SNV" in block and "SPLICEAI_INDEL" in block
        assert "SPLICEAI_MODE" not in block, (
            "the requirement is now mode-dependent; the page's Cost caveat should be removed "
            "rather than this test relaxed")
