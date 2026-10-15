"""A document must not name a script, function or constant that does not exist.

`build_config.py` was documented as the tool that generates `Config.yaml`. There is no
`build_config.py` in this repository and there never was one here: `create_config/README.md`
named it six times, `scripts/post_analysis.py`'s header seven, and both told the reader to
edit `TIER2_FIELDS` / `TIER3_FIELDS` / `DROP_FIELDS` **in it** to change which columns reach
the clinical TSV. The real generator is `create_config/update_confg.py` — the typo is in the
filename. A maintainer following those instructions had nothing to open.

The same header called the file `merge_maf.py`, gave `USAGE: python merge_maf.py`, and told
the reader to "edit the entry point at the bottom of the script" to set `input_dir` /
`output_dir` / `yaml_config`; the file is `post_analysis.py`, `entrypoint.sh` runs it as
`python3 post_analysis.py --config ...`, and the config is an argument. `create_config/README.md`
also described `expand_field_patterns()`, a function that exists nowhere — the patterns are
never expanded, `fnmatch` resolves them at run time — and left a bare `**TODO**` where the
reader needed the location of the tier sets.

Two omissions of the same kind: `utils/README.md` documented 4 of 15 scripts, leaving the
five builders of the Tier 4 classifier's optional reference tables undocumented, and
`reference_lists/README.md` omitted 5 of its 12 files, all of them SVIG-UK resources.

So these tests do three things, for the documents above:

  1. every `.py` / `.sh` file they name must exist, and the two vanished names must appear
     only where the text says they are gone;
  2. every constant and function they tell a reader to go and edit must be defined where
     they say it is;
  3. every script in `utils/` and every file in `reference_lists/` must be documented,
     which is what catches the next omission rather than this one.

One group is different in kind: the tier semantics. `create_config/README.md` used to say
"Tier 1 = primary clinical fields shown by default ... Tier 3 = reserved for future use". The
code says the opposite — in `post_analysis.py`, `tier == 3` reaches both the Tier 3 and the
Tier 2 sheet, `tier == 2` reaches only the Tier 2 sheet, and every other value, tier 1
included, reaches neither TSV. A reader acting on the old text would assign clinical fields
to tier 1 and silently keep them out of both TSVs, which no test and no error message would
report. Those tests read the routing out of the code, so that if the routing ever changes it
is the documentation that is revisited rather than this file.
"""

import os
import re

import pytest

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))

# The documents under test. post_analysis.py is here for its header comment block, which is
# the module's documentation; only the first 120 lines are read, so the parser never sees code.
CREATE_CONFIG_README = "create_config/README.md"
UTILS_README = "utils/README.md"
REFERENCE_LISTS_README = "reference_lists/README.md"
POST_ANALYSIS = "scripts/post_analysis.py"

DOCS = (CREATE_CONFIG_README, UTILS_README, REFERENCE_LISTS_README)

# Names the documentation is allowed to mention *because* it says they are gone. Each is
# checked below: it must not exist, and every mention must sit beside a denial.
ABSENT_BY_DESIGN = ("build_config.py", "merge_maf.py")
DENIALS = ("there is no", "no build_config.py", "does not exist", "no longer",
           "used to", "nowhere")

# Named in the docs but written by a run rather than committed, so there is nothing to find
# on disk. `pipeline_config.sh` is the path/version record `utils/get_ref_files.sh` writes.
GENERATED_AT_RUN_TIME = ("pipeline_config.sh",)


def read(rel, limit=None):
    with open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
        text = fh.read()
    if limit is None:
        return text
    return "\n".join(text.split("\n")[:limit])


def documents():
    """(label, directory the doc lives in, text) for every document under test."""
    out = [(rel, os.path.dirname(rel), read(rel)) for rel in DOCS]
    out.append((POST_ANALYSIS + " (header)", "scripts", read(POST_ANALYSIS, limit=120)))
    return out


def repository_files():
    """Every filename in the checkout, and every path, for resolving names in the docs."""
    names, paths = set(), set()
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d != ".git"]
        for name in filenames:
            names.add(name)
            paths.add(os.path.relpath(os.path.join(dirpath, name), ROOT))
    return names, paths


FILES, PATHS = repository_files()
SCRIPT_TOKEN = re.compile(r"[\w./-]*\w\.(?:py|sh)\b")


def resolves(token, doc_dir):
    """A name in a doc resolves if it points at a file, by path or by basename."""
    if os.path.basename(token) in GENERATED_AT_RUN_TIME:
        return True
    if "/" in token:
        candidates = [os.path.normpath(os.path.join(doc_dir, token)), os.path.normpath(token)]
        if any(c in PATHS for c in candidates):
            return True
        # A bare path that resolves nowhere is still wrong if its basename is unknown too.
        return False
    return token in FILES


def denied_near(text, index):
    window = text[max(0, index - 240):index + 240].lower()
    return any(d in window for d in DENIALS)


def definitions(rel):
    """The text of a script, for asserting that a symbol is defined in it."""
    return read(rel)


class TestEveryScriptTheDocsNameExists:
    """The `build_config.py` class of error, generalised: a reader must be able to open
    every file the documentation tells them to open."""

    @pytest.mark.parametrize("label,doc_dir,text", documents(),
                             ids=[d[0] for d in documents()])
    def test_the_names_resolve_to_files(self, label, doc_dir, text):
        unresolved = []
        for match in SCRIPT_TOKEN.finditer(text):
            token = match.group(0)
            if os.path.basename(token) in ABSENT_BY_DESIGN:
                continue          # checked in TestTheVanishedNamesAreOnlyMentionedAsGone
            if not resolves(token, doc_dir):
                unresolved.append(token)
        assert not unresolved, (
            f"{label} names {sorted(set(unresolved))}, which do not exist. This is how "
            f"`build_config.py` survived in six places: nothing checked.")


class TestTheVanishedNamesAreOnlyMentionedAsGone:
    """`build_config.py` and `merge_maf.py` may be named only to say they are not here --
    both names are still searched for in issues and commit messages, so a note is useful."""

    @pytest.mark.parametrize("name", ABSENT_BY_DESIGN)
    def test_the_file_really_does_not_exist(self, name):
        """The premise. If one of these ever appears, these tests are the wrong ones."""
        assert name not in FILES, (
            f"{name} now exists; the documentation that was corrected for its absence "
            f"should be revisited rather than this test")

    @pytest.mark.parametrize("label,doc_dir,text", documents(),
                             ids=[d[0] for d in documents()])
    def test_every_mention_sits_beside_a_denial(self, label, doc_dir, text):
        for name in ABSENT_BY_DESIGN:
            for match in re.finditer(re.escape(name), text):
                assert denied_near(text, match.start()), (
                    f"{label} names {name} as if it existed, at offset {match.start()}")

    def test_nobody_is_told_to_edit_it(self):
        """The specific instruction that wasted a maintainer's time, named so it cannot
        creep back: `python build_config.py`, and "edit X in build_config.py"."""
        for label, _, text in documents():
            assert not re.search(r"python3?\s+build_config\.py", text), label
            assert not re.search(r"in\s+`?build_config\.py`?", text), label

    def test_the_post_analysis_header_names_its_own_file(self):
        header = read(POST_ANALYSIS, limit=120)
        assert "post_analysis.py" in header, "the header does not name the file it heads"
        assert "python merge_maf.py" not in header


class TestTheHeaderDescribesHowTheScriptIsActuallyRun:
    """The header told the reader to edit `input_dir` / `output_dir` / `yaml_config` at the
    bottom of the script. The config is a required argument; the two directories are not
    arguments at all but relative paths, which is why `entrypoint.sh` has to cd first. Both
    halves of that are worth stating, and both are read out of the code here."""

    @pytest.fixture(scope="class")
    def header(self):
        return read(POST_ANALYSIS, limit=120)

    @pytest.fixture(scope="class")
    def code(self):
        return read(POST_ANALYSIS)

    def test_config_is_a_required_argument(self, code, header):
        assert re.search(r'add_argument\("--config",\s*required=True', code), (
            "--config is no longer a required argument; the header documents it as one")
        assert "--config" in header

    def test_the_directories_are_still_hard_coded(self, code, header):
        """If these ever become arguments, the header's "NOT configurable" has to go."""
        assert 'input_dir="./FINAL_Table/"' in code and 'output_dir="./output/"' in code
        assert "./FINAL_Table/" in header and "./output/" in header
        assert "entry point at the bottom" not in header

    def test_the_entrypoint_invocation_matches(self, header):
        entrypoint = read("entrypoint.sh")
        assert re.search(r'python3 "\$SCRIPT_DIR/post_analysis\.py" --config', entrypoint), (
            "entrypoint.sh no longer invokes post_analysis.py this way")
        assert "python3 post_analysis.py --config" in header

    def test_the_optional_arguments_named_in_the_header_exist(self, code, header):
        for flag in ("--smart-version", "--vcf-metadata", "--oncokb-version", "--liftover"):
            if flag in header:
                assert f'add_argument("{flag}"' in code, (
                    f"the header documents {flag}, which post_analysis.py does not accept")


class TestTheSymbolsTheDocsTellYouToEditExist:
    """Curated on purpose. These are the names a reader is sent to go and change -- the tier
    sets, the paths at the top of the generator, the function that resolves a pattern. A test
    that checked every capitalised word in the documents would flag column names and report
    nothing useful."""

    # (symbol, the file the documentation says it is in)
    CONSTANTS = [
        ("DROP_FIELDS", "create_config/update_confg.py"),
        ("TIER2_FIELDS", "create_config/update_confg.py"),
        ("TIER3_FIELDS", "create_config/update_confg.py"),
        ("TIER2_ONCOKB_TX_PATTERNS", "create_config/update_confg.py"),
        ("TIER2_ONCOKB_DIAG_PATTERNS", "create_config/update_confg.py"),
        ("INPUT_YAML", "create_config/update_confg.py"),
        ("OUTPUT_YAML", "create_config/update_confg.py"),
        ("CLINVAR_VCF", "create_config/update_confg.py"),
        ("HOTSPOTS_VCF", "create_config/update_confg.py"),
        ("COLS_TO_REMOVE", POST_ANALYSIS),
    ]
    FUNCTIONS = [
        ("resolve_field_meta", POST_ANALYSIS),
        ("load_canonical_genomic", POST_ANALYSIS),
        ("load_functional_assays", "scripts/svig_uk_classifier.py"),
    ]

    @pytest.mark.parametrize("symbol,path", CONSTANTS, ids=[c[0] for c in CONSTANTS])
    def test_the_constant_is_assigned_where_the_docs_say(self, symbol, path):
        assert re.search(rf"^{symbol}\s*=", definitions(path), re.M), (
            f"{symbol} is not defined in {path}; the documentation sends the reader there")

    @pytest.mark.parametrize("symbol,path", FUNCTIONS, ids=[f[0] for f in FUNCTIONS])
    def test_the_function_is_defined_where_the_docs_say(self, symbol, path):
        assert re.search(rf"^\s*def {symbol}\(", definitions(path), re.M), (
            f"{symbol}() is not defined in {path}")

    def test_the_constants_are_named_in_the_readme_that_tells_you_to_edit_them(self):
        """Without this, the fix could be "stop naming them", which leaves the reader with
        "update the sets in the script" and no way to find them."""
        readme = read(CREATE_CONFIG_README)
        for symbol in ("DROP_FIELDS", "TIER2_FIELDS", "TIER3_FIELDS"):
            assert symbol in readme
        assert "update_confg.py" in readme

    def test_no_doc_describes_expand_field_patterns(self):
        """It never existed. `update_confg.py` assigns a tier to each pattern and leaves it a
        pattern; `fnmatch` resolves it against the run's real columns."""
        assert not re.search(r"^\s*def expand_field_patterns\(", read("create_config/update_confg.py"), re.M)
        for label, _, text in documents():
            assert "expand_field_patterns" not in text, (
                f"{label} describes expand_field_patterns(), which exists nowhere")

    def test_the_patterns_are_resolved_by_fnmatch_at_run_time(self):
        """The premise of the sentence that replaced it."""
        assert "fnmatch.fnmatch" in read(POST_ANALYSIS)
        assert "fnmatch.fnmatch" in read("scripts/vcf_metadata.py")

    def test_no_document_leaves_a_bare_todo(self):
        """`create_config/README.md` ended its tier section with a bare `**TODO**`, exactly
        where the reader needed to be told where the sets live."""
        for label, _, text in documents():
            assert not re.search(r"\*\*TODO\*\*", text), f"{label} still has a bare TODO"


class TestTheTierSemanticsMatchTheRouting:
    """The documented meaning of a tier was inverted. These read the routing out of
    `post_analysis.py` so that a change in the code fails here, not silently in a report."""

    @pytest.fixture(scope="class")
    def routing(self):
        code = read(POST_ANALYSIS)
        start = code.index("Resolve final column sets")
        return code[start:start + 1200]

    def test_tier_three_reaches_both_sheets(self, routing):
        block = routing[routing.index("elif tier == 3"):]
        block = block[:block.index("print(")]
        assert "tier3_cols.append" in block and "tier2_cols.append" in block, (
            "tier 3 no longer reaches the Tier 2 sheet as well; the READMEs say it does")

    def test_tier_two_reaches_only_the_tier_two_sheet(self, routing):
        block = routing[routing.index("if tier == 2"):routing.index("elif tier == 3")]
        assert "tier2_cols.append" in block and "tier3_cols.append" not in block

    def test_no_other_tier_reaches_a_tsv(self, routing):
        """Including tier 1. This is the trap: a clinical field given tier 1 is in the MAF
        and in neither TSV, and nothing says so at run time."""
        assert len(re.findall(r"tier(?:2|3)_cols\.append", routing)) == 3, (
            "a third route into a TSV appeared; the tier tables in the READMEs and in "
            "post_analysis.py's header describe exactly two")

    def test_the_generator_names_the_sets_for_their_audience(self):
        """Where the documentation's "Tier 3 is the clinical sheet" comes from."""
        generator = read("create_config/update_confg.py")
        assert "# --- Tier 3 (clinician) ---" in generator
        assert "# --- Tier 2 (bioinformatician) ---" in generator

    @pytest.mark.parametrize("path", [CREATE_CONFIG_README, POST_ANALYSIS])
    def test_the_inverted_wording_is_gone(self, path):
        text = read(path, limit=120 if path == POST_ANALYSIS else None).lower()
        assert "reserved for future use" not in text, (
            f"{path} still calls a tier reserved for future use; tier 3 is the clinical sheet")
        for line in text.split("\n"):
            if "tier 1" in line:
                assert "shown by default" not in line, (
                    f"{path} still presents tier 1 as the fields shown by default; tier 1 "
                    f"reaches the MAF and neither TSV")

    def test_the_readme_says_where_a_clinical_field_belongs(self):
        readme = read(CREATE_CONFIG_README)
        assert "Tier 3 is the clinical sheet" in readme, (
            "the one sentence a reader needs before assigning a tier")


class TestEveryUtilsScriptIsDocumented:
    """`utils/README.md` documented 4 of 15 scripts. The five builders of the Tier 4
    classifier's optional tables were among the eleven missing, so a maintainer had no way
    to learn that those tables exist, let alone how to build them."""

    def scripts(self):
        utils = os.path.join(ROOT, "utils")
        return sorted(name for name in os.listdir(utils)
                      if name.endswith((".py", ".sh"))
                      and os.path.isfile(os.path.join(utils, name)))

    def test_there_is_a_section_for_each(self):
        readme = read(UTILS_README)
        undocumented = [name for name in self.scripts()
                        if not re.search(rf"^#+ .*`{re.escape(name)}`", readme, re.M)]
        assert not undocumented, (
            f"utils/README.md has no section for {undocumented}. A script nobody documents "
            f"is a script nobody finds: the five Tier 4 table builders sat here unmentioned.")

    def test_the_tier_four_builders_are_the_ones_that_matter_most(self):
        """Named individually, because their absence is what made the omission costly:
        without these tables the classifier silently awards fewer evidence codes."""
        readme = read(UTILS_README)
        for name in ("build_svig_uk_resources.py", "build_genie_o4.py",
                     "build_functional_assays.py", "build_o5_references.py",
                     "build_gnomad_chip_ab.py"):
            assert name in readme
            assert os.path.exists(os.path.join(ROOT, "utils", name))

    def test_the_options_the_doc_shows_are_real_options(self):
        """Every `--flag` in the two reference documents must be accepted by something in
        the repository -- a util, a pipeline script or the entrypoint."""
        accepted = ""
        for directory in ("utils", "scripts"):
            path = os.path.join(ROOT, directory)
            for name in sorted(os.listdir(path)):
                if name.endswith((".py", ".sh")):
                    accepted += read(os.path.join(directory, name))
        accepted += read("entrypoint.sh")
        for rel in (UTILS_README, REFERENCE_LISTS_README):
            flags = set(re.findall(r"(?<![\w-])--[a-z][a-z0-9-]+", read(rel)))
            missing = sorted(f for f in flags if f not in accepted)
            assert not missing, f"{rel} documents {missing}, which nothing accepts"

    def test_the_reference_manifest_is_documented(self):
        """`get_ref_files.sh` STEP 13 writes `smart_reference_manifest.json`, and that file
        is what every later reference check consumes. The README described twelve download
        steps and not the one that records what was downloaded."""
        assert "smart_reference_manifest.json" in read(UTILS_README)
        assert "smart_reference_manifest.json" in read("utils/get_ref_files.sh")
        assert "verify_reference_files.py" in read("utils/get_ref_files.sh")


class TestEveryReferenceListIsDocumented:
    """`reference_lists/README.md` described 7 of its 12 files. The five it omitted were the
    SVIG-UK resources the Tier 4 classifier reads, one of which (`svig_uk_canonical_variants.tsv`)
    `post_analysis.py` reads as well."""

    def files(self):
        directory = os.path.join(ROOT, "reference_lists")
        return sorted(name for name in os.listdir(directory)
                      if os.path.isfile(os.path.join(directory, name)) and name != "README.md")

    def test_every_file_is_named_in_the_readme(self):
        readme = read(REFERENCE_LISTS_README)
        undocumented = [name for name in self.files() if name not in readme]
        assert not undocumented, f"reference_lists/README.md does not mention {undocumented}"

    def test_every_file_the_readme_names_exists(self):
        """The other direction: a documented file that was deleted or renamed."""
        readme = read(REFERENCE_LISTS_README)
        present = set(self.files())
        for name in set(re.findall(r"`(\w[\w.-]*\.(?:tsv|txt|md))`", readme)):
            assert name in present or name == "README.md", (
                f"reference_lists/README.md documents {name}, which is not in the directory")

    def test_the_classifier_options_it_names_are_real(self):
        classifier = read("scripts/svig_uk_classifier.py")
        for flag in ("--canonical", "--vcep-thresholds", "--gene-specific-splice",
                     "--functional-assays"):
            assert f'add_argument("{flag}"' in classifier, (
                f"{flag} is documented in reference_lists/README.md but not accepted")

    def test_the_readme_lists_exactly_the_files_the_dockerfile_copies(self):
        """Which lists reach the image decides whether a container run applies O10/B6 and
        the TP53 synonymous provision, so the README's claim has to match the Dockerfile.

        When this test was written the image carried two of the five and the README said
        so accurately. The Dockerfile now copies all five -- the two that were missing
        made a container classify the same variant differently from a checkout, with
        nothing in the output to say why. Rather than pin the expected set, this compares
        the two sources, so either one changing alone fails."""
        copied = re.findall(r"reference_lists/(\S+?\.tsv)", read("Dockerfile"))
        assert copied, "the Dockerfile no longer COPYs any reference list"
        in_image = set(copied)
        readme = read(os.path.join("reference_lists", "README.md"))
        claimed = set(re.findall(r"(svig_uk_\S+?\.tsv|oncokb_\S+?\.tsv)", readme))
        missing_from_readme = sorted(in_image - claimed)
        assert not missing_from_readme, (
            f"the Dockerfile copies {missing_from_readme} but reference_lists/README.md "
            f"does not mention them")

    def test_post_analysis_really_reads_the_canonical_list(self):
        """The README tells a reader this file has a second consumer outside the classifier."""
        code = read(POST_ANALYSIS)
        assert "svig_uk_canonical_variants.tsv" in code
        assert "svig_canonical" in code
