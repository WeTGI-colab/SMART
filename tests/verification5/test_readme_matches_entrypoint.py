"""README must not tell a new user something `entrypoint.sh` contradicts.

Three claims in README broke a fresh install, each in the same way: a user follows the
documented setup, the run aborts, and nothing in the message points back at the doc.

  - the 91 GB SpliceAI tables were documented as needed "ONLY for --spliceai lookup",
    and the default `compute` mode's reference data as "10 MB". Both paths are in
    `REQUIRED_PATHS`, checked with no reference to the mode, so a user who downloaded
    only the 10 MB file got `ERROR: Required resource not found`.
  - the VEP cache was documented as `homo_sapiens/`. The pipeline runs `vep --merged`,
    which needs `homo_sapiens_merged/`. On a machine that happens to have both -- the
    developer's -- this never shows.
  - `--output-dir` was absent from the "Required" block although the run aborts without
    it, while the token and the transcript list were listed as required although
    `--vep-only` needs neither.

These tests read both files and compare. A doc claim nobody can check is how all three
survived.
"""

import os
import re

import pytest

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))


@pytest.fixture(scope="module")
def readme():
    return open(os.path.join(ROOT, "README.md"), encoding="utf-8").read()


@pytest.fixture(scope="module")
def entrypoint():
    return open(os.path.join(ROOT, "entrypoint.sh"), encoding="utf-8").read()


def required_paths(entrypoint):
    m = re.search(r"REQUIRED_PATHS=\((.*?)\)", entrypoint, re.S)
    assert m, "REQUIRED_PATHS not found in entrypoint.sh"
    return set(re.findall(r"\$\{?(\w+)\}?", m.group(1)))


class TestTheSpliceaiTablesAreNotDocumentedAsOptional:
    def test_they_are_in_required_paths(self, entrypoint):
        """The premise. If this ever stops being true, the README claim becomes correct
        and these tests should be revisited rather than the README."""
        req = required_paths(entrypoint)
        assert "SPLICEAI_SNV" in req and "SPLICEAI_INDEL" in req

    def test_the_check_does_not_depend_on_the_mode(self, entrypoint):
        block = entrypoint[entrypoint.index("REQUIRED_PATHS=("):]
        loop_end = block.index("done")
        assert "SPLICEAI_MODE" not in block[:loop_end], (
            "the requirement is now mode-dependent; update the README and this test")

    def test_the_readme_does_not_call_them_lookup_only(self, readme):
        assert "ONLY for --spliceai lookup" not in readme, (
            "they are required in both modes")

    def test_the_readme_says_so_where_the_10_mb_figure_appears(self, readme):
        """The figure itself is right -- compute's own reference data is 10 MB. What
        misled was the absence of a note that the 91 GB tables are still checked for."""
        i = readme.find("**10 MB**")
        assert i != -1
        window = readme[i:i + 400]
        assert "still checked" in window or "still required" in window


class TestTheCacheDirectoryMatchesTheVepInvocation:
    def test_vep_runs_with_merged(self, entrypoint):
        assert "--merged" in entrypoint

    def test_the_readme_names_the_merged_cache(self, readme):
        assert "homo_sapiens_merged" in readme

    def test_the_readme_does_not_name_the_plain_cache_as_the_requirement(self, readme):
        """`homo_sapiens/` alone will not satisfy `--merged`, and the developer machine
        has both, so this failure only ever appears on someone else's install."""
        for claim in ("(`homo_sapiens/114_GRCh38/`)",
                      "└── homo_sapiens/ "):
            assert claim not in readme, f"README still points at {claim}"


class TestTheRequiredArgumentsMatch:
    def test_output_dir_aborts_the_run(self, entrypoint):
        assert "--output-dir is required" in entrypoint

    def test_output_dir_is_listed_as_required(self, readme):
        i = readme.index("Required:")
        block = readme[i:i + 600]
        assert "--output-dir" in block, (
            "the run aborts without it, so it belongs in the Required block")

    def test_ref_dir_is_listed_as_required(self, readme):
        i = readme.index("Required:")
        assert "--ref-dir" in readme[i:i + 600]

    def test_the_conditional_ones_are_marked_conditional(self, readme):
        """The token and the transcript list are required only without --vep-only
        (entrypoint.sh checks them inside that branch), so listing them flatly as
        required overstates what a VEP-only run needs."""
        i = readme.index("Required:")
        block = readme[i:i + 600]
        assert "vep-only" in block, (
            "say which of these are conditional, or a --vep-only user is told to supply "
            "a token they do not need")
