"""The numbers in the READMEs must match what the repository actually contains.

Every headline figure in the test documentation had drifted, some by a factor of two:
"632 tests" against 1,033; "722 tests" in a second README for the same suite; the
cohort as "22 records" against 45; the classifier's "146 unit tests" against 211.

The worst was not a count but a claim: **"219 columns against the sources' own APIs"**.
219 is a count of source attributions in Config.yaml, not of anything verified. The API
verifier compares **47** distinct columns. It overstated API-verified coverage by 4.7x,
in the README section a reader consults to decide how much to trust the output.

Counts are compared with a tolerance, because a suite that grows by three tests should
not fail its own documentation; staleness of the kind above is far outside it. The
cohort and column figures are exact, because they come from committed files.
"""

import csv
import os
import re
import subprocess
import sys

import pytest

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
TOLERANCE = 0.10          # a tenth: catches 632-vs-1033, forgives adding a few tests
csv.field_size_limit(10 ** 7)


def documented_numbers(path, pattern):
    """Every integer the docs state for `pattern`, commas stripped."""
    text = open(os.path.join(ROOT, path), encoding="utf-8").read()
    return [int(m.replace(",", "")) for m in re.findall(pattern, text)]


def collected_tests(target):
    out = subprocess.run([sys.executable, "-m", "pytest", target, "--collect-only", "-q"],
                         cwd=ROOT, capture_output=True, text=True).stdout
    m = re.search(r"(\d+) tests? collected", out)
    return int(m.group(1)) if m else None


def close(documented, actual):
    return abs(documented - actual) <= max(1, actual * TOLERANCE)


class TestTheSuiteSizeIsStatedHonestly:
    @pytest.mark.parametrize("path,pattern", [
        ("README.md", r"\*\*([\d,]+) tests\*\*"),
        ("tests/README.md", r"the logic, ([\d,]+) tests"),
        ("tests/verification5/README.md", r"\*\*([\d,]+) tests, all passing"),
    ])
    def test_the_documented_count_is_current(self, path, pattern):
        claimed = documented_numbers(path, pattern)
        assert claimed, f"{path} no longer states a test count in the expected form"
        actual = collected_tests("tests/verification5")
        assert actual, "could not collect the suite"
        for n in claimed:
            assert close(n, actual), (
                f"{path} says {n:,} tests; the suite collects {actual:,}. "
                f"Update the doc -- a reader uses this to judge coverage.")

    def test_the_two_readmes_do_not_disagree_with_each_other(self):
        """`tests/README.md` said 722 while `tests/verification5/README.md` said 632,
        for the same suite. Two wrong numbers are worse than one, because neither can
        be trusted even as an order of magnitude."""
        a = documented_numbers("tests/README.md", r"the logic, ([\d,]+) tests")
        b = documented_numbers("tests/verification5/README.md",
                               r"\*\*([\d,]+) tests, all passing")
        if a and b:
            assert close(a[0], b[0]), f"{a[0]:,} vs {b[0]:,} for the same suite"


class TestTheVerifiedColumnCountsAreNotInflated:
    """These come from the committed results files, so they are exact."""

    def distinct_columns(self, results):
        path = os.path.join(ROOT, "tests", "verification1", results)
        if not os.path.exists(path):
            pytest.skip(f"{results} has not been generated")
        with open(path, encoding="utf-8") as fh:
            return {r["MAF_Field"] for r in csv.DictReader(fh, delimiter="\t")}

    @pytest.mark.parametrize("path", ["README.md", "tests/README.md"])
    def test_the_api_column_count_matches_the_results_file(self, path):
        claimed = documented_numbers(path, r"\*\*(\d+) columns\*\* against the sources")
        actual = len(self.distinct_columns("results.tsv"))
        assert claimed, (
            f"{path} no longer states the API column count. It previously said 219, "
            f"which was a Config.yaml attribution count, not a verified-column count.")
        for n in claimed:
            assert n == actual, f"{path} says {n}, the results file holds {actual}"

    @pytest.mark.parametrize("path", ["README.md", "tests/README.md"])
    def test_the_file_verified_count_matches(self, path):
        claimed = documented_numbers(path, r"a further \*\*(\d+) columns\*\*")
        actual = len(self.distinct_columns("local_source_results.tsv"))
        assert claimed, f"{path} no longer states the file-verified column count"
        for n in claimed:
            assert n == actual, f"{path} says {n}, the results file holds {actual}"

    def test_219_appears_nowhere_as_a_verified_count(self):
        """The specific inflated figure, named so it cannot creep back."""
        for path in ("README.md", "tests/README.md", "tests/verify_local_sources.py"):
            text = open(os.path.join(ROOT, path), encoding="utf-8").read()
            for m in re.finditer(r"219", text):
                window = text[max(0, m.start() - 90):m.start() + 90]
                assert "against the sources" not in window and "APIs" not in window, (
                    f"{path} still presents 219 as a verified-column count")


class TestTheCohortSizeIsStatedExactly:
    def actual_records(self):
        path = os.path.join(ROOT, "tests", "verification1", "verification1.vcf")
        return sum(1 for l in open(path, encoding="utf-8")
                   if l.strip() and not l.startswith("#"))

    @pytest.mark.parametrize("path", ["tests/verification1/README.md"])
    def test_the_record_count_is_right(self, path):
        actual = self.actual_records()
        claimed = documented_numbers(path, r"\| (\d+) variants \(")
        claimed += documented_numbers(path, r"the same (\d+) records")
        assert claimed, f"{path} no longer states the cohort size"
        for n in claimed:
            assert n == actual, f"{path} says {n} records; the fixture holds {actual}"

    def test_no_readme_still_says_twenty_two_records(self):
        """The cohort grew from 22 to 45 in stages, and each stage left the figure
        behind in a different file."""
        for path in ("README.md", "tests/README.md", "tests/verification1/README.md"):
            text = open(os.path.join(ROOT, path), encoding="utf-8").read()
            for line in text.split("\n"):
                if not re.search(r"\b22 (records|synthetic|variants)\b", line):
                    continue
                # A dated note about what the fixture *had* is history, not staleness --
                # the #19 write-up needs to say 22 to explain what changed.
                if line.lstrip().startswith(">") or re.search(r"\bhad\b", line):
                    continue
                pytest.fail(f"{path} still says 22 in the present tense: {line.strip()[:90]!r}; "
                            f"the fixture holds {self.actual_records()}")
