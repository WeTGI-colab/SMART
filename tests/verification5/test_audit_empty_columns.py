"""
Unit tests for scripts/audit_empty_columns.py, plus the audit itself against the real output.

`ONCOKB_VUS` carried no value on any row of any run for years, because the code read the API's
`vus` field as `VUS`. Nothing caught it: an absent key, a null value and an inapplicable field
are all an empty cell by the time they reach the MAF, and "." reads as a legitimate blank.

The audit turns that into something a test can hold. It requires every always-empty column to
be listed with a reason, and it fails in both directions:

  - a column that is empty and unlisted may be the next ONCOKB_VUS;
  - a listed column that now carries data means the entry is stale. This is the direction
    that fires when a variant is added to the cohort, which is how a blank gets examined
    instead of accumulating.

The logic is tested on synthetic tables, so these run with no Docker, no refs and no token.
The last class runs the audit over verification1's real output when it is present.
"""
import importlib.util
import os
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "..", "scripts")
_spec = importlib.util.spec_from_file_location(
    "audit_empty_columns", os.path.join(SCRIPTS, "audit_empty_columns.py"))
aec = importlib.util.module_from_spec(_spec)
sys.modules["audit_empty_columns"] = aec
_spec.loader.exec_module(aec)

REAL_MAF = os.path.join(os.path.dirname(__file__), "..", "verification1",
                        "output", "FINAL_Table", "verification1.oncoKB.maf")
ALLOWLIST = os.path.join(os.path.dirname(__file__), "..", "verification1",
                         "empty_columns_allowed.tsv")


def write_maf(tmp_path, header, rows, lead="#SMART_VERSION 2.1.0\n"):
    p = tmp_path / "t.maf"
    with open(p, "w") as fh:
        fh.write(lead)
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return str(p)


def write_allow(tmp_path, entries):
    p = tmp_path / "allow.tsv"
    with open(p, "w") as fh:
        fh.write("pattern\tkind\treason\n")
        for e in entries:
            fh.write("\t".join(e) + "\n")
    return str(p)


class TestWhichColumnsCountAsEmpty:
    def test_a_column_with_a_value_is_not_empty(self, tmp_path):
        p = write_maf(tmp_path, ["A", "B"], [["1", ""], ["2", ""]])
        h, r = aec.read_table(p)
        assert aec.always_empty(h, r) == ["B"]

    @pytest.mark.parametrize("sentinel", ["", ".", "-", "nan", "NA", "NaN", "None", "null"])
    def test_every_sentinel_counts_as_empty(self, tmp_path, sentinel):
        """'.' is MafAnnotator's blank and '-' is VEP's; a column of either is not data."""
        p = write_maf(tmp_path, ["A", "B"], [["x", sentinel], ["y", sentinel]])
        h, r = aec.read_table(p)
        assert "B" in aec.always_empty(h, r)

    def test_whitespace_only_counts_as_empty(self, tmp_path):
        p = write_maf(tmp_path, ["A", "B"], [["x", "   "], ["y", ""]])
        h, r = aec.read_table(p)
        assert "B" in aec.always_empty(h, r)

    def test_one_value_among_sentinels_is_enough(self, tmp_path):
        """The question is whether the column *ever* carries a value, not usually."""
        p = write_maf(tmp_path, ["A", "B"], [["x", "."], ["y", "True"], ["z", "."]])
        h, r = aec.read_table(p)
        assert "B" not in aec.always_empty(h, r)

    def test_false_is_a_value_not_a_blank(self, tmp_path):
        """The ONCOKB_VUS case: 'False' on every row is the column working, not failing."""
        p = write_maf(tmp_path, ["A", "ONCOKB_VUS"], [["x", "False"], ["y", "False"]])
        h, r = aec.read_table(p)
        assert "ONCOKB_VUS" not in aec.always_empty(h, r)

    def test_a_short_row_does_not_crash(self, tmp_path):
        p = tmp_path / "ragged.maf"
        p.write_text("#v\nA\tB\tC\nx\ty\nx\ty\tz\n")
        h, r = aec.read_table(str(p))
        assert aec.always_empty(h, r) == []

    def test_the_hash_inside_a_value_does_not_truncate_the_row(self, tmp_path):
        """Issue #4: a '#' mid-line is not a comment. ClinVar writes them inside values."""
        p = write_maf(tmp_path, ["A", "B"],
                      [["UniProtKB:P01111#VAR_006847", "kept"]])
        h, r = aec.read_table(p)
        assert aec.always_empty(h, r) == [], "a value after a '#' was lost"


class TestTheTwoDirectionsOfFailure:
    def test_an_unlisted_empty_column_is_reported(self, tmp_path):
        maf = write_maf(tmp_path, ["A", "MYSTERY"], [["x", "."]])
        allow = write_allow(tmp_path, [])
        h, r = aec.read_table(maf)
        unexplained, stale, _ = aec.audit(h, r, aec.load_allowlist(allow))
        assert unexplained == ["MYSTERY"]
        assert not stale

    def test_a_listed_column_that_now_has_data_is_reported_stale(self, tmp_path):
        """Fires when a variant is added that fills a previously blank column."""
        maf = write_maf(tmp_path, ["A", "LEVEL_Dx1"], [["x", "LEVEL_Dx1"]])
        allow = write_allow(tmp_path, [("LEVEL_Dx1", "cohort", "a Dx1 variant")])
        h, r = aec.read_table(maf)
        unexplained, stale, _ = aec.audit(h, r, aec.load_allowlist(allow))
        assert not unexplained
        assert [p for p, _, _ in stale] == ["LEVEL_Dx1"]

    def test_a_listed_column_absent_from_the_output_is_not_stale(self, tmp_path):
        """A column dropped by the tier split is not evidence the entry is wrong."""
        maf = write_maf(tmp_path, ["A"], [["x"]])
        allow = write_allow(tmp_path, [("GONE", "cohort", "something")])
        h, r = aec.read_table(maf)
        _, stale, _ = aec.audit(h, r, aec.load_allowlist(allow))
        assert not stale

    def test_a_glob_covers_the_indexed_columns(self, tmp_path):
        """OncoKB emits one set per treatment entry, so the allowlist needs patterns."""
        maf = write_maf(tmp_path, ["ONCOKB_TX_0_levelAssociatedCancerType.name",
                                   "ONCOKB_TX_17_levelAssociatedCancerType.name"],
                        [[".", "."]])
        allow = write_allow(tmp_path, [
            ("ONCOKB_TX_*_levelAssociatedCancerType.name", "structural", "null pan-cancer")])
        h, r = aec.read_table(maf)
        unexplained, stale, _ = aec.audit(h, r, aec.load_allowlist(allow))
        assert not unexplained and not stale

    def test_a_glob_is_case_sensitive(self, tmp_path):
        """The whole point is case: ONCOKB_VUS and ONCOKB_vus are different columns."""
        maf = write_maf(tmp_path, ["ONCOKB_vus"], [["."]])
        allow = write_allow(tmp_path, [("ONCOKB_VUS", "structural", "wrong case")])
        h, r = aec.read_table(maf)
        unexplained, _, _ = aec.audit(h, r, aec.load_allowlist(allow))
        assert unexplained == ["ONCOKB_vus"]

    def test_comments_and_the_header_are_ignored(self, tmp_path):
        p = tmp_path / "a.tsv"
        p.write_text("# a note\npattern\tkind\treason\nA\tcohort\tr\n\n# another\n")
        assert aec.load_allowlist(str(p)) == [("A", "cohort", "r")]

    def test_a_kind_other_than_the_two_is_still_loaded(self, tmp_path):
        """The audit must not silently drop a row because its kind was mistyped."""
        allow = write_allow(tmp_path, [("A", "struct", "typo in the kind")])
        assert aec.load_allowlist(allow) == [("A", "struct", "typo in the kind")]


class TestTheRealOutput:
    """The audit over verification1's actual output, which is what guards the pipeline."""

    @pytest.mark.skipif(not os.path.exists(REAL_MAF),
                        reason="verification1 has not been run (needs Docker, refs, token)")
    def test_every_empty_column_is_accounted_for(self):
        h, r = aec.read_table(REAL_MAF)
        unexplained, stale, empty = aec.audit(h, r, aec.load_allowlist(ALLOWLIST))
        assert not unexplained, (
            f"{len(unexplained)} column(s) carry no value on any of the {len(r)} variants and "
            f"are not listed in empty_columns_allowed.tsv. Each is either inapplicable -- add "
            f"it with a reason -- or broken like ONCOKB_VUS was: {unexplained}")
        assert not stale, (
            f"{len(stale)} allowlist row(s) no longer match an empty column, so those columns "
            f"now carry data and the entries must be removed: {[p for p, _, _ in stale]}")
        assert len(empty) > 0, "no empty columns at all: the audit is probably misreading"

    @pytest.mark.skipif(not os.path.exists(REAL_MAF),
                        reason="verification1 has not been run")
    def test_oncokb_vus_carries_a_value(self):
        """The specific column this whole file exists for."""
        h, r = aec.read_table(REAL_MAF)
        assert "ONCOKB_VUS" in h
        assert "ONCOKB_VUS" not in aec.always_empty(h, r), (
            "ONCOKB_VUS is empty again; oncokb2.0.py is probably reading the API's 'vus' "
            "field as 'VUS' once more")

    @pytest.mark.skipif(not os.path.exists(REAL_MAF),
                        reason="verification1 has not been run")
    def test_the_allowlist_reasons_are_written(self):
        """An entry without a reason defeats the purpose: it records a blank, not a judgement."""
        for pattern, kind, reason in aec.load_allowlist(ALLOWLIST):
            assert kind in ("structural", "cohort"), f"{pattern}: unknown kind {kind!r}"
            assert len(reason) > 15, f"{pattern}: the reason {reason!r} says too little"


class TestSeveralAllowlists:
    """The two MAFs audited do not have the same columns.

    post_analysis pads the OncoKB indexed columns out to a fixed count, so
    ONCOKB_TX_*_abstracts is wholly populated in FINAL_Table and partly empty in the
    Tier 1 MAF. One list describing both would be reported stale for the first table
    while being required by the second -- and "remove the entry" would then be wrong
    advice. So the lists combine rather than replace.
    """

    def test_patterns_from_several_lists_combine(self, tmp_path):
        a = tmp_path / "a.tsv"
        a.write_text("ALPHA\tcohort\tin the shared list\n")
        b = tmp_path / "b.tsv"
        b.write_text("BETA\tstructural\tin the table-specific list\n")
        allow = aec.load_allowlist(str(a)) + aec.load_allowlist(str(b))
        header = ["Hugo_Symbol", "ALPHA", "BETA"]
        rows = [["TP53", "", ""]]
        unexplained, stale, empty = aec.audit(header, rows, allow)
        assert unexplained == []
        assert stale == []
        assert set(empty) == {"ALPHA", "BETA"}

    def test_one_list_alone_still_reports_the_other_column(self, tmp_path):
        """The split must not become a way to lose a check: each list on its own still
        reports what it does not cover."""
        a = tmp_path / "a.tsv"
        a.write_text("ALPHA\tcohort\tonly alpha\n")
        header = ["Hugo_Symbol", "ALPHA", "BETA"]
        rows = [["TP53", "", ""]]
        unexplained, _, _ = aec.audit(header, rows, aec.load_allowlist(str(a)))
        assert unexplained == ["BETA"]

    def test_the_cli_accepts_the_option_twice(self):
        import inspect
        src = inspect.getsource(aec.main)
        assert 'action="append"' in src, "--allowlist must be repeatable"
        assert "for path in a.allowlist" in src, "every list given must be loaded"


class TestTheTierOneListIsNotARedundantCopy:
    """A table-specific list earns its existence only by holding patterns the shared
    list does not. If it ever becomes a subset, it should be folded back in."""

    def test_the_two_committed_lists_do_not_overlap(self):
        shared = os.path.join(os.path.dirname(__file__), "..", "verification1",
                              "empty_columns_allowed.tsv")
        tier1 = os.path.join(os.path.dirname(__file__), "..", "verification1",
                             "empty_columns_allowed_tier1.tsv")
        if not os.path.exists(tier1):
            pytest.skip("no Tier 1 allowlist committed")
        sp = {p for p, _, _ in aec.load_allowlist(shared)}
        tp = {p for p, _, _ in aec.load_allowlist(tier1)}
        assert tp, "the Tier 1 list is empty; delete it rather than keep an empty file"
        assert not (sp & tp), f"patterns in both lists: {sorted(sp & tp)}"

    def test_every_tier_one_entry_carries_a_reason(self):
        tier1 = os.path.join(os.path.dirname(__file__), "..", "verification1",
                             "empty_columns_allowed_tier1.tsv")
        if not os.path.exists(tier1):
            pytest.skip("no Tier 1 allowlist committed")
        for pat, kind, reason in aec.load_allowlist(tier1):
            assert kind in ("cohort", "structural"), (pat, kind)
            assert len(reason.strip()) > 20, f"{pat} has no usable reason: {reason!r}"


class TestTheRunnerDefinesItsPathsUnconditionally:
    """`--only <suite>` aborted on the runner's last line.

    V1_DIR was set inside the verification1 branch but used by the report section at the
    end, so any single-suite run reached it with the variable unset and died under
    `set -u` -- after the summary had printed PASS. The run looked successful and then
    failed, which is the worst of both.
    """

    RUNNER = os.path.join(os.path.dirname(__file__), "..", "run_all_verifications.sh")

    def test_the_paths_are_set_outside_any_conditional(self):
        """Unindented is the testable form of "unconditional": an assignment indented
        inside an `if` is one a `--only` run can skip. Asserting it comes before the
        first branch is too strong -- earlier suites have branches of their own."""
        lines = open(self.RUNNER).read().splitlines()
        for var in ("V1_DIR", "V1_OUTPUT", "V1_RESULTS"):
            assigns = [(i + 1, l) for i, l in enumerate(lines)
                       if l.split("=")[0].strip() == var and "=" in l]
            assert assigns, f"{var} is never assigned"
            top = [(n, l) for n, l in assigns if not l.startswith((" ", "\t"))]
            assert top, (
                f"every assignment of {var} is indented, so it sits inside a conditional "
                f"a --only run can skip: {assigns}")

    def test_every_use_of_them_is_after_they_are_set(self):
        src = open(self.RUNNER).read()
        for var in ("V1_DIR", "V1_OUTPUT", "V1_RESULTS"):
            decl = src.index(var + "=")
            uses = [i for i in range(len(src)) if src.startswith('"$' + var, i)
                    or src.startswith("$" + var + "/", i)]
            early = [i for i in uses if i < decl]
            assert not early, f"${var} is used at offset(s) {early} before being set"

    def test_the_runner_still_runs_under_set_u(self):
        """`set -euo pipefail` is what turns an unset variable into an abort; the guard
        above only matters while it is in force."""
        src = open(self.RUNNER).read()
        assert "set -euo pipefail" in src
