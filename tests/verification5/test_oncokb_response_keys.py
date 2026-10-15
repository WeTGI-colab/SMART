"""
Every key scripts/oncokb2.0.py reads off an OncoKB response must exist in one.

`out["vus"] = clean_value(data.get("VUS"))` read a key the API does not have -- it is `vus`,
lower-case -- so `ONCOKB_VUS` was `.` on every variant of every run. Nothing caught it, and
nothing could have: `dict.get` returns None for a key that is absent exactly as it does for
one whose value is null, `clean_value` turns both into ".", and "." is indistinguishable from
a legitimately empty annotation. The column that held the right answer all along,
`ONCOKB_vus`, comes from the flattening path, which uses the API's own key names and so cannot
make this mistake.

The test is static: it reads the source, extracts the literal keys, and checks them against a
frozen real response. It therefore needs no token and no network, and it fails on the next
typo rather than on the next release.

Frozen rather than live for the usual reasons -- the endpoint is rate limited and needs a
token, and a dated capture makes a later difference meaningful: if a key disappears, either we
renamed it or OncoKB did. `tests/verification5/data/oncokb_response_frozen.json` records the
endpoint, the variant (BRAF V600E, public) and the data version.
"""
import json
import os
import re

import pytest

HERE = os.path.dirname(__file__)
SRC = os.path.join(HERE, "..", "..", "scripts", "oncokb2.0.py")
FROZEN = os.path.join(HERE, "data", "oncokb_response_frozen.json")


@pytest.fixture(scope="module")
def frozen():
    with open(FROZEN) as fh:
        return json.load(fh)["response"]


@pytest.fixture(scope="module")
def source():
    with open(SRC) as fh:
        return fh.read()


def _code_only(source):
    """The source with comments and docstrings removed.

    Necessary, not fastidious: the comment explaining this bug quotes `data.get("VUS")`
    verbatim, and a regex over the raw file finds it there and reports a defect that is
    only a sentence about the defect. Tokenising is the reliable way to ask what the code
    does rather than what it says.
    """
    import io, tokenize
    out = []
    prev_end, prev_type = (1, 0), tokenize.INDENT
    for tok in tokenize.generate_tokens(io.StringIO(source).readline):
        if tok.type == tokenize.COMMENT:
            continue
        # a string that is the whole statement is a docstring, not a value we care about
        if tok.type == tokenize.STRING and prev_type in (
                tokenize.INDENT, tokenize.NEWLINE, tokenize.NL, tokenize.DEDENT):
            prev_end, prev_type = tok.end, tok.type
            continue
        out.append(tok.string)
        prev_end, prev_type = tok.end, tok.type
    return " ".join(out)


def _keys_read_from(source, obj):
    """The literal keys the code reads off `obj`, e.g. data.get("oncogenic")."""
    return sorted(set(re.findall(rf'\b{obj}\s*\.\s*get\s*\(\s*"([^"]+)"', _code_only(source))))


class TestEveryKeyExists:
    def test_the_frozen_response_is_a_real_one(self, frozen):
        """Guard the guard: a fixture trimmed to nothing would make every test below pass."""
        assert len(frozen) >= 20, f"only {len(frozen)} top-level keys; the capture looks wrong"
        for k in ("oncogenic", "geneExist", "treatments", "mutationEffect", "query"):
            assert k in frozen, f"{k} missing from the frozen response"

    def test_the_api_key_for_vus_is_lower_case(self, frozen):
        """The specific bug, named, so a revert is caught rather than merely detected."""
        assert "vus" in frozen, "the API no longer returns 'vus'"
        assert "VUS" not in frozen, (
            "the API now returns 'VUS' too; check which one oncokb2.0.py should read")

    def test_the_source_reads_vus_lower_case(self, source):
        read = _keys_read_from(source, "data")
        assert "vus" in read, "oncokb2.0.py must read the lower-case 'vus'"
        assert "VUS" not in read, (
            "oncokb2.0.py is reading 'VUS' again; that key does not exist, so ONCOKB_VUS "
            "will be '.' on every row")

    def test_every_key_read_off_the_response_exists(self, source, frozen):
        read = _keys_read_from(source, "data")
        assert read, "found no data.get(...) calls; the extraction pattern changed"
        missing = [k for k in read if k not in frozen]
        assert not missing, (
            f"oncokb2.0.py reads {missing} off the response, and the frozen response has no "
            f"such key(s). A key that is absent returns None, which becomes '.', which looks "
            f"like an empty annotation -- so this is silent. Did the case change? "
            f"Candidates: "
            + ", ".join(f"{m!r}->{next((k for k in frozen if k.lower() == m.lower()), '?')!r}"
                        for m in missing))

    def test_no_key_differs_only_by_case(self, source, frozen):
        """The failure mode that produced this test, stated as its own assertion."""
        lower = {k.lower(): k for k in frozen}
        wrong = [(k, lower[k.lower()]) for k in _keys_read_from(source, "data")
                 if k not in frozen and k.lower() in lower]
        assert not wrong, f"case mismatches between source and API: {wrong}"


class TestNestedKeys:
    def test_mutation_effect_keys_exist(self, source, frozen):
        me = frozen.get("mutationEffect") or {}
        for k in _keys_read_from(source, "me"):
            assert k in me, f"me.get({k!r}) but mutationEffect has {sorted(me)}"

    def test_query_keys_exist(self, source, frozen):
        q = frozen.get("query") or {}
        for k in _keys_read_from(source, "query"):
            assert k in q, f"query.get({k!r}) but query has {sorted(q)}"
