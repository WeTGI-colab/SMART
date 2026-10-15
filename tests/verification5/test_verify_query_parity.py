"""The API verifier must query OncoKB the way the pipeline does.

verification1 compares the pipeline's OncoKB columns against a fresh API call. That is
only a valid comparison if both sides ask the same question. They did not: the verifier
sent `byProteinChange` with the alteration alone, while MafAnnotator also sends the
consequence, read from `CONSEQUENCE_HEADERS = ['VARIANT_CLASSIFICATION', 'MUTATION_TYPE']`.

Measured against the live API for TP53 `T125=`:

    alteration=T125=                                  -> Unknown / Unknown
    alteration=T125=&consequence=Splice_Region        -> Likely Oncogenic / Likely Loss-of-function
    alteration=T125=&consequence=splice_region_variant-> Likely Oncogenic / Likely Loss-of-function
    alteration=T125 (no suffix)                       -> Likely Oncogenic / Unknown

The pipeline's output held the second of those, which is correct. The verifier reported
it as a MISMATCH for two fields on every run, because its own query was weaker. A check
that fails correct output costs more than one that passes broken output: it teaches the
reader to ignore the red lines that are always there.

No network: these tests assert the request that is built, not the answer that comes back.
"""

import importlib.util
import os
import sys

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_VERIFY = os.path.join(_HERE, "..", "verification1", "verify.py")

pytestmark = pytest.mark.skipif(not os.path.exists(_VERIFY),
                                reason="verification1/verify.py not present")


def load_verify(monkeypatch):
    """Import verify.py with its network entry point replaced by a recorder."""
    spec = importlib.util.spec_from_file_location("v1_verify", _VERIFY)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["v1_verify"] = mod
    spec.loader.exec_module(mod)
    calls = []

    class Resp:
        status_code = 200
        def raise_for_status(self):
            return None
        def json(self):
            return {"oncogenic": "Likely Oncogenic"}

    def recorder(url, **kw):
        calls.append((url, kw))
        return Resp()

    monkeypatch.setattr(mod, "http_get", recorder)
    return mod, calls


class TestTheConsequenceIsSent:
    def test_the_consequence_reaches_the_request(self, monkeypatch):
        mod, calls = load_verify(monkeypatch)
        mod.query_oncokb_mutation("TP53", "T125=", "tok", "Splice_Region")
        assert len(calls) == 1
        params = calls[0][1]["params"]
        assert params["consequence"] == "Splice_Region"
        assert params["alteration"] == "T125="
        assert params["hugoSymbol"] == "TP53"

    def test_an_empty_consequence_is_omitted_rather_than_sent_blank(self, monkeypatch):
        """A blank consequence is not the same request as no consequence, and OncoKB
        need not treat them alike."""
        mod, calls = load_verify(monkeypatch)
        mod.query_oncokb_mutation("TP53", "R175H", "tok", "")
        assert "consequence" not in calls[0][1]["params"]

    def test_it_is_optional_so_older_callers_still_work(self, monkeypatch):
        mod, calls = load_verify(monkeypatch)
        mod.query_oncokb_mutation("TP53", "R175H", "tok")
        assert "consequence" not in calls[0][1]["params"]

    def test_the_caller_passes_variant_classification_not_consequence(self):
        """MafAnnotator reads VARIANT_CLASSIFICATION. Sending VEP's `Consequence`
        instead would be a different string -- and for this variant a worse one, since
        its Consequence is the compound `splice_region_variant&synonymous_variant` whose
        synonymous half returns Unknown."""
        src = open(_VERIFY).read()
        assert 'row.get("Variant_Classification"' in src, (
            "the verifier must take the consequence from Variant_Classification, the "
            "header MafAnnotator reads")
        i = src.index("api_data = query_oncokb_mutation(")
        window = src[max(0, i - 500):i]
        assert "Variant_Classification" in window, (
            "the consequence must be derived at the call site, not left unset")


class TestTheParityIsDocumented:
    def test_the_function_explains_why_the_consequence_matters(self):
        src = open(_VERIFY).read()
        i = src.index("def query_oncokb_mutation")
        doc = src[i:i + 1400]
        assert "T125" in doc, "the measurement that proved it should be recorded"
        assert "CONSEQUENCE_HEADERS" in doc, (
            "say where MafAnnotator gets the consequence, so the parity can be rechecked")
