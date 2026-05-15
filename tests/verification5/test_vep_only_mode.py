"""
Tests for the no-token → VEP-only behaviour added in entrypoint.sh.

Two levels of testing:

1. Bash-level  — executes the exact token-detection snippet from entrypoint.sh
   to verify that missing/empty token correctly sets VEP_ONLY=1, without
   needing Docker or any reference files.

2. Python-level — verifies that the OncoKB query functions return None when
   called with an empty token (HTTP 401), documenting why the pipeline must
   skip OncoKB when no token is available.

Neither group requires Docker, an OncoKB token, or reference files.
"""
import os
import subprocess
from unittest.mock import MagicMock, patch

import pytest

from conftest import oncokb

ENTRYPOINT = os.path.join(
    os.path.dirname(__file__), "..", "..", "entrypoint.sh"
)

# ---------------------------------------------------------------------------
# The bash snippet that mirrors the token-detection block in entrypoint.sh.
# We test it in isolation so this suite does not need Docker.
# ---------------------------------------------------------------------------
_TOKEN_DETECT = r"""
set -euo pipefail
VEP_ONLY=0
ONCOKB_TOKEN=""
if [[ -z "${1:-}" || "${1:-}" == --* ]]; then
    ONCOKB_TOKEN=""
    VEP_ONLY=1
else
    ONCOKB_TOKEN="$1"
fi
echo "VEP_ONLY=${VEP_ONLY}"
echo "ONCOKB_TOKEN=${ONCOKB_TOKEN}"
"""


def _detect(first_arg: str) -> dict:
    """Run the token-detection snippet with first_arg and return parsed output."""
    result = subprocess.run(
        ["bash", "-c", _TOKEN_DETECT, "--", first_arg],
        capture_output=True, text=True, check=True,
    )
    output = {}
    for line in result.stdout.strip().splitlines():
        k, _, v = line.partition("=")
        output[k] = v
    return output


# ===========================================================================
# Bash-level: token detection logic
# ===========================================================================

class TestTokenDetection:
    """Verify the bash snippet that auto-sets VEP_ONLY when no token given."""

    def test_empty_string_triggers_vep_only(self):
        out = _detect("")
        assert out["VEP_ONLY"] == "1"
        assert out["ONCOKB_TOKEN"] == ""

    def test_flag_as_first_arg_triggers_vep_only(self):
        # User runs: smart --ref-dir /refs  (forgot the token)
        out = _detect("--ref-dir")
        assert out["VEP_ONLY"] == "1"

    def test_any_double_dash_flag_triggers_vep_only(self):
        for flag in ["--no-liftover", "--vep-only", "--keep-tmp", "--help"]:
            out = _detect(flag)
            assert out["VEP_ONLY"] == "1", f"Expected VEP_ONLY=1 for flag {flag!r}"

    def test_valid_token_does_not_trigger_vep_only(self):
        out = _detect("mytoken_abc123")
        assert out["VEP_ONLY"] == "0"
        assert out["ONCOKB_TOKEN"] == "mytoken_abc123"

    def test_token_is_preserved_when_given(self):
        token = "eyJhbGciOiJSUzI1NiJ9.test"
        out = _detect(token)
        assert out["ONCOKB_TOKEN"] == token
        assert out["VEP_ONLY"] == "0"


# ===========================================================================
# Bash-level: help text describes the no-token behaviour
# ===========================================================================

class TestEntrypointHelp:
    """Verify the usage text correctly documents the optional-token behaviour."""

    @pytest.fixture(scope="class")
    def help_text(self):
        result = subprocess.run(
            ["bash", ENTRYPOINT, "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, "entrypoint.sh --help should exit 0"
        return result.stdout.lower()

    def test_help_exits_cleanly(self):
        result = subprocess.run(
            ["bash", ENTRYPOINT, "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0

    def test_help_says_token_is_omittable(self, help_text):
        assert "omitted" in help_text or "optional" in help_text

    def test_help_says_vep_only_is_auto_triggered(self, help_text):
        assert "automatically" in help_text

    def test_help_mentions_vep_only_flag(self, help_text):
        assert "--vep-only" in help_text


# ===========================================================================
# Python-level: OncoKB functions fail gracefully with empty token
# ===========================================================================

class TestOncokbRequiresToken:
    """
    Document that calling OncoKB query functions with an empty token returns
    None (HTTP 401), confirming why the pipeline must skip OncoKB in
    VEP-only mode rather than calling with an empty token.
    """

    def _mock_response(self, status_code: int):
        resp = MagicMock()
        resp.status_code = status_code
        resp.json.return_value = {"error": "Unauthorized"}
        return resp

    def test_mutation_query_returns_none_on_401(self):
        with patch("requests.get", return_value=self._mock_response(401)):
            result = oncokb.query_oncokb_mutation(
                gene="BRAF", alteration="V600E",
                tumor_type="", token="", cache={},
            )
        assert result is None

    def test_cna_query_returns_none_on_401(self):
        with patch("requests.get", return_value=self._mock_response(401)):
            result = oncokb.query_oncokb_cna(
                gene="CDKN2A", cna_type="DELETION",
                tumor_type="", token="", cache={},
            )
        assert result is None

    def test_mutation_query_returns_none_for_empty_gene(self):
        # Guard: empty gene/alteration returns None before any HTTP call
        result = oncokb.query_oncokb_mutation(
            gene="", alteration="V600E",
            tumor_type="", token="sometoken", cache={},
        )
        assert result is None

    def test_mutation_query_uses_cache_on_second_call(self):
        """Second call with same args should NOT hit the network."""
        cache = {}
        cached_value = {"oncogenic": "Oncogenic"}
        key = ("mut", "BRAF", "V600E", "")
        cache[key] = cached_value

        with patch("requests.get") as mock_get:
            result = oncokb.query_oncokb_mutation(
                gene="BRAF", alteration="V600E",
                tumor_type="", token="", cache=cache,
            )
            mock_get.assert_not_called()

        assert result == cached_value

    def test_mutation_query_successful_with_valid_response(self):
        """With a real token and 200 response, data is returned and cached."""
        expected = {"oncogenic": "Oncogenic", "hotspot": True}
        resp = MagicMock()
        resp.status_code = 200
        resp.json.return_value = expected

        cache = {}
        with patch("requests.get", return_value=resp):
            result = oncokb.query_oncokb_mutation(
                gene="BRAF", alteration="V600E",
                tumor_type="", token="valid_token", cache=cache,
            )
        assert result == expected
        # Result is cached for next call
        assert ("mut", "BRAF", "V600E", "") in cache
