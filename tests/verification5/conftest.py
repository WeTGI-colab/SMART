"""
Shared fixtures and import helpers for verification5 unit tests.

Scripts live in scripts/ and one of them (oncokb2.0.py) is not a valid
Python identifier, so all three are loaded via importlib rather than a
normal import statement.
"""
import importlib.util
import os
import sys
from unittest.mock import MagicMock

import pytest

# cyvcf2 requires compiled C extensions unavailable in plain Python environments.
# Mock it so the scripts can be imported for unit-testing their pure-Python logic.
if "cyvcf2" not in sys.modules:
    sys.modules["cyvcf2"] = MagicMock()

SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "scripts")


def _load(module_name: str, filename: str):
    path = os.path.join(SCRIPTS_DIR, filename)
    spec = importlib.util.spec_from_file_location(module_name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = mod
    spec.loader.exec_module(mod)
    return mod


# Expose the three scripts as module-level names so test files can import them.
vcf2table = _load("vcf2table", "vcf2table.py")
oncokb = _load("oncokb2", "oncokb2.0.py")
post = _load("post_analysis", "post_analysis.py")


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def nm_file(tmp_path):
    """Transcript whitelist file with version suffixes."""
    p = tmp_path / "transcripts.txt"
    p.write_text(
        "NM_000077.5\n"
        "NM_002524.3\n"
        "NM_005228.5\n"
        "# comment line\n"
        "\n"
        "NM_000546.6\n"
    )
    return str(p)


@pytest.fixture
def minimal_yaml(tmp_path):
    """Minimal Config.yaml for testing post_analysis helpers."""
    content = """\
fields:
  SYMBOL:
    description: Gene symbol
    source: VEP
    version: "114.2"
    tier: 3
  ONCOKB_ONCOGENIC:
    description: Oncogenicity classification
    source: OncoKB
    version: "5.4"
    tier: 3
  UNKNOWN_FIELD:
    description: Internal placeholder
    source: Internal
    version: "1.0"
    tier: drop
field_patterns:
  - pattern: "ONCOKB_TX_*"
    description: OncoKB treatment implication entry
    source: OncoKB
    version: "5.4"
    tier: 2
"""
    p = tmp_path / "Config.yaml"
    p.write_text(content)
    return str(p)


@pytest.fixture
def vep_vcf(tmp_path):
    """Minimal VEP-annotated VCF with a valid CSQ header."""
    header = (
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. '
        'Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|HGVSc|HGVSp|'
        'MANE_SELECT|MANE_PLUS_CLINICAL">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    p = tmp_path / "test.vcf"
    p.write_text(header)
    return str(p)


@pytest.fixture
def vep_vcf_no_csq(tmp_path):
    """VCF without a CSQ header — should cause parse_csq_format to raise."""
    p = tmp_path / "no_csq.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    return str(p)
