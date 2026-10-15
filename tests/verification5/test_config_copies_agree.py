"""Every Config.yaml copy must describe a shared column the same way as the root.

There are nine Config.yaml files in this repository. They legitimately differ in *which*
columns they list -- each cohort run keeps the config it was run with -- but not in what a
column *means*. That distinction had been lost: corrections landed in the root and nowhere
else, so for months the copies carried

    VAF: 'Frequency of existing variant in 1000 Genomes'

while the root said "computed by SMART as ALT/(REF+ALT) ... not a population frequency".
Eight of the nine copies had the wrong one, including `tests/verification1/Config.yaml` --
the file the test suite loads and the one README's own example command mounts. The
description is printed in the Tier 2/3 metadata row, so a reporting scientist read a 0.33
somatic allele fraction annotated as a population frequency, in the sheet described as
"optimised for clinical scientists".

76 (column, file) pairs had drifted across 23 columns, among them LOEUF's version and
citation, the gnomAD AF versions, and `hotspot`'s and `FractionInformativeReads`'
attribution to the wrong tool -- every one of them a correction that had been made once
and propagated nowhere.

Hand-patching nine files is how it happened. This test is the part that stops it.
"""

import glob
import os
import re

import pytest
import yaml

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
COMPARED = ("description", "source", "version")


def field_entries(node, out=None):
    """Every column entry, wherever it sits in the document."""
    out = {} if out is None else out
    if isinstance(node, dict):
        for key, value in node.items():
            if isinstance(value, dict) and (
                    "description" in value or "source" in value or "tier" in value):
                out[key] = value
            elif isinstance(value, dict):
                field_entries(value, out)
    return out


def load(path):
    with open(path, encoding="utf-8") as fh:
        return field_entries(yaml.safe_load(fh))


def config_copies():
    found = set()
    for pattern in ("*/Config.yaml", "*/*/Config.yaml", "*/*/*/Config.yaml",
                    "*/*/*/*/Config.yaml"):
        for p in glob.glob(os.path.join(ROOT, pattern)):
            if os.path.normpath(p) != os.path.join(ROOT, "Config.yaml"):
                found.add(os.path.normpath(p))
    return sorted(found)


@pytest.fixture(scope="module")
def root_entries():
    return load(os.path.join(ROOT, "Config.yaml"))


class TestTheCopiesAgreeWithTheRoot:
    def test_there_is_more_than_one_config(self):
        """If the copies ever go away this test is free to go with them -- but while they
        exist, nobody should have to guess which one is authoritative."""
        assert config_copies(), "no Config.yaml copies found; delete this test"

    @pytest.mark.parametrize("path", config_copies(),
                             ids=lambda p: os.path.relpath(p, ROOT))
    def test_shared_columns_are_described_identically(self, path, root_entries):
        copy = load(path)
        shared = set(root_entries) & set(copy)
        assert shared, f"{os.path.relpath(path, ROOT)} shares no column with the root"
        drift = []
        for col in sorted(shared):
            for attr in COMPARED:
                if attr not in root_entries[col]:
                    continue
                want = str(root_entries[col].get(attr, "")).strip()
                got = str(copy[col].get(attr, "")).strip()
                if want != got:
                    drift.append(f"{col}.{attr}:\n      root: {want[:150]}\n      here: {got[:150]}")
        assert not drift, (
            f"{len(drift)} entr(ies) in {os.path.relpath(path, ROOT)} disagree with the "
            f"root Config. The root is authoritative; propagate it rather than editing "
            f"one file.\n    " + "\n    ".join(drift))


class TestTheSpecificClaimsThatDrifted:
    """Named individually because each was wrong in a clinically visible way, and a
    generic comparison would not say which."""

    @pytest.mark.parametrize("path", [os.path.join(ROOT, "Config.yaml")] + config_copies(),
                             ids=lambda p: os.path.relpath(p, ROOT))
    def test_vaf_is_not_called_a_population_frequency(self, path):
        entries = load(path)
        if "VAF" not in entries:
            pytest.skip("no VAF column in this config")
        desc = str(entries["VAF"].get("description", ""))
        assert "1000 Genomes" not in desc, (
            "VAF is the tumour allele fraction this pipeline computes, not a population "
            "frequency. This text reaches the Tier 3 clinical sheet.")
        assert "computed by SMART" in desc

    @pytest.mark.parametrize("path", [os.path.join(ROOT, "Config.yaml")] + config_copies(),
                             ids=lambda p: os.path.relpath(p, ROOT))
    def test_loeuf_is_not_called_gene_level(self, path):
        entries = load(path)
        if "LOEUF" not in entries:
            pytest.skip("no LOEUF column in this config")
        desc = str(entries["LOEUF"].get("description", ""))
        assert "gene intolerance" not in desc, (
            "LOEUF is per transcript: 93% of multi-transcript genes carry different "
            "values, which is why the plugin runs with match_by=transcript")

    @pytest.mark.parametrize("path", [os.path.join(ROOT, "Config.yaml")] + config_copies(),
                             ids=lambda p: os.path.relpath(p, ROOT))
    def test_no_entry_claims_a_vep_version_the_image_does_not_install(self, path):
        """Dockerfile pins `ARG VEP_VERSION=114.0` and the annotated VCF header of a real
        run reads `##VEP="v114.0"`. 531 Config entries across the nine files claimed
        114.2, a version this image has never contained -- and the Tier 2/3 metadata row
        printed it next to a banner saying 114.0."""
        dockerfile = open(os.path.join(ROOT, "Dockerfile"), encoding="utf-8").read()
        m = re.search(r"ARG VEP_VERSION=(\S+)", dockerfile)
        assert m, "Dockerfile no longer pins VEP_VERSION; this test needs updating"
        installed = m.group(1)
        entries = load(path)
        wrong = sorted({c for c, e in entries.items()
                        if str(e.get("version", "")).strip().replace("VEP_", "").startswith("114.")
                        and str(e.get("version", "")).strip().replace("VEP_", "") != installed})
        assert not wrong, (
            f"{len(wrong)} entr(ies) claim a VEP version other than the installed "
            f"{installed}: {wrong[:10]}")

    @pytest.mark.parametrize("path", [os.path.join(ROOT, "Config.yaml")] + config_copies(),
                             ids=lambda p: os.path.relpath(p, ROOT))
    def test_no_gnomad_entry_claims_v4_0(self, path):
        """The VEP cache's own info.txt records source_gnomADe v4.1 / source_gnomADg v4.1,
        and the constraint table is verified as v4.1. v4.0 is a different release with a
        different column set."""
        entries = load(path)
        wrong = [c for c, e in entries.items()
                 if str(e.get("version", "")).strip() == "gnomAD_v4.0"]
        assert not wrong, f"{len(wrong)} entr(ies) still claim gnomAD_v4.0: {wrong[:12]}"
