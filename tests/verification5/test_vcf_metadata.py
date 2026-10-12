"""
Unit tests for scripts/vcf_metadata.py

Field metadata is taken from the VCFs rather than from a hand-kept list, because a hand-kept
list drifts: `hotspot` was attributed to Cancer Hotspots when the variant caller had written
it, and `FractionInformativeReads` to Manta when it came from DRAGEN. Every expectation below
is about what the files actually say.

Covers:
  - reading ##INFO / ##FORMAT headers, and keeping the two kinds apart
  - run-level provenance: the caller from ##source, VEP version and cache, the reference
  - who owns a field: the caller, VEP, a VEP custom annotation, a plugin, OncoKB
  - precedence: the VCF wins for what it describes, Config.yaml keeps the rest
  - the checks: name collisions, lost fields, OncoKB versions, liftover
"""
import importlib.util
import os
import re
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "..", "scripts")
_spec = importlib.util.spec_from_file_location("vcf_metadata",
                                               os.path.join(SCRIPTS, "vcf_metadata.py"))
vm = importlib.util.module_from_spec(_spec)
sys.modules["vcf_metadata"] = vm
_spec.loader.exec_module(vm)

DRAGEN = os.path.join(os.path.dirname(__file__), "data", "dragen_like.vcf")

ANNOTATED = '''##fileformat=VCFv4.2
##reference=GRCh38
##VEP="v114.0" API="v114" time="2026-06-30 12:36:27" cache="/refs/homo_sapiens_merged/114_GRCh38"
##SpliceAI_pred_DS_AG=SpliceAI (v1.3) predicted effect on splicing. Delta score for acceptor gain
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##INFO=<ID=hotspot,Number=0,Type=Flag,Description="Known somatic site, used to increase confidence in call">
##INFO=<ID=ClinVar_CLNSIG,Number=.,Type=String,Description="CLNSIG field from /refs/ClinVar/clinvar.vcf.gz">
##INFO=<ID=COSMIC_CNT,Number=.,Type=String,Description="CNT field from /refs/COSMIC/cosmic_cmc_grch38.vcf.gz">
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|AF|SpliceAI_pred_DS_AG|REVEL|ClinVar_CLNSIG|COSMIC_CNT">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t.\tPASS\tDP=8;hotspot
'''

ONCOKB = '''##fileformat=VCFv4.2
##INFO=<ID=ClinVar_CLNSIG,Number=.,Type=String,Description="CLNSIG field from /refs/ClinVar/clinvar.vcf.gz">
##INFO=<ID=ONCOKB_ONCOGENIC,Number=1,Type=String,Description="OncoKB oncogenic classification">
##INFO=<ID=ONCOKB_DATA_VERSION,Number=1,Type=String,Description="OncoKB data version">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t.\tPASS\tONCOKB_ONCOGENIC=Oncogenic;ONCOKB_DATA_VERSION=v7.6
'''


@pytest.fixture
def vcfs(tmp_path):
    a = tmp_path / "s_annotated.vcf"; a.write_text(ANNOTATED)
    o = tmp_path / "s.oncoKB.vcf"; o.write_text(ONCOKB)
    return str(a), str(o)


@pytest.fixture
def meta(vcfs):
    a, o = vcfs
    return vm.harvest(annotated=a, oncokb=o, input_vcf=DRAGEN, oncokb_version="v7.6")


# ===========================================================================
# Reading the headers
# ===========================================================================

class TestReadHeaders:
    def test_info_and_format_are_read_with_their_description(self):
        fields, _, _ = vm.read_headers(DRAGEN)
        assert fields["hotspot"]["kind"] == "INFO"
        assert fields["hotspot"]["type"] == "Flag"
        assert fields["hotspot"]["description"] == \
            "Known somatic site, used to increase confidence in call"
        assert fields["AF"]["kind"] == "FORMAT"
        assert fields["AF"]["description"] == "Allele fraction in tumor"

    def test_run_level_provenance(self, vcfs):
        _, run, _ = vm.read_headers(DRAGEN)
        assert run["source"] == "DRAGEN_SOMATIC 4.5.1"
        assert run["reference"] == "file:///ref/GRCh38.fa"
        _, run_a, csq = vm.read_headers(vcfs[0])
        assert run_a["vep_version"] == "v114.0"
        assert run_a["vep_cache"].endswith("114_GRCh38")
        assert csq[:4] == ["Allele", "Consequence", "IMPACT", "SYMBOL"]

    def test_missing_file_is_not_an_error(self):
        assert vm.read_headers("/nowhere/at/all.vcf") == ({}, {}, [])


# ===========================================================================
# Who owns a field
# ===========================================================================

class TestOwnership:
    def test_caller_fields_are_attributed_to_the_input_vcf(self, meta):
        """The defect this module exists to prevent: these two were attributed to
        Cancer Hotspots and to Manta, when DRAGEN wrote both."""
        for name in ("hotspot", "FractionInformativeReads", "SoftClipRatio", "STR", "RU", "RPA"):
            d = meta["fields"][name]
            assert d["origin"] == "input"
            assert "DRAGEN_SOMATIC 4.5.1" in d["source"]
        assert meta["fields"]["hotspot"]["description"].startswith("Known somatic site")

    def test_a_custom_annotation_names_its_resource(self, meta):
        """VEP names the fields of a custom annotation after the annotation, and its
        Description names the file — which is the version actually used in this run."""
        assert meta["fields"]["ClinVar_CLNSIG"]["source"] == "ClinVar"
        assert meta["fields"]["ClinVar_CLNSIG"]["version"] == "clinvar.vcf.gz"
        assert meta["fields"]["COSMIC_CNT"]["source"] == "COSMIC"
        assert meta["fields"]["COSMIC_CNT"]["version"] == "cosmic_cmc_grch38.vcf.gz"

    def test_a_plugin_announces_its_own_version(self, meta):
        d = meta["fields"]["SpliceAI_pred_DS_AG"]
        assert (d["source"], d["version"]) == ("SpliceAI", "v1.3")

    def test_core_vep_fields_carry_the_vep_version(self, meta):
        d = meta["fields"]["Consequence"]
        assert (d["source"], d["version"], d["origin"]) == ("VEP", "v114.0", "vep")

    def test_only_fields_oncokb_adds_are_attributed_to_oncokb(self, meta):
        """The OncoKB VCF is the annotated VCF with OncoKB's fields added, so everything it
        inherits keeps the provenance it already had."""
        assert meta["fields"]["ONCOKB_ONCOGENIC"]["source"] == "OncoKB"
        assert meta["fields"]["ONCOKB_ONCOGENIC"]["version"] == "v7.6"
        assert meta["fields"]["ClinVar_CLNSIG"]["source"] == "ClinVar"   # not OncoKB

    def test_a_format_field_does_not_claim_the_column_name(self, meta):
        """post_analysis expands FORMAT and renames it (AF becomes VAF), so the column called
        AF holds VEP's population frequency. Letting the caller's FORMAT AF claim the name
        would label a population frequency as a tumour allele fraction."""
        assert meta["fields"]["AF"]["source"] == "VEP"
        assert "AF" in meta["input_format_fields"]
        assert meta["input_format_fields"]["AF"] == "Allele fraction in tumor"


# ===========================================================================
# Precedence
# ===========================================================================

class TestPrecedence:
    CFG = {
        "hotspot": {"description": "Binary flag indicating presence in any Cancer Hotspots "
                                   "annotation", "source": "Cancer Hotspots",
                    "version": "changv2", "tier": 3},
        "ClinVar_CLNSIG": {"description": "Clinical significance for this variant",
                           "source": "ClinVar", "version": "ClinVar_2024-12", "tier": 3},
        "Hugo_Symbol": {"description": "Gene symbol", "source": "VEP", "version": "114", "tier": 3},
        "ONCOKB_TX_LEVEL_1": {"description": "Level 1 therapy", "source": "OncoKB",
                              "version": "stale", "tier": 3},
    }

    def test_the_vcf_wins_for_a_field_it_describes(self, meta):
        r = vm.resolve("hotspot", meta, self.CFG)
        assert "DRAGEN" in r["source"]            # not Cancer Hotspots
        assert r["version"] != "changv2"
        assert r["description"].startswith("Known somatic site")
        assert r["tier"] == 3                      # the tier still comes from Config.yaml

    def test_config_keeps_fields_no_vcf_describes(self, meta):
        r = vm.resolve("Hugo_Symbol", meta, self.CFG)
        assert (r["description"], r["source"]) == ("Gene symbol", "VEP")

    def test_mafannotator_columns_take_the_run_oncokb_version(self, meta):
        """They have no VCF header, but the analysis cannot span two OncoKB releases, so the
        version is the one observed in this run."""
        r = vm.resolve("ONCOKB_TX_LEVEL_1", meta, self.CFG, oncokb_version="v7.6")
        assert r["version"] == "v7.6"
        assert r["description"] == "Level 1 therapy"

    def test_a_pattern_in_config_still_applies(self, meta):
        patterns = [{"pattern": "ONCOKB_TX_*", "description": "from a pattern",
                     "source": "OncoKB", "tier": 2}]
        r = vm.resolve("ONCOKB_TX_LEVEL_9", meta, {}, patterns, oncokb_version="v7.6")
        assert r["description"] == "from a pattern" and r["version"] == "v7.6"

    def test_conflicts_are_reported_on_the_source(self, meta):
        found = {name: (claimed, actual) for name, claimed, actual
                 in vm.conflicts(meta, self.CFG)}
        assert "hotspot" in found
        assert found["hotspot"][0] == "Cancer Hotspots"
        assert "DRAGEN" in found["hotspot"][1]
        # the same source worded the same way is not a conflict
        assert "ClinVar_CLNSIG" not in found


# ===========================================================================
# The checks
# ===========================================================================

class TestChecks:
    def level(self, results, text):
        return [lv for lv, msg in results if text in msg]

    def test_an_info_collision_is_an_error(self, meta):
        """An INFO entry is applied after the CSQ fields, so a shared name replaces the
        annotation's value — a wrong value, not a missing one."""
        meta = dict(meta, input_info_ids=["REVEL"], csq_fields=["REVEL", "Consequence"])
        assert self.level(vm.check(meta), "name collision") == ["error"]

    def test_a_format_collision_is_not_a_collision(self, meta):
        """AF is FORMAT in the input and a CSQ field in VEP, and the pipeline renames one."""
        assert not self.level(vm.check(meta), "name collision")

    def test_lost_fields_are_reported(self, meta):
        res = vm.check(meta, table_columns=["DP", "MQ"])
        lost = [msg for lv, msg in res if "did not" in msg]
        assert lost and "hotspot" in lost[0]

    def test_nothing_lost_when_every_field_arrives(self, meta):
        cols = meta["input_info_ids"]
        assert not [m for lv, m in vm.check(meta, table_columns=cols) if "did not" in m]

    def test_two_oncokb_versions_in_one_run_is_an_error(self, meta):
        res = vm.check(meta, oncokb_versions=["v7.5", "v7.6", ""])
        assert self.level(res, "more than one OncoKB data version") == ["error"]
        assert not self.level(vm.check(meta, oncokb_versions=["v7.6", "v7.6", ""]),
                              "more than one OncoKB")

    def test_liftover_warns_about_caller_fields(self, meta):
        assert self.level(vm.check(meta, liftover=True), "liftover enabled") == ["warning"]
        assert not self.level(vm.check(meta, liftover=False), "liftover enabled")

    def test_a_missing_vep_header_is_reported(self):
        bare = vm.harvest(input_vcf=DRAGEN)
        assert self.level(vm.check(bare), "VEP version unknown") == ["warning"]


# ===========================================================================
# The sidecar
# ===========================================================================

class TestSidecar:
    def test_round_trip_through_the_cli(self, vcfs, tmp_path, capsys):
        a, o = vcfs
        out = tmp_path / "meta.json"
        assert vm.main(["--annotated", a, "--oncokb", o, "--input", DRAGEN,
                        "--oncokb-version", "v7.6", "--out", str(out)]) == 0
        assert out.exists()
        assert vm.main(["--sidecar", str(out), "--summary"]) == 0
        printed = capsys.readouterr().out
        assert "fields described by the VCFs" in printed
        assert "DRAGEN_SOMATIC 4.5.1" in printed


# ===========================================================================
# Against a real DRAGEN VCF
#   tests/verification5/data/dragen_real_subset.vcf — a de-identified extract of a run from
#   the laboratory: five somatic calls, the headers intact. Every expectation here is read
#   from that file, which is the point: when the same fields were written by hand into a
#   synthetic fixture, all eight descriptions came out subtly wrong.
# ===========================================================================

REAL = os.path.join(os.path.dirname(__file__), "data", "dragen_real_subset.vcf")


class TestRealDragenVcf:
    @pytest.fixture
    def meta(self):
        return vm.harvest(input_vcf=REAL)

    def test_the_caller_and_its_version_come_from_the_file(self, meta):
        assert meta["run"]["input_caller"] == "DRAGEN_SNV"
        assert meta["run"]["reference"].startswith("file:///")

    def test_every_info_field_dragen_writes_is_described(self, meta):
        """Six of these had no entry in Config.yaml at all, and two were attributed to the
        wrong source. The file describes all of them."""
        for name in ("DP", "MQ", "MQRankSum", "ReadPosRankSum", "FractionInformativeReads",
                     "SoftClipRatio", "STR", "RU", "RPA", "hotspot"):
            d = meta["fields"][name]
            assert d["origin"] == "input", name
            assert "DRAGEN_SNV" in d["source"], name
            assert d["description"], f"{name} has no description"

    def test_the_hotspot_flag_is_the_callers_own(self, meta):
        """The defect that started this: Config.yaml attributed it to Cancer Hotspots."""
        d = meta["fields"]["hotspot"]
        assert d["description"] == "Known somatic site, used to increase confidence in call"
        assert "DRAGEN" in d["source"] and "Cancer Hotspots" not in d["source"]

    def test_format_fields_are_kept_apart_from_info(self, meta):
        """AF is FORMAT here — a tumour allele fraction — while the AF column of the table is
        VEP's population frequency. Keeping them apart is what stops one being labelled as
        the other."""
        assert "AF" not in meta["input_info_ids"]
        assert meta["input_format_fields"]["AF"] == \
            "Allele fractions for alt alleles in the order listed"
        for name in ("GT", "AD", "SQ", "F1R2", "F2R1", "SB", "MB", "PS"):
            assert name in meta["input_format_fields"], name

    def test_no_collision_with_the_annotation(self, meta):
        """None of DRAGEN's INFO names clashes with a VEP field, so nothing is overwritten."""
        meta = dict(meta, csq_fields=["Consequence", "IMPACT", "SYMBOL", "AF", "REVEL"])
        assert not [m for lv, m in vm.check(meta) if "collision" in m]

    def test_the_extract_carries_no_identifying_content(self):
        """It is kept in the repository, so the de-identification is part of the contract."""
        text = open(REAL).read()
        for pattern in (r"Date=", r"referenceSexKaryotype", r"bham",
                        r"\d{2}/\d{2}/\d{4}"):
            assert not re.search(pattern, text), f"{pattern} is still in the fixture"
        # a run of heterozygous germline genotypes is a genetic fingerprint; these are somatic
        fractions = [float(l.split("\t")[9].split(":")[3])
                     for l in text.splitlines() if not l.startswith("#")]
        assert fractions and max(fractions) < 0.4
