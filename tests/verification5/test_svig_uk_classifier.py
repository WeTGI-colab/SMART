"""
Unit tests for scripts/svig_uk_classifier.py

Every expectation is taken from the ACGS / SVIG-UK 2025 guideline or its
Supplementary Material — never from observed classifier output — so these tests
document the rules the scorer is required to follow.

Covers:
  - normalise_protein()        — HGVS protein normalisation used for all lookups
  - variant_class()            — VEP consequence -> SVIG-UK variant category
  - stand-alone codes          — O1, B1 (>1%), B2 (null variant in an oncogene)
  - O2 PVS1 proxy              — NMD / final exon / first 100 bp / +2T>C
  - O3 / B1 population rules   — absent, <=0.001%, >1 allele, >0.1%
  - O6 / B3 / B4 computational — REVEL and SpliceAI thresholds, LOF restriction
  - O7 hotspot tiers, O8 missense constraint, O4 case counting, O9 / B5
  - combination rules          — Supplementary Figures 1 and 2
  - minimum-two-codes rule     — guideline §6.3
"""
import importlib.util
import json
import os
import sys

import pytest

SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "scripts")
_spec = importlib.util.spec_from_file_location(
    "svig_uk_classifier", os.path.join(SCRIPTS_DIR, "svig_uk_classifier.py"))
svig = importlib.util.module_from_spec(_spec)
sys.modules["svig_uk_classifier"] = svig
_spec.loader.exec_module(svig)

CANONICAL = os.path.join(SCRIPTS_DIR, "..", "reference_lists", "svig_uk_canonical_variants.tsv")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_resources(roles=None, hotspots=None, mis_z=None, cosmic_rows=None, tmp_path=None):
    res = svig.Resources()
    res.load_canonical(CANONICAL)
    res.gene_roles = roles or {}
    res.gene_roles_source = "test"
    for e in hotspots or []:
        res.hotspots_single[(e["hugoSymbol"], e["aminoAcidPosition"]["start"])] = e
    res.mis_z = mis_z or {}
    if cosmic_rows is not None:
        path = tmp_path / "cosmic.tsv"
        header = "gene\ttranscript\tprotein\tmut_class\taa_start\taa_end\tn_samples\n"
        path.write_text(header + "".join("\t".join(map(str, r)) + "\n" for r in cosmic_rows))
        res.cosmic = svig.CosmicO4(str(path))
    return res


def row(**kw):
    base = {
        "Hugo_Symbol": "GENE1", "Chromosome": "chr1", "Start_Position": "100",
        "Reference_Allele": "C", "Tumor_Seq_Allele2": "T",
        "Consequence": "missense_variant", "HGVSp": "ENSP1:p.Arg100His", "HGVSc": "c.299G>A",
        "EXON": "3/10", "CDS_position": "299", "DOMAINS": "",
        "REVEL": "", "SpliceAI_pred_DS_AG": "", "SpliceAI_pred_DS_AL": "",
        "SpliceAI_pred_DS_DG": "", "SpliceAI_pred_DS_DL": "",
        "gnomADe_AF": "", "gnomADg_AF": "", "VARIANT_IN_ONCOKB": "False", "MUTATION_EFFECT": "",
    }
    base.update(kw)
    return base


def spliceai(v):
    return {k: str(v) for k in ("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
                                "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")}


def codes(result):
    return {e.code: e for e in result.codes}


TSG = {"GENE1": {"TSG"}}
ONCOGENE = {"GENE1": {"ONC"}}


# ===========================================================================
# normalise_protein / variant_class
# ===========================================================================

class TestNormalisation:
    @pytest.mark.parametrize("hgvsp,expected", [
        ("ENSP00000358548.4:p.Val600Glu", "V600E"),
        ("p.Trp288CysfsTer12", "W288Cfs"),
        ("p.W288Cfs*12", "W288Cfs"),
        ("p.T15GlufsTer", "T15Efs"),           # malformed legacy SMART HGVSp_Short
        ("p.Pro95_Arg102del", "P95_R102del"),
        ("p.(Arg201Ter)", "R201*"),
        ("p.Thr477%3D", "T477="),
        ("p.?", ""),
        ("", ""),
    ])
    def test_normalise_protein(self, hgvsp, expected):
        assert svig.normalise_protein(hgvsp) == expected

    @pytest.mark.parametrize("consequence,expected", [
        ("stop_gained&frameshift_variant", "truncating"),
        ("frameshift_variant&splice_region_variant", "truncating"),
        ("splice_donor_variant&intron_variant", "canonical_splice"),
        ("missense_variant&splice_region_variant", "missense"),
        ("inframe_deletion", "inframe"),
        ("synonymous_variant", "synonymous"),
        ("3_prime_UTR_variant", "utr"),
        ("intron_variant", "intronic"),
    ])
    def test_variant_class(self, consequence, expected):
        assert svig.variant_class(consequence) == expected


# ===========================================================================
# Stand-alone codes
# ===========================================================================

class TestStandalone:
    def test_o1_canonical_protein_match(self):
        res = make_resources()
        r = svig.classify_row(row(Hugo_Symbol="DNMT3A", HGVSp="ENSP1:p.Arg882His"), res)
        assert r.classification == "Oncogenic"
        assert r.codes[0].code == "O1"

    def test_o1_carries_likely_oncogenic_assessment(self):
        # Supplementary Table 3 lists PDGFRA D842V as Likely Oncogenic
        res = make_resources()
        r = svig.classify_row(row(Hugo_Symbol="PDGFRA", HGVSp="p.Asp842Val"), res)
        assert r.classification == "Likely Oncogenic"

    def test_o1_frameshift_matches_regardless_of_tail(self):
        res = make_resources()
        r = svig.classify_row(row(Hugo_Symbol="NPM1", Consequence="frameshift_variant",
                                  HGVSp="p.Trp288CysfsTer12"), res)
        assert r.classification == "Oncogenic"

    def test_o1_tert_promoter_genomic_match(self):
        res = make_resources()
        r = svig.classify_row(row(Hugo_Symbol="TERT", Chromosome="chr5", Start_Position="1295113",
                                  Reference_Allele="G", Tumor_Seq_Allele2="A",
                                  Consequence="upstream_gene_variant", HGVSp=""), res)
        assert r.classification == "Oncogenic"

    def test_b1_standalone_above_1_percent(self):
        r = svig.classify_row(row(gnomADe_NFE_AF="0.02"), make_resources(roles=TSG))
        assert r.classification == "Benign"
        assert r.codes[0].code == "B1"

    def test_b2_null_variant_in_oncogene_is_standalone_vus(self):
        r = svig.classify_row(row(Consequence="stop_gained", HGVSp="p.Arg100Ter"),
                              make_resources(roles=ONCOGENE))
        assert r.classification == "VUS"
        assert r.codes[0].code == "B2"

    def test_b2_not_applied_to_final_exon_truncation_in_oncogene(self):
        # Final-exon truncation in an oncogene is the O9 gain-of-function route [S T1 O9]
        r = svig.classify_row(row(Consequence="frameshift_variant", HGVSp="p.Arg100fs", EXON="10/10"),
                              make_resources(roles=ONCOGENE))
        assert "B2" not in codes(r)
        assert codes(r)["O9"].points == 2

    def test_ezh2_mode_of_action_depends_on_tumour_type(self):
        # [G §7.2] EZH2: TSG in myeloid, activating in lymphoid
        res = make_resources(roles={"EZH2": {"ONC", "TSG"}})
        lof = row(Hugo_Symbol="EZH2", Consequence="stop_gained", HGVSp="p.Arg100Ter")
        assert "O2" in codes(svig.classify_row(dict(lof, TUMOUR_TYPE="Myeloid"), res))
        assert svig.classify_row(dict(lof, TUMOUR_TYPE="Lymphoid"), res).codes[0].code == "B2"


# ===========================================================================
# O2 — PVS1 proxy
# ===========================================================================

class TestO2:
    def test_nmd_predicted_is_very_strong(self):
        r = svig.classify_row(row(Consequence="stop_gained", HGVSp="p.Arg100Ter", EXON="3/10"),
                              make_resources(roles=TSG))
        assert codes(r)["O2"].points == 8

    def test_final_exon_is_moderate(self):
        r = svig.classify_row(row(Consequence="stop_gained", HGVSp="p.Arg100Ter", EXON="10/10"),
                              make_resources(roles=TSG))
        assert codes(r)["O2"].points == 2

    def test_first_100bp_of_exon1_is_moderate(self):
        r = svig.classify_row(row(Consequence="stop_gained", HGVSp="p.Arg20Ter", EXON="1/10",
                                  CDS_position="58"), make_resources(roles=TSG))
        assert codes(r)["O2"].points == 2

    def test_plus2_t_to_c_needs_spliceai_0_8(self):
        base = row(Consequence="splice_donor_variant", HGVSp="", HGVSc="c.300+2T>C", EXON="")
        res = make_resources(roles=TSG)
        assert "O2" not in codes(svig.classify_row(dict(base, **spliceai(0.5)), res))
        assert "O2" in codes(svig.classify_row(dict(base, **spliceai(0.9)), res))

    def test_o2_not_applied_when_mode_of_action_unknown(self):
        r = svig.classify_row(row(Consequence="stop_gained", HGVSp="p.Arg100Ter"), make_resources())
        assert "O2" not in codes(r)


# ===========================================================================
# O3 / B1 — population frequency
# ===========================================================================

class TestPopulation:
    def test_absent_is_o3_moderate(self):
        assert codes(svig.classify_row(row(), make_resources()))["O3"].points == 2

    def test_rare_in_large_subpopulation_is_supporting(self):
        # 2 alleles in >1.1M NFE exome alleles -> AF ~1.8e-6 <= 0.001%
        r = svig.classify_row(row(gnomADe_AF="1.4e-06", gnomADe_NFE_AF="1.8e-06"), make_resources())
        assert codes(r)["O3"].points == 1

    def test_more_than_one_allele_in_small_subpopulation_blocks_o3(self):
        # 3 alleles in ~33,600 AFR exome alleles
        r = svig.classify_row(row(gnomADe_AF="2e-06", gnomADe_AFR_AF="8.9e-05"), make_resources())
        assert "O3" not in codes(r)

    def test_b1_strong_between_0_1_and_1_percent(self):
        r = svig.classify_row(row(gnomADg_AFR_AF="0.005"), make_resources())
        assert codes(r)["B1"].points == -4


# ===========================================================================
# O6 / B3 / B4 — computational evidence
# ===========================================================================

class TestComputational:
    def test_revel_0_7_gives_o6(self):
        r = svig.classify_row(row(REVEL="0.72"), make_resources(roles=TSG))
        assert codes(r)["O6"].points == 1

    def test_revel_benign_threshold_is_0_4_not_0_7(self):
        res = make_resources(roles=TSG)
        assert "B3" not in codes(svig.classify_row(row(REVEL="0.55"), res))
        assert codes(svig.classify_row(row(REVEL="0.30"), res))["B3"].points == -1

    def test_revel_not_used_for_benignity_in_activating_genes(self):
        r = svig.classify_row(row(REVEL="0.10"), make_resources(roles=ONCOGENE))
        assert "B3" not in codes(r)

    def test_o6_in_gof_gene_requires_functional_region(self):
        res = make_resources(roles=ONCOGENE)
        assert "O6" not in codes(svig.classify_row(row(REVEL="0.9"), res))
        assert "O6" in codes(svig.classify_row(row(REVEL="0.9", DOMAINS="Pfam:PF00001"), res))

    def test_spliceai_0_2_gives_o6(self):
        r = svig.classify_row(row(Consequence="intron_variant", HGVSp="", HGVSc="c.300+5G>A",
                                  **spliceai(0.25)), make_resources(roles=TSG))
        assert "O6" in codes(r)

    def test_synonymous_without_splicing_impact_is_b3_plus_b4(self):
        # "Synonymous variants with no predicted effect on splicing can be classified as
        #  likely benign (B3 and B4)" [S T2 B4]; both supporting in silico strength
        r = svig.classify_row(row(Consequence="synonymous_variant", HGVSp="p.Arg100=",
                                  gnomADe_NFE_AF="0.0003", **spliceai(0.01)), make_resources(roles=TSG))
        c = codes(r)
        assert c["B3"].points == -1 and c["B4"].points == -1

    def test_b4_intronic_needs_plus7_or_minus21(self):
        res = make_resources(roles=TSG)
        near = row(Consequence="intron_variant", HGVSp="", HGVSc="c.300+5G>A", **spliceai(0.0))
        deep = row(Consequence="intron_variant", HGVSp="", HGVSc="c.300+25G>A", **spliceai(0.0))
        assert "B4" not in codes(svig.classify_row(near, res))
        assert "B4" in codes(svig.classify_row(deep, res))

    def test_synonymous_in_splice_region_does_not_get_b4(self):
        r = svig.classify_row(row(Consequence="synonymous_variant&splice_region_variant",
                                  HGVSp="p.Arg100=", **spliceai(0.0)), make_resources(roles=TSG))
        assert "B4" not in codes(r)


# ===========================================================================
# O7 / O8 / O4 / O9 / B5
# ===========================================================================

def hotspot(total, same, alt="H"):
    return {"hugoSymbol": "GENE1", "residue": "R100", "type": "single residue",
            "aminoAcidPosition": {"start": 100, "end": 100}, "tumorCount": total,
            "variantAminoAcid": {alt: same}}


class TestHotspotsConstraintCosmic:
    @pytest.mark.parametrize("total,same,points", [
        (60, 12, 4),   # >=50 at residue and >=10 same change
        (40, 12, 2),   # <50 at residue, >=10 same change
        (40, 5, 1),    # 2-9 same change
        (40, 1, None),
    ])
    def test_o7_tiers(self, total, same, points):
        r = svig.classify_row(row(), make_resources(roles=TSG, hotspots=[hotspot(total, same)]))
        assert (codes(r)["O7"].points if "O7" in codes(r) else None) == points

    def test_o8_threshold_3_09(self):
        res_hi = make_resources(mis_z={"GENE1": {"_mane": 3.2}})
        res_lo = make_resources(mis_z={"GENE1": {"_mane": 2.9}})
        assert "O8" in codes(svig.classify_row(row(), res_hi))
        assert "O8" not in codes(svig.classify_row(row(), res_lo))

    @pytest.mark.parametrize("n,points", [
        (51, 4),     # >50 COSMIC samples (GENIE >10 x5; the guideline's COSMIC example)
        (25, 2),     # 25-50
        (24, None),
    ])
    def test_o4_missense_cosmic_thresholds(self, tmp_path, n, points):
        rows = [("GENE1", "T1", "R100H", "missense", 100, 100, n)]
        res = make_resources(roles=TSG, cosmic_rows=rows, tmp_path=tmp_path)
        r = svig.classify_row(row(), res)
        assert (codes(r)["O4"].points if "O4" in codes(r) else None) == points

    def test_o4_counts_amino_acid_change_across_transcripts_without_double_counting(self, tmp_path):
        # same change reported on two COSMIC transcripts: take one transcript, not the sum
        rows = [("GENE1", "T1", "R100H", "missense", 100, 100, 30),
                ("GENE1", "T2", "R100H", "missense", 100, 100, 30)]
        res = make_resources(roles=TSG, cosmic_rows=rows, tmp_path=tmp_path)
        assert codes(svig.classify_row(row(), res))["O4"].points == 2

    def test_o4_truncating_counts_same_position_or_downstream(self, tmp_path):
        rows = [("GENE1", "T1", "R100*", "truncating", 100, 100, 150),
                ("GENE1", "T1", "Q120*", "truncating", 120, 120, 120),
                ("GENE1", "T1", "S90*", "truncating", 90, 90, 5000)]   # upstream: not counted
        res = make_resources(roles=TSG, cosmic_rows=rows, tmp_path=tmp_path)
        r = svig.classify_row(row(Consequence="stop_gained", HGVSp="p.Arg100Ter"), res)
        assert codes(r)["O4"].points == 4       # 270 samples > 250

    def test_o4_requires_o3(self, tmp_path):
        rows = [("GENE1", "T1", "R100H", "missense", 100, 100, 999)]
        res = make_resources(roles=TSG, cosmic_rows=rows, tmp_path=tmp_path)
        r = svig.classify_row(row(gnomADe_AFR_AF="0.0005"), res)
        assert "O4" not in codes(r)

    def test_o9_inframe_multi_residue_in_domain(self):
        r = svig.classify_row(row(Consequence="inframe_deletion", HGVSp="p.Glu746_Ala750del",
                                  DOMAINS="Pfam:PF07714"), make_resources(roles=ONCOGENE))
        assert codes(r)["O9"].points == 2

    def test_b5_inframe_in_low_complexity(self):
        r = svig.classify_row(row(Consequence="inframe_deletion", HGVSp="p.Gln50del",
                                  DOMAINS="Low_complexity_(Seg):seg"), make_resources(roles=ONCOGENE))
        assert codes(r)["B5"].points == -1


# ===========================================================================
# Combination rules and minimum-two-codes
# ===========================================================================

class TestCombinations:
    def test_o6_removed_with_o2(self):
        r = svig.classify_row(row(Consequence="stop_gained", HGVSp="p.Arg100Ter", **spliceai(0.5)),
                              make_resources(roles=TSG))
        assert "O2" in codes(r) and "O6" not in codes(r)

    def test_o8_removed_with_o7(self):
        res = make_resources(roles=TSG, hotspots=[hotspot(60, 12)], mis_z={"GENE1": {"_mane": 4.0}})
        r = svig.classify_row(row(), res)
        assert "O7" in codes(r) and "O8" not in codes(r)

    def test_o6_plus_o7_capped_at_4(self):
        res = make_resources(roles=TSG, hotspots=[hotspot(60, 12)])
        r = svig.classify_row(row(REVEL="0.95"), res)
        assert codes(r)["O7"].points == 4 and "O6" not in codes(r)

    def test_o7_removed_when_b1_applies(self):
        res = make_resources(roles=TSG, hotspots=[hotspot(60, 12)])
        r = svig.classify_row(row(gnomADe_NFE_AF="0.003"), res)
        assert "O7" not in codes(r)

    # gnomAD AFR exome AF 1e-4 ~ 3 alleles in a sub-population <100,000 alleles:
    # blocks O3 without reaching the 0.1% B1 threshold, isolating the code under test.
    NOT_RARE_NOT_COMMON = {"gnomADe_AF": "3e-06", "gnomADe_AFR_AF": "0.0001"}

    def test_single_very_strong_code_is_still_vus(self):
        # "a single line of 'very strong' evidence achieving +8 points will still be
        #  classified as a VUS unless a second evidence code ... can also be applied" [G §6.3]
        r = svig.classify_row(row(Consequence="stop_gained", HGVSp="p.Arg100Ter",
                                  **self.NOT_RARE_NOT_COMMON), make_resources(roles=TSG))
        assert list(codes(r)) == ["O2"] and r.score == 8
        assert r.classification == "VUS"

    def test_single_benign_code_stays_vus(self):
        # "a single evidence code supporting benignity (-1 point to -4 points) should still be
        #  classified as a VUS unless corroborating evidence from a second criteria" [G §6.3]
        r = svig.classify_row(row(REVEL="0.2", **self.NOT_RARE_NOT_COMMON), make_resources(roles=TSG))
        assert list(codes(r)) == ["B3"] and r.score == -1
        assert r.classification == "VUS"

    def test_two_benign_codes_reach_likely_benign(self):
        r = svig.classify_row(row(REVEL="0.2", gnomADe_NFE_AF="0.005"), make_resources(roles=TSG))
        assert set(codes(r)) == {"B1", "B3"} and r.score == -5
        assert r.classification == "Likely Benign"

    def test_o3_does_not_rescue_likely_benign_into_vus(self):
        # B3 + B4 = -2 (Likely Benign); absent-from-gnomAD O3 (+2) must not lift it to VUS
        r = svig.classify_row(row(Consequence="synonymous_variant", HGVSp="p.Arg100=", **spliceai(0.0)),
                              make_resources(roles=TSG))
        assert "O3" not in codes(r)
        assert r.classification == "Likely Benign"
