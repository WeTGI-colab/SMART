"""
Unit tests for scripts/svig_uk_classifier.py

Every expectation is taken from the ACGS / SVIG-UK 2025 guideline or its
Supplementary Material — never from observed classifier output — so these tests
document the rules the scorer is required to follow.

Covers:
  - normalise_protein()        — HGVS protein normalisation used for all lookups
  - variant_class()            — VEP consequence -> SVIG-UK variant category
  - stand-alone codes          — O1, B1 (>1%), B2 (null variant in an oncogene)
  - O2 PVS1 tree              — NMD position, protein fraction, exon skipping, +2T>C
  - O3 / B1 population rules   — absent, <=0.001%, >1 allele, >0.1%, CHIP correction
  - O6 / B3 / B4 computational — REVEL and SpliceAI thresholds, LOF restriction
  - O5 missense route, O7 hotspot tiers, O8 missense constraint, O4 case counting (COSMIC, GENIE), O9 / B5
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
# O2 — PVS1 tree with transcript structure [S T1 O2; Abou Tayoun et al. 2018]
# ===========================================================================
#
# Synthetic transcript (plus strand; the minus strand is its mirror image):
#   exon  genomic        coding bp   CDS nt    residues   mRNA nt
#   1     1000-1100      90 (+10 UTR) 1-90     1-30       1-100
#   2     2000-2102      102          91-192   31-64      101-202
#   3     3000-3100      100          193-292  65-97      203-302
#   4     4000-4099      99           293-391  98-130     303-401
#   5     5000-5300      209          392-600  131-199+*  402-701
# Protein length 199 aa. Exon 2 and 4 skipping keep the frame; exon 3 does not.

MIRROR = 6300
EXONS = [(1000, 1100), (2000, 2102), (3000, 3100), (4000, 4099), (5000, 5300)]
CDS = (1010, 5209)


def model_line(strand):
    exons, cds = EXONS, CDS
    if strand == "-":
        exons = sorted((MIRROR - e, MIRROR - s) for s, e in exons)
        cds = (MIRROR - CDS[1], MIRROR - CDS[0])
    starts = ",".join(str(s) for s, _ in exons) + ","
    ends = ",".join(str(e) for _, e in exons) + ","
    return "\t".join(map(str, [0, "NM_000001.1", "chr1", strand, exons[0][0], exons[-1][1],
                               cds[0], cds[1], len(exons), starts, ends])) + "\n"


def tx_resources(tmp_path, strand="+", genome_stop=False):
    res = make_resources(roles=TSG)
    path = tmp_path / "ncbiRefSeq.txt"
    path.write_text(model_line(strand) + model_line(strand).replace("chr1", "chr1_alt"))
    res.load_transcript_models(str(path))
    if genome_stop:
        seq = ["C"] * 6400
        if strand == "+":
            seq[4000:4003] = "TAA"  # first codon of exon 4, in frame once exon 3 is skipped
        else:
            seq[MIRROR - 4003:MIRROR - 4000] = "TTA"  # reverse complement of TAA
        fasta = tmp_path / "genome.fa"
        body = "".join("".join(seq)[i:i + 60] + "\n" for i in range(0, 6400, 60))
        fasta.write_text(">chr1\n" + body)
        (tmp_path / "genome.fa.fai").write_text("chr1\t6400\t6\t60\t61\n")
        res.genome = svig.Genome(str(fasta))
    return res


def tx_row(**kw):
    return row(Feature="NM_000001.1", Start_Position="3050", **kw)


def o2(result):
    return codes(result).get("O2")


class TestTranscriptModel:
    @pytest.mark.parametrize("strand", ["+", "-"])
    def test_structure_is_read_in_transcript_order(self, tmp_path, strand):
        tx = tx_resources(tmp_path, strand).transcript_for(tx_row())
        assert tx.protein_length == 199 and tx.utr5 == 10
        assert tx.exon_lengths == [100, 102, 100, 99, 300]
        assert tx.coding == [90, 102, 100, 99, 209]
        assert tx.exon_has_start(1) and tx.exon_has_stop(5)

    def test_versioned_accession_must_match(self, tmp_path):
        res = tx_resources(tmp_path)
        assert res.transcript_for(row(Feature="NM_000001.2", Start_Position="3050")) is None
        assert res.transcript_for(row(NM_Transcript="NM_000001", Start_Position="3050")) is not None

    def test_alt_contigs_are_ignored(self, tmp_path):
        assert len(tx_resources(tmp_path).transcripts["NM_000001"]) == 1

    @pytest.mark.parametrize("hgvsp,expected", [
        ("p.Arg100Ter", (100, 100)),
        ("ENSP1:p.Phe126LeufsTer5", (126, 130)),   # stop counted from the first changed residue
        ("p.Phe126LeufsTer?", (126, None)),
        ("p.Arg100*", (100, 100)),
        ("p.Missense100His", None),
    ])
    def test_ptc_residues(self, hgvsp, expected):
        assert svig.ptc_residues(hgvsp) == expected


class TestO2TranscriptStructure:
    # "if the premature termination codon (PTC) occurs in the last 50 bp's of the
    # penultimate exon or in the final exon, then NMD would generally not be expected"
    @pytest.mark.parametrize("strand", ["+", "-"])
    def test_ptc_upstream_of_penultimate_exon_is_very_strong(self, tmp_path, strand):
        r = svig.classify_row(tx_row(Consequence="stop_gained", HGVSp="p.Arg40Ter", EXON="2/5"),
                              tx_resources(tmp_path, strand))
        assert o2(r).strength == "VSTR"

    def test_ptc_early_in_penultimate_exon_is_very_strong(self, tmp_path):
        # residue 105 -> mRNA 323, 78 bp before the end of exon 4
        r = svig.classify_row(tx_row(Consequence="stop_gained", HGVSp="p.Arg105Ter", EXON="4/5"),
                              tx_resources(tmp_path))
        assert o2(r).strength == "VSTR"

    def test_ptc_in_last_50bp_of_penultimate_exon_escapes_nmd(self, tmp_path):
        # residue 125 -> mRNA 383, 18 bp before the end of exon 4; 75/199 residues lost (>10%)
        r = svig.classify_row(tx_row(Consequence="stop_gained", HGVSp="p.Arg125Ter", EXON="4/5"),
                              tx_resources(tmp_path))
        assert o2(r).strength == "STR"

    @pytest.mark.parametrize("residue,strength", [(140, "STR"), (195, "MOD")])
    def test_final_exon_weighted_by_protein_fraction(self, tmp_path, residue, strength):
        # Abou Tayoun 2018: NMD escape -> removes >10% of protein STR, otherwise MOD
        r = svig.classify_row(tx_row(Consequence="stop_gained", HGVSp=f"p.Arg{residue}Ter", EXON="5/5"),
                              tx_resources(tmp_path))
        assert o2(r).strength == strength

    def test_frameshift_uses_new_stop_position(self, tmp_path):
        res = tx_resources(tmp_path)
        fs = dict(Consequence="frameshift_variant", EXON="4/5", Reference_Allele="AC", Tumor_Seq_Allele2="A")
        # stop at residue 109: 65 bp before the end of exon 4 -> NMD
        assert o2(svig.classify_row(tx_row(HGVSp="p.Arg100GlyfsTer10", **fs), res)).strength == "VSTR"
        # stop at residue 119: 35 bp before the end of exon 4 -> escape, 85/199 residues altered
        assert o2(svig.classify_row(tx_row(HGVSp="p.Arg115GlyfsTer5", **fs), res)).strength == "STR"

    def test_exon_count_disagreeing_with_vep_falls_back(self, tmp_path):
        r = svig.classify_row(tx_row(Consequence="stop_gained", HGVSp="p.Arg195Ter", EXON="5/9"),
                              tx_resources(tmp_path))
        assert o2(r).strength == "VSTR"  # structure ignored: exon 5 of 9 predicted NMD

    def test_penultimate_exon_without_model_uses_lowest_strength(self):
        r = svig.classify_row(row(Consequence="stop_gained", HGVSp="p.Arg100Ter", EXON="9/10"),
                              make_resources(roles=TSG))
        assert o2(r).strength == "MOD"

    # Canonical splice: "it is assumed that the exon it is closest to is skipped"
    def splice(self, hgvsc, intron, **kw):
        return tx_row(Consequence="splice_donor_variant" if "+" in hgvsc else "splice_acceptor_variant",
                      HGVSp="", HGVSc=hgvsc, EXON="", INTRON=intron, **kw)

    @pytest.mark.parametrize("hgvsc,intron", [("c.192+1G>A", "2/4"), ("c.293-1G>A", "3/4")])
    def test_in_frame_exon_skipping_over_10_percent_is_strong(self, tmp_path, hgvsc, intron):
        # donor of intron 2 skips exon 2 (34 aa, 17%); acceptor of intron 3 skips exon 4 (33 aa, 17%)
        r = svig.classify_row(self.splice(hgvsc, intron), tx_resources(tmp_path))
        assert o2(r).strength == "STR" and "17%" in o2(r).reason

    def test_in_frame_exon_skipping_under_10_percent_is_moderate(self, tmp_path):
        res = tx_resources(tmp_path)
        tx = res.transcripts["NM_000001"][0]
        tx.protein_length = 400  # exon 2 = 34 aa = 8.5%
        assert o2(svig.classify_row(self.splice("c.192+1G>A", "2/4"), res)).strength == "MOD"

    @pytest.mark.parametrize("strand", ["+", "-"])
    def test_out_of_frame_skipping_with_nmd_is_very_strong(self, tmp_path, strand):
        r = svig.classify_row(self.splice("c.193-2A>G", "2/4"), tx_resources(tmp_path, strand, genome_stop=True))
        assert o2(r).strength == "VSTR"

    def test_out_of_frame_skipping_without_stop_escapes_nmd(self, tmp_path):
        res = tx_resources(tmp_path, genome_stop=True)
        res.genome.fetch = lambda chrom, s, e: "C" * (e - s)
        assert o2(svig.classify_row(self.splice("c.193-2A>G", "2/4"), res)).strength == "STR"

    def test_out_of_frame_skipping_without_genome_is_not_very_strong(self, tmp_path):
        r = svig.classify_row(self.splice("c.193-2A>G", "2/4"), tx_resources(tmp_path))
        assert o2(r).strength == "STR"

    def test_skipping_the_initiation_codon_exon_is_supporting(self, tmp_path):
        r = svig.classify_row(self.splice("c.90+1G>A", "1/4"), tx_resources(tmp_path))
        assert o2(r).strength == "SUPP"

    def test_final_intron_uses_lowest_strength(self, tmp_path):
        r = svig.classify_row(self.splice("c.392-1G>A", "4/4"), tx_resources(tmp_path))
        assert o2(r).strength == "MOD"

    @pytest.mark.parametrize("strand", ["+", "-"])
    def test_frameshift_without_protein_annotation_is_translated(self, tmp_path, strand):
        # Reference CDS on the synthetic genome reads CCC... (plus) or GGG... (minus); a
        # 1-bp deletion in exon 3 shifts the frame onto the TAA placed at exon 4's start
        res = tx_resources(tmp_path, strand, genome_stop=True)
        tx = res.transcripts["NM_000001"][0]
        first, stop = svig.translate_exonic_indel("c.200del", tx, res.genome)
        assert stop is not None and first <= stop
        r = svig.classify_row(self.splice("c.200del", "2/4"), res)  # VEP intron from the unshifted allele
        assert o2(r) is not None and "PTC" in o2(r).reason

    @pytest.mark.parametrize("hgvsc,expected", [
        ("c.192+1G>A", True), ("c.193-2A>G", True), ("c.192+2dup", True),
        ("c.192+5_192+8del", False),          # repeat-shifted: donor dinucleotide intact
        ("c.191_192+3del", True),             # spans the exon-intron boundary
        ("c.192+2_192+6del", True), ("c.192+5G>C", False), ("c.1799T>A", False),
    ])
    def test_canonical_site_from_hgvs(self, hgvsc, expected):
        assert svig.touches_canonical_splice_site(hgvsc) is expected

    def test_shifted_deletion_outside_canonical_site_gets_no_o2(self, tmp_path):
        r = svig.classify_row(self.splice("c.192+5_192+8del", "2/4"), tx_resources(tmp_path))
        assert o2(r) is None

    def test_nearby_cryptic_site_uses_lowest_strength(self, tmp_path):
        # "the lowest strength of O2 should be applied when different in-frame and
        # out-of-frame scenarios are possible"
        r = svig.classify_row(self.splice("c.193-2A>G", "2/4", SpliceAI_pred_DS_AG="0.45",
                                          SpliceAI_pred_DP_AG="-9"),
                              tx_resources(tmp_path, genome_stop=True))
        assert o2(r).strength == "MOD"


# ===========================================================================
# O4 from GENIE [S T1 O4] — synthetic transcript above; residue r of exon 2 starts at
# genomic 1-based 2001 + 3*(r-31), exon 5 residue r at 5001 + 3*(r-131) + 2
# ===========================================================================

# Minimal OncoTree: TISSUE > MYELOID > MNM > AML > AMLMRC; MNM > MDS; LYMPHOID > DLBCLNOS; LUNG > LUAD
ONCOTREE = {"tumorTypes": [
    {"code": "MYELOID", "tissue": "Myeloid", "parent": "TISSUE"},
    {"code": "MNM", "tissue": "Myeloid", "parent": "MYELOID"},
    {"code": "AML", "tissue": "Myeloid", "parent": "MNM"},
    {"code": "AMLMRC", "tissue": "Myeloid", "parent": "AML"},
    {"code": "MDS", "tissue": "Myeloid", "parent": "MNM"},
    {"code": "LYMPHOID", "tissue": "Lymphoid", "parent": "TISSUE"},
    {"code": "DLBCLNOS", "tissue": "Lymphoid", "parent": "LYMPHOID"},
    {"code": "LUNG", "tissue": "Lung", "parent": "TISSUE"},
    {"code": "LUAD", "tissue": "Lung", "parent": "LUNG"},
]}
DEFAULT_CODE = {"Myeloid": ("AML", "Myeloid"), "Lymphoid": ("DLBCLNOS", "Lymphoid"), "Solid": ("LUAD", "Lung")}


def genie_resources(tmp_path, entries, strand="+", level="lineage", old_format=False):
    res = tx_resources(tmp_path, strand)
    res.gene_roles = {"GENE1": {"TSG"}, "EZH2": {"TSG", "ONC"}}
    path = tmp_path / "genie.tsv"
    header = "gene\tvariant_class\tchrom\tstart\tend\tref\talt\tref_aa\talt_aa\tlineage\tn_patients"
    lines = [header if old_format else header + "\toncotree_code\ttissue"]
    for e in entries:
        e = tuple(e)
        if len(e) == 11 and not old_format:
            e = e + DEFAULT_CODE[e[9]]
        lines.append("\t".join(map(str, e)))
    path.write_text("\n".join(lines) + "\n")
    res.genie = svig.GenieO4(str(path))
    (tmp_path / "oncotree.json").write_text(json.dumps(ONCOTREE))
    res.oncotree = svig.OncoTree(str(tmp_path / "oncotree.json"))
    res.o4_level = level
    return res


def mis(pos, lineage, n, ref_aa="R", alt_aa="H", gene="GENE1"):
    return (gene, "missense", "chr1", pos, pos, "G", "A", ref_aa, alt_aa, lineage, n)


def trunc(pos, n, lineage="Myeloid"):
    return ("GENE1", "truncating", "chr1", pos, pos, "C", "T", "", "", lineage, n)


def genie_o4(res, **kw):
    base = dict(TUMOUR_TYPE="Myeloid")
    base.update(kw)
    r = svig.classify_row(tx_row(**base), res)
    return codes(r).get("O4"), r


class TestTranscriptCoordinates:
    @pytest.mark.parametrize("strand,pos", [("+", 2001), ("-", MIRROR - 2000)])
    def test_residue_at(self, tmp_path, strand, pos):
        tx = tx_resources(tmp_path, strand).transcripts["NM_000001"][0]
        assert tx.residue_at(pos) == 31
        assert tx.residue_at(1500 if strand == "+" else MIRROR - 1500) is None  # intronic

    def test_maf_allele_from_vcf_style(self):
        assert svig.maf_allele(row(Chromosome="1", POS="100", REF="CAG", ALT="C")) == ("chr1", 101, 102, "AG", "-")
        assert svig.maf_allele(row(POS="100", REF="C", ALT="CGG")) == ("chr1", 100, 101, "-", "GG")
        assert svig.maf_allele(row(POS="100", REF="C", ALT="T")) == ("chr1", 100, 100, "C", "T")

    def test_left_align(self, tmp_path):
        seq = "ACGTTTTG"  # positions 1-8
        fasta = tmp_path / "g.fa"
        fasta.write_text(">chr1\n" + seq + "\n")
        (tmp_path / "g.fa.fai").write_text(f"chr1\t{len(seq)}\t6\t{len(seq)}\t{len(seq) + 1}\n")
        g = svig.Genome(str(fasta))
        # deleting the T at position 7 is the same as deleting the T at position 4
        assert svig.left_align_maf(g, "chr1", 7, 7, "T", "-") == (4, 4, "T", "-")
        # inserting TT after position 7 -> after position 3
        assert svig.left_align_maf(g, "chr1", 7, 8, "-", "TT") == (3, 4, "-", "TT")


class TestO4Genie:
    def test_missense_same_amino_acid_change_haematological_entries(self, tmp_path):
        # myeloid + lymphoid = 11 > 10 -> strong; solid entries are off-target
        res = genie_resources(tmp_path, [mis(2001, "Myeloid", 6), mis(2001, "Lymphoid", 5),
                                         mis(2001, "Solid", 400)])
        o4, _ = genie_o4(res, HGVSp="p.Arg31His")
        assert o4.strength == "STR" and "GENIE 11 patients" in o4.reason

    @pytest.mark.parametrize("n,strength", [(4, None), (5, "MOD"), (10, "MOD"), (11, "STR")])
    def test_missense_thresholds(self, tmp_path, n, strength):
        res = genie_resources(tmp_path, [mis(2001, "Myeloid", n)])
        o4, _ = genie_o4(res, HGVSp="p.Arg31His")
        assert (o4.strength if o4 else None) == strength

    def test_missense_counts_any_nucleotide_change_but_not_other_amino_acids(self, tmp_path):
        # "regardless of the underlying nucleotide change": another base of codon 31
        res = genie_resources(tmp_path, [mis(2002, "Myeloid", 6), mis(2003, "Myeloid", 5),
                                         mis(2001, "Myeloid", 50, alt_aa="C"),
                                         mis(2004, "Myeloid", 50, ref_aa="R", alt_aa="H")])
        o4, _ = genie_o4(res, HGVSp="p.Arg31His")
        assert "GENIE 11 patients" in o4.reason

    def test_solid_tumour_counts_solid_entries(self, tmp_path):
        res = genie_resources(tmp_path, [mis(2001, "Myeloid", 3), mis(2001, "Solid", 11)])
        o4, _ = genie_o4(res, HGVSp="p.Arg31His", TUMOUR_TYPE="solid")
        assert o4.strength == "STR" and "GENIE 11 patients" in o4.reason

    @pytest.mark.parametrize("level,n,label", [
        ("code", 5, "AML and subtypes"),       # AML 3 + AMLMRC 2
        ("tissue", 9, "tissue Myeloid"),       # + MDS 4
        ("lineage", 14, "Lymphoid+Myeloid"),   # + DLBCLNOS 5
        ("all", 34, "all tumour types"),       # + LUAD 20 (off-target)
    ])
    def test_on_target_levels_for_an_oncotree_code(self, tmp_path, level, n, label):
        # "'on-target' entries may be specific to AML, to myeloid disorders or to all
        # haematopoietic disorders (myeloid and lymphoid)"
        entries = [mis(2001, "Myeloid", 3) + ("AML", "Myeloid"), mis(2001, "Myeloid", 2) + ("AMLMRC", "Myeloid"),
                   mis(2001, "Myeloid", 4) + ("MDS", "Myeloid"), mis(2001, "Lymphoid", 5) + ("DLBCLNOS", "Lymphoid"),
                   mis(2001, "Solid", 20) + ("LUAD", "Lung")]
        res = genie_resources(tmp_path, entries, level=level)
        o4, _ = genie_o4(res, HGVSp="p.Arg31His", TUMOUR_TYPE="AML")
        assert f"GENIE {n} patients" in o4.reason and label in o4.reason

    def test_level_finer_than_tumour_type_falls_back(self, tmp_path):
        entries = [mis(2001, "Myeloid", 3) + ("AML", "Myeloid"), mis(2001, "Myeloid", 4) + ("MDS", "Myeloid"),
                   mis(2001, "Lymphoid", 5) + ("DLBCLNOS", "Lymphoid")]
        res = genie_resources(tmp_path, entries, level="code")
        # "myeloid" is also the OncoTree tissue node MYELOID: its subtree is the myeloid tissue
        assert "GENIE 7 patients" in genie_o4(res, HGVSp="p.Arg31His", TUMOUR_TYPE="myeloid")[0].reason
        assert "Lymphoid+Myeloid" in genie_o4(res, HGVSp="p.Arg31His", TUMOUR_TYPE="haem")[0].reason

    def test_patient_with_on_target_and_off_target_samples_counts_once(self, tmp_path):
        # "multiple samples from the same patient are only counted as a single entry"
        entries = [mis(2001, "Myeloid", 4) + ("AML", "Myeloid"),
                   ("GENE1", "missense", "chr1", 2001, 2001, "G", "A", "R", "H", "Myeloid|Solid", 2,
                    "AML|LUAD", "Myeloid|Lung")]
        res = genie_resources(tmp_path, entries, level="code")
        assert "GENIE 6 patients" in genie_o4(res, HGVSp="p.Arg31His", TUMOUR_TYPE="AML")[0].reason
        res.o4_level = "all"
        assert "GENIE 6 patients" in genie_o4(res, HGVSp="p.Arg31His", TUMOUR_TYPE="AML")[0].reason

    def test_unknown_tumour_type_counts_all_entries(self, tmp_path):
        res = genie_resources(tmp_path, [mis(2001, "Myeloid", 3), mis(2001, "Solid", 8)])
        o4, _ = genie_o4(res, HGVSp="p.Arg31His", TUMOUR_TYPE="TP53 [Other]")
        assert "GENIE 11 patients" in o4.reason and "all tumour types" in o4.reason

    def test_old_table_without_codes_supports_tissue_level(self, tmp_path):
        res = genie_resources(tmp_path, [mis(2001, "Myeloid", 6), mis(2001, "Lymphoid", 50)],
                              level="tissue", old_format=True)
        assert "GENIE 6 patients" in genie_o4(res, HGVSp="p.Arg31His")[0].reason

    @pytest.mark.parametrize("level", ["lineage", "all"])
    def test_lineage_dependent_gene_never_counts_the_other_lineage(self, tmp_path, level):
        res = genie_resources(tmp_path, [mis(2001, "Myeloid", 6, gene="EZH2"),
                                         mis(2001, "Lymphoid", 50, gene="EZH2")], level=level)
        assert genie_o4(res, Hugo_Symbol="EZH2", HGVSp="p.Arg31His")[0].strength == "MOD"
        # an OncoTree code resolves the lineage too: DLBCLNOS is lymphoid
        o4, _ = genie_o4(res, Hugo_Symbol="EZH2", HGVSp="p.Arg31His", TUMOUR_TYPE="DLBCLNOS")
        assert "GENIE 50 patients" in o4.reason

    def test_lineage_dependent_gene_counts_own_lineage_only(self, tmp_path):
        # EZH2: tumour suppressor in myeloid, activating in lymphoid disorders [G §7.2]
        res = genie_resources(tmp_path, [mis(2001, "Myeloid", 6, gene="EZH2"),
                                         mis(2001, "Lymphoid", 50, gene="EZH2")])
        o4, _ = genie_o4(res, Hugo_Symbol="EZH2", HGVSp="p.Arg31His")
        assert o4.strength == "MOD"

    def test_truncating_same_position_or_downstream(self, tmp_path):
        # residue 35 (upstream) not counted; residues 40, 60 and 150 counted; 3'UTR ignored
        res = genie_resources(tmp_path, [trunc(2013, 100), trunc(2028, 5), trunc(2088, 20),
                                         trunc(5059, 26), trunc(5250, 100)])
        o4, _ = genie_o4(res, Consequence="stop_gained", HGVSp="p.Arg40Ter", EXON="2/5")
        assert o4.strength == "STR" and "GENIE 51 patients" in o4.reason

    @pytest.mark.parametrize("n,strength", [(9, None), (10, "SUPP"), (19, "SUPP"), (20, "MOD"),
                                            (50, "MOD"), (51, "STR")])
    def test_truncating_thresholds(self, tmp_path, n, strength):
        res = genie_resources(tmp_path, [trunc(2088, n)])
        o4, _ = genie_o4(res, Consequence="stop_gained", HGVSp="p.Arg40Ter", EXON="2/5")
        assert (o4.strength if o4 else None) == strength

    def test_inframe_deletion_counts_equal_or_smaller_within_region(self, tmp_path):
        dels = [("GENE1", "inframe_del", "chr1", 3011, 3016, "AAAAAA", "-", "", "", "Myeloid", 12),
                ("GENE1", "inframe_del", "chr1", 3012, 3014, "AAA", "-", "", "", "Myeloid", 9),
                ("GENE1", "inframe_del", "chr1", 3008, 3013, "CCCAAA", "-", "", "", "Myeloid", 100)]
        res = genie_resources(tmp_path, dels)
        o4, _ = genie_o4(res, Consequence="inframe_deletion", HGVSp="p.Lys2_Lys3del",
                         POS="3010", REF="CAAAAAA", ALT="C")
        assert o4.strength == "MOD" and "GENIE 21 patients" in o4.reason

    def test_inframe_insertion_exact_match_only(self, tmp_path):
        ins = [("GENE1", "inframe_ins", "chr1", 3010, 3011, "-", "GGG", "", "", "Myeloid", 51),
               ("GENE1", "inframe_ins", "chr1", 3010, 3011, "-", "GGA", "", "", "Myeloid", 100)]
        res = genie_resources(tmp_path, ins)
        o4, _ = genie_o4(res, Consequence="inframe_insertion", HGVSp="p.Lys2dup",
                         POS="3010", REF="C", ALT="CGGG")
        assert o4.strength == "STR" and "GENIE 51 patients" in o4.reason

    def test_splice_exact_nucleotide_change(self, tmp_path):
        spl = [("GENE1", "splice", "chr1", 3000, 3000, "G", "A", "", "", "Myeloid", 11),
               ("GENE1", "splice", "chr1", 3000, 3000, "G", "T", "", "", "Myeloid", 100)]
        res = genie_resources(tmp_path, spl)
        o4, _ = genie_o4(res, Consequence="splice_acceptor_variant", HGVSp="", HGVSc="c.193-1G>A",
                         INTRON="2/4", EXON="", POS="3000", REF="G", ALT="A")
        assert o4.strength == "STR" and "GENIE 11 patients" in o4.reason

    def test_genie_preferred_over_cosmic(self, tmp_path):
        res = genie_resources(tmp_path, [mis(2001, "Myeloid", 3)])
        path = tmp_path / "cosmic.tsv"
        path.write_text("gene\ttranscript\tprotein\tmut_class\taa_start\taa_end\tn_samples\n"
                        "GENE1\tENST1\tR31H\tmissense\t31\t31\t500\n")
        res.cosmic = svig.CosmicO4(str(path))
        o4, _ = genie_o4(res, HGVSp="p.Arg31His")
        assert o4 is None


# ===========================================================================
# O5 missense route [S T1 O5] — residue 31 of the synthetic transcript = genomic 2001-2003
# ===========================================================================

class FakeRevel:
    def __init__(self, scores):  # {(pos, aaalt): revel}
        self.scores = scores

    def at(self, chrom, pos):
        return [("G", "A", "R", alt, sc) for (p, alt), sc in self.scores.items() if p == pos]


def o5_resources(tmp_path, references, revel):
    res = tx_resources(tmp_path)
    for transcript, protein, source in references:
        res.add_o5_reference("GENE1", transcript, protein, source)
    res.revel = FakeRevel(revel)
    return res


def o5(res, vua_revel, protein="p.Arg31His"):
    r = tx_row(HGVSp=protein, REVEL=str(vua_revel))
    return svig.rule_o5(r, res, "missense", svig.normalise_protein(protein))


class TestO5:
    @pytest.mark.parametrize("vua,ref,strength", [
        (0.79, 0.80, "MOD"),    # equivalent: difference < 0.02
        (0.85, 0.80, "MOD"),    # more deleterious than the reference
        (0.78, 0.95, "MOD"),    # both > 0.773
        (0.72, 0.90, "SUPP"),   # > 0.7 and < 0.773, less deleterious
        (0.65, 0.90, None),
    ])
    def test_revel_comparison(self, tmp_path, vua, ref, strength):
        res = o5_resources(tmp_path, [("", "p.Arg31Cys", "CGI")], {(2001, "C"): ref})
        e = o5(res, vua)
        assert (e.strength if e else None) == strength

    def test_reference_revel_is_highest_nucleotide_change_for_that_amino_acid(self, tmp_path):
        res = o5_resources(tmp_path, [("", "p.Arg31Cys", "CGI")], {(2001, "C"): 0.80, (2003, "C"): 0.95})
        assert o5(res, 0.79) is not None and "REVEL 0.950" in o5(res, 0.79).reason

    def test_same_amino_acid_change_is_not_o5(self, tmp_path):
        res = o5_resources(tmp_path, [("", "p.Arg31His", "CGI")], {(2001, "H"): 0.9})
        assert o5(res, 0.9) is None

    @pytest.mark.parametrize("protein", ["p.Arg32Cys", "p.Lys31Cys"])
    def test_reference_without_transcript_needs_same_residue_and_amino_acid(self, tmp_path, protein):
        res = o5_resources(tmp_path, [("", protein, "CGI")], {(2001, "C"): 0.9})
        assert o5(res, 0.9) is None

    def test_reference_on_another_transcript_is_mapped_through_the_genome(self, tmp_path):
        res = o5_resources(tmp_path, [("NM_000002.1", "p.Arg41Cys", "SVIG-UK canonical list")],
                           {(2001, "C"): 0.9})
        # NM_000002 starts translation in an extra 20-bp exon and codes all of exon 1 (+10 bp):
        # 30 extra coding bp, so its residue 41 is residue 31 here
        path = tmp_path / "second.txt"
        exons = [(500, 520)] + EXONS
        path.write_text("\t".join(map(str, [0, "NM_000002.1", "chr1", "+", 500, 5300, 500, CDS[1],
                                            len(exons), ",".join(str(a) for a, _ in exons) + ",",
                                            ",".join(str(b) for _, b in exons) + ","])) + "\n")
        res.load_transcript_models(str(path))
        assert o5(res, 0.9).strength == "MOD"
        res.transcripts.pop("NM_000002")
        assert o5(res, 0.9) is None  # cannot map without the reference transcript

    def test_combinations(self, tmp_path):
        res = o5_resources(tmp_path, [("", "p.Arg31Cys", "CGI")], {(2001, "C"): 0.8})
        res.mis_z = {"GENE1": {"_mane": 4.0}}
        res.hotspots_single[("GENE1", 31)] = {"hugoSymbol": "GENE1", "aminoAcidPosition": {"start": 31},
                                              "variantAminoAcid": {"H": 3}, "type": "single residue"}
        r = svig.classify_row(tx_row(HGVSp="p.Arg31His", REVEL="0.85"), res)
        c = codes(r)
        # O5 (MOD +2) beats O7 (SUPP +1); O8 not with O5
        assert "O5" in c and "O7" not in c and "O8" not in c


# ===========================================================================
# O3 CHIP correction [S T1 O3]
# ===========================================================================

def chip_resources(tmp_path, hist_exomes, an_exomes=1_400_000, genomes=None):
    res = make_resources(roles={"TET2": {"TSG"}, "GENE1": {"TSG"}})
    lines = ["chrom\tpos\tref\talt\tdataset\tfilter\tAC\tAN\tnhomalt\tab_hist_alt_bin_freq"]
    lines.append(f"chr4\t105243568\tA\tG\texomes\tPASS\t{sum(hist_exomes)}\t{an_exomes}\t0\t"
                 + "|".join(map(str, hist_exomes)))
    if genomes:
        hist, an = genomes
        lines.append(f"chr4\t105243568\tA\tG\tgenomes\tPASS\t{sum(hist)}\t{an}\t0\t" + "|".join(map(str, hist)))
    path = tmp_path / "chip.tsv"
    path.write_text("\n".join(lines) + "\n")
    res.load_chip_allele_balance(str(path))
    return res


def chip_row(**kw):
    # total AF 5e-5 with NFE 6e-5 (>100,000 alleles): O3 fails on frequency alone
    base = dict(Hugo_Symbol="TET2", Chromosome="chr4", POS="105243568", REF="A", ALT="G",
                gnomADe_AF="5e-05", gnomADe_NFE_AF="6e-05")
    base.update(kw)
    return row(**base)


def hist(**bins):
    h = [0] * 20
    for k, v in bins.items():
        h[int(k[1:])] = v
    return h


class TestO3Chip:
    def test_low_vaf_alleles_only_restores_o3(self, tmp_path):
        # every gnomAD carrier at VAF 0.05-0.20: all excluded as likely CHIP
        res = chip_resources(tmp_path, hist(b1=40, b2=20, b3=10))
        assert codes(svig.classify_row(chip_row(), res))["O3"].strength == "MOD"

    def test_recalculated_frequency_at_or_below_0_001_percent_is_supporting(self, tmp_path):
        # 5 alleles at VAF 0.30-0.65 over 1,400,000 + 150,000 alleles = 3.2e-6
        res = chip_resources(tmp_path, hist(b2=60, b8=3, b12=1), genomes=(hist(b10=1), 150_000))
        assert codes(svig.classify_row(chip_row(), res))["O3"].strength == "SUPP"

    def test_vaf_bins_are_0_30_to_0_65(self, tmp_path):
        # bin 5 (0.25-0.30) and bin 13 (0.65-0.70) are outside the germline-like window
        res = chip_resources(tmp_path, hist(b5=500, b13=500))
        assert codes(svig.classify_row(chip_row(), res))["O3"].strength == "MOD"

    def test_germline_like_alleles_keep_o3_off(self, tmp_path):
        res = chip_resources(tmp_path, hist(b9=40, b10=30))
        r = svig.classify_row(chip_row(), res)
        assert "O3" not in codes(r) and any("still not applicable" in n for n in r.notes)

    def test_only_chip_genes(self, tmp_path):
        res = chip_resources(tmp_path, hist(b1=70))
        assert "O3" not in codes(svig.classify_row(chip_row(Hugo_Symbol="GENE1"), res))

    def test_requires_total_af_above_0_001_percent(self, tmp_path):
        # O3 fails through the sub-population rule, but total AF is not >0.001%
        res = chip_resources(tmp_path, hist(b1=3))
        r = svig.classify_row(chip_row(gnomADe_AF="8e-06", gnomADe_NFE_AF="",
                                       gnomADe_AFR_AF="0.0001"), res)
        assert "O3" not in codes(r)

    def test_not_applied_when_b1_applies(self, tmp_path):
        res = chip_resources(tmp_path, hist(b1=5000))
        r = svig.classify_row(chip_row(gnomADe_AF="0.002", gnomADe_NFE_AF="0.002"), res)
        assert "O3" not in codes(r) and "B1" in codes(r)

    def test_missing_table_is_reported(self):
        res = make_resources(roles={"TET2": {"TSG"}})
        r = svig.classify_row(chip_row(), res)
        assert "O3" not in codes(r) and any("allele-balance table not provided" in n for n in r.notes)


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


# ===========================================================================
# Tier 4 output — two-row header (field names; "description | source | version")
# ===========================================================================

class TestTier4:
    def write_maf(self, tmp_path, rows):
        cols = sorted({k for r in rows for k in r})
        path = tmp_path / "tier1.maf"
        path.write_text("#SMART_VERSION 2.0.0\n" + "\t".join(cols) + "\n" +
                        "".join("\t".join(str(r.get(c, "")) for c in cols) + "\n" for r in rows))
        return path

    def run(self, tmp_path, rows, *extra):
        maf = self.write_maf(tmp_path, rows)
        out = tmp_path / "tier4.tsv"
        svig.main(["--maf", str(maf), "--tier4", str(out), "--canonical", CANONICAL, *extra])
        lines = out.read_text().splitlines()
        header, meta = lines[0].split("\t"), lines[1].split("\t")
        records = [dict(zip(header, l.split("\t"))) for l in lines[2:]]
        return header, meta, records

    def test_two_row_header_documents_every_column(self, tmp_path):
        header, meta, _ = self.run(tmp_path, [row(Hugo_Symbol="GENE1")])
        assert len(header) == len(meta)
        assert all(m.count(" | ") == 2 for m in meta)
        for code in ["O1", "O2", "O11", "B1", "B7"]:
            assert code in header
        assert "ACGS/SVIG-UK" in meta[header.index("SVIG_UK_CLASS")]
        assert "SMART 2.0.0" in meta[header.index("SVIG_UK_CLASS")]
        assert "not automated" in meta[header.index("O11")]

    def test_one_column_per_code_with_strength_and_points(self, tmp_path):
        null = row(Hugo_Symbol="TP53", Consequence="stop_gained", HGVSp="p.Arg100Ter", EXON="3/10",
                   TUMOUR_TYPE="solid")
        _, _, records = self.run(tmp_path, [null])
        r = records[0]
        assert r["O3"] == "MOD +2" and r["B1"] == "" and r["TUMOUR_GROUP"] == "solid"
        assert r["SVIG_UK_CODES"].startswith(("O2", "O3"))

    def test_standalone_codes(self, tmp_path):
        common = row(gnomADe_AF="0.2")
        _, _, records = self.run(tmp_path, [common])
        assert records[0]["B1"] == "stand-alone" and records[0]["SVIG_UK_CLASS"] == "Benign"
        assert records[0]["SVIG_UK_SCORE"] == ""

    def test_out_or_tier4_required(self, tmp_path):
        maf = self.write_maf(tmp_path, [row()])
        with pytest.raises(SystemExit):
            svig.main(["--maf", str(maf)])


# ===========================================================================
# Gene-specific population thresholds for hereditary cancer genes [S T1 O3; S T2 B1]
# (ClinGen VCEP rules from reference_lists/svig_uk_vcep_population_thresholds.tsv)
# ===========================================================================

VCEP = os.path.join(SCRIPTS_DIR, "..", "reference_lists", "svig_uk_vcep_population_thresholds.tsv")


def vcep_resources(gene):
    res = make_resources(roles={gene: {"TSG"}})
    res.load_vcep_thresholds(VCEP)
    return res


def classify_af(gene, **afs):
    r = svig.classify_row(row(Hugo_Symbol=gene, **{k: str(v) for k, v in afs.items()}), vcep_resources(gene))
    return codes(r), r


class TestVcepThresholds:
    def test_table_loads_expected_genes(self):
        res = vcep_resources("TP53")
        assert {"TP53", "RUNX1", "BRCA1", "BRCA2", "MLH1", "MSH2", "MSH6", "PMS2", "APC", "PTEN",
                "CDH1", "ATM", "PALB2", "VHL", "DICER1"} <= set(res.vcep)
        assert res.vcep["TP53"].label == "ClinGen TP53 VCEP GN009 v2.4"

    def test_faf95_is_below_the_point_estimate(self):
        assert 0 < svig.faf95(28, 1_112_004) < 28 / 1_112_004
        assert svig.faf95(0, 1000) == 0

    # O3 — "ClinGen VCEP gene-specific guidelines must be consulted (if they exist) to establish a cut-off"
    def test_tp53_o3_uses_vcep_cut_off(self):
        # 0.002%: above the generic 0.001% but below TP53 PM2_Supporting (<0.003%)
        c, _ = classify_af("TP53", gnomADe_AF=2e-5, gnomADe_NFE_AF=2.5e-5)
        assert c["O3"].strength == "SUPP" and "GN009" in c["O3"].reason
        c, _ = classify_af("TP53", gnomADe_AF=5e-5, gnomADe_NFE_AF=6e-5)
        assert "O3" not in c

    def test_multiallelic_group_threshold(self):
        # TP53: a group with multiple alleles must stay <0.004% (AFR 1e-4 with ~3 alleles fails)
        c, _ = classify_af("TP53", gnomADe_AF=2e-5, gnomADe_AFR_AF=1e-4)
        assert "O3" not in c

    def test_brca1_o3_only_when_absent(self):
        assert "O3" not in classify_af("BRCA1", gnomADe_AF=1e-6, gnomADe_NFE_AF=1e-6)[0]
        assert classify_af("BRCA1")[0]["O3"].strength == "MOD"

    def test_apc_single_allele_rule(self):
        assert "O3" in classify_af("APC", gnomADe_AF=6e-7)[0]          # ~1 allele: <0.001%
        assert "O3" not in classify_af("APC", gnomADe_AF=5e-6)[0]      # >1 allele: must be <=0.0003%

    def test_dicer1_no_more_than_one_allele_per_group(self):
        assert "O3" not in classify_af("DICER1", gnomADe_AF=2e-6, gnomADe_SAS_AF=3e-5)[0]

    # B1 — ">1% (or as specified by respective expert panel group)"
    def test_pms2_ba1_is_standalone(self):
        c, r = classify_af("PMS2", gnomADe_AF=0.004, gnomADe_NFE_AF=0.004)
        assert c["B1"].strength == "STANDALONE" and r.classification == "Benign"

    def test_msh6_bs1_is_strong_below_generic_threshold(self):
        # 0.1% is not >0.1% (generic), but meets MSH6 BS1 (GrpMax FAF >=0.022%)
        c, _ = classify_af("MSH6", gnomADe_AF=0.001, gnomADe_NFE_AF=0.001)
        assert c["B1"].strength == "STR"

    def test_vhl_ba1_at_low_frequency(self):
        c, _ = classify_af("VHL", gnomADe_AF=0.0002, gnomADe_NFE_AF=0.0002)
        assert c["B1"].strength == "STANDALONE" and "GN078" in c["B1"].reason

    def test_runx1_requires_five_alleles(self):
        # AMR genomes (15,294 alleles) at 0.02% is ~3 alleles: below the RUNX1 minimum of 5
        assert "B1" not in classify_af("RUNX1", gnomADg_AF=0.0002, gnomADg_AMR_AF=0.0002)[0]

    def test_genes_without_vcep_keep_generic_rules(self):
        c, _ = classify_af("GENE1", gnomADe_AF=0.002, gnomADe_NFE_AF=0.002)
        assert c["B1"].strength == "STR" and "0.1%" in c["B1"].reason

    def test_chip_correction_uses_gene_cut_off(self, tmp_path):
        res = vcep_resources("TP53")
        path = tmp_path / "chip.tsv"
        path.write_text("chrom\tpos\tref\talt\tdataset\tfilter\tAC\tAN\tnhomalt\tab_hist_alt_bin_freq\n"
                        "chr17\t7675237\tC\tT\texomes\tPASS\t80\t1400000\t0\t" + "|".join(["0"] * 8 + ["1"] + ["0"] * 11) + "\n")
        res.load_chip_allele_balance(str(path))
        r = svig.classify_row(row(Hugo_Symbol="TP53", Chromosome="chr17", POS="7675237", REF="C", ALT="T",
                                  gnomADe_AF="5e-05", gnomADe_NFE_AF="6e-05"), res)
        assert codes(r)["O3"].strength == "SUPP" and "GN009" in codes(r)["O3"].reason


class TestGenieIndelResidues:
    def genome(self, tmp_path, seq):
        fasta = tmp_path / "g.fa"
        fasta.write_text(">chr1\n" + seq + "\n")
        (tmp_path / "g.fa.fai").write_text(f"chr1\t{len(seq)}\t6\t{len(seq)}\t{len(seq) + 1}\n")
        return svig.Genome(str(fasta))

    def test_right_align(self, tmp_path):
        g = self.genome(tmp_path, "ACGTTTTG")
        assert svig.right_align_maf(g, "chr1", 4, 4, "T", "-") == (7, 7, "T", "-")
        assert svig.right_align_maf(g, "chr1", 3, 4, "-", "TT") == (7, 8, "-", "TT")

    def test_left_aligned_insertion_counts_at_hgvs_residue(self, tmp_path):
        # "at the same position or further downstream" [S T1 O4]: a repeat-shifted insertion that
        # GENIE stores left-aligned must count at the residue HGVS reports (3'-shifted)
        tx = tx_resources(tmp_path).transcripts["NM_000001"][0]
        seq = ["C"] * 6400
        seq[2087:2093] = "AGAGAG"  # exon 2, 1-based 2088-2093 (residues 60-61)
        g = self.genome(tmp_path, "".join(seq))
        e = svig.GenieEntry("truncating", "chr1", 2087, 2088, "-", "AG", "", "", "Myeloid", 1)
        # left-aligned: first altered base 2088 = CDS nt 178 = residue 60; shifted 3' past the
        # AGAGAG repeat: base 2094 = CDS nt 184 = residue 62
        assert svig.entry_residue(tx, e, None) == 60
        assert svig.entry_residue(tx, e, g) == 62


class TestTier4Review:
    def test_review_status_and_flags(self, tmp_path):
        rows = [row(Hugo_Symbol="GENE1", gnomADe_AF="0.002", REVEL="0.9")]
        header, meta, records = TestTier4().run(tmp_path, rows)
        r = records[0]
        assert "clinical scientist" in r["SVIG_UK_REVIEW_STATUS"]
        flags = r["SVIG_UK_REVIEW_FLAGS"]
        assert "conflicting oncogenic and benign evidence" in flags      # B1 strong + O6
        assert "gene mode of action unknown" in flags
        assert "O11 and B7" in flags
        assert "professional" in meta[header.index("SVIG_UK_REVIEW_STATUS")]

    def test_boundary_flag(self):
        result = svig.Result("Likely Oncogenic", 6, [svig.onc("O2", "STR", "x"), svig.onc("O3", "MOD", "y")], [],
                             gene_role="tumour suppressor")
        assert "score within 1 point of a class boundary" in svig.review_flags(result)

    def test_json_export(self, tmp_path):
        maf = TestTier4().write_maf(tmp_path, [row(Hugo_Symbol="GENE1", gnomADe_AF="0.2")])
        out = tmp_path / "t4.json"
        svig.main(["--maf", str(maf), "--tier4-json", str(out), "--canonical", CANONICAL])
        data = json.loads(out.read_text())
        v = data["variants"][0]
        assert data["smart_version"] == "2.0.0" and "SVIG-UK" in data["classification_system"]
        assert v["classification"] == "Benign" and v["score"] is None
        assert v["codes"][0]["code"] == "B1" and v["codes"][0]["strength"] == "STANDALONE"
        assert "O11" in v["codes_not_automated"] and "O2" in v["codes_not_met"]


# ===========================================================================
# Gene-specific splice provisions and the B4 exon-edge caution
#   [G §7.11; S T2 B4 p.16] — synthetic 5-exon transcript: exon 2 spans CDS 91-192
#   (102 coding bp, in frame), exon 3 spans 193-292 (100 bp, out of frame)
# ===========================================================================

GSS_TABLE = os.path.join(SCRIPTS_DIR, "..", "reference_lists", "svig_uk_gene_specific_splice.tsv")


def gss_resources(tmp_path, rows=None, strand="+"):
    """Transcript resources with a gene-specific splice table (the shipped one by default)."""
    res = tx_resources(tmp_path, strand)
    if rows is None:
        res.load_gene_specific_splice(GSS_TABLE)
    else:
        path = tmp_path / "gss.tsv"
        header = "gene\ttranscript\texon\tcds_position\tapplies_to\tassessment\tsource\n"
        path.write_text(header + "".join("\t".join(map(str, r)) + "\n" for r in rows))
        res.load_gene_specific_splice(str(path))
    return res


class TestExonEdgePositions:
    """B4 "should be used with caution for synonymous variants ... located at the first base
    or last three bases of the exon" [S T2 B4 p.16]."""

    @pytest.mark.parametrize("cds_pos,expected", [
        (91, "first base of the exon"),
        (192, "last base of the exon"),
        (191, "second to last base of the exon"),
        (190, "third to last base of the exon"),
        (189, None),          # fourth to last: the caution does not reach it
        (150, None),          # middle of the exon
    ])
    def test_edge_detection(self, tmp_path, cds_pos, expected):
        tx = gss_resources(tmp_path).transcript_for(tx_row())
        assert svig.exon_edge_position(tx, cds_pos)[1] == expected

    def test_span_is_none_for_untranslated_exon(self, tmp_path):
        tx = gss_resources(tmp_path).transcript_for(tx_row())
        assert svig.exon_cds_span(tx, 1) == (1, 90)
        assert svig.exon_cds_span(tx, 99) is None


class TestB4ExonEdgeCaution:
    """B4 is withheld at the exon edges, where the risk of splice-motif disruption is higher."""

    def syn_row(self, cds_pos):
        return tx_row(Consequence="synonymous_variant", HGVSp="p.Thr50=",
                      CDS_position=str(cds_pos), EXON="2/5", **spliceai(0.0))

    def test_b4_withheld_at_the_last_base_of_an_exon(self, tmp_path):
        r = svig.classify_row(self.syn_row(192), gss_resources(tmp_path))
        assert "B4" not in codes(r)
        assert "B3" in codes(r)          # the computational code itself still applies
        assert any("B4 not applied" in n for n in r.notes)
        assert any("exon edge" in f for f in svig.review_flags(r))

    def test_b4_withheld_at_the_first_base_of_an_exon(self, tmp_path):
        r = svig.classify_row(self.syn_row(91), gss_resources(tmp_path))
        assert "B4" not in codes(r)

    def test_b4_applied_in_the_middle_of_the_exon(self, tmp_path):
        r = svig.classify_row(self.syn_row(150), gss_resources(tmp_path))
        assert codes(r)["B4"].strength == "SUPP"
        assert not any("B4 not applied" in n for n in r.notes)

    def test_caution_does_not_touch_intronic_or_utr_b4(self, tmp_path):
        r = svig.classify_row(tx_row(Consequence="intron_variant", HGVSp="", HGVSc="c.192+30A>G",
                                     CDS_position="", EXON="", **spliceai(0.0)),
                              gss_resources(tmp_path))
        assert codes(r)["B4"].strength == "SUPP"


class TestGeneSpecificSpliceProvision:
    """§7.11: synonymous changes at the last base of TP53 exons 4, 6 and 9 "should be classed as
    driver variants as they affect splicing of these exons"."""

    def test_shipped_table_matches_the_guideline(self):
        table = svig.load_gene_specific_splice(GSS_TABLE)
        entries = table[("TP53", "NM_000546")]
        # the last coding base of exons 4, 6 and 9 of NM_000546.6
        assert {(e["exon"], e["cds_position"]) for e in entries} == {(4, 375), (6, 672), (9, 993)}
        assert all(e["applies_to"] == "synonymous" and e["assessment"] == "driver" for e in entries)

    def test_lookup_requires_gene_transcript_position_and_class(self, tmp_path):
        res = gss_resources(tmp_path)
        hit = row(Hugo_Symbol="TP53", NM_Transcript="NM_000546.6", CDS_position="375")
        assert svig.gene_specific_splice(hit, res, "synonymous") is not None
        assert svig.gene_specific_splice(hit, res, "missense") is None          # applies_to
        assert svig.gene_specific_splice({**hit, "CDS_position": "374"}, res, "synonymous") is None
        assert svig.gene_specific_splice({**hit, "Hugo_Symbol": "TP63"}, res, "synonymous") is None
        assert svig.gene_specific_splice({**hit, "NM_Transcript": "NM_999999.1"}, res, "synonymous") is None

    def test_position_read_from_hgvsc_when_cds_position_is_absent(self, tmp_path):
        res = gss_resources(tmp_path)
        r = row(Hugo_Symbol="TP53", NM_Transcript="NM_000546.6", CDS_position="", HGVSc="c.375G>A")
        assert svig.gene_specific_splice(r, res, "synonymous") is not None

    def provision_row(self, cds_pos=192):
        """A synonymous variant at the position the synthetic table marks as splice-disrupting."""
        return tx_row(Consequence="synonymous_variant", HGVSp="p.Thr64=",
                      CDS_position=str(cds_pos), EXON="2/5", **spliceai(0.6))

    def synthetic(self, tmp_path):
        return gss_resources(tmp_path, rows=[
            ("GENE1", "NM_000001.1", 2, 192, "synonymous", "driver", "test provision")])

    def test_o2_strength_comes_from_exon_skipping(self, tmp_path):
        """Exon 2 is 102 coding bp: skipping it is in frame and removes 34 of 199 residues
        (17% > 10%), so the PVS1 tree gives strong rather than a fixed weight."""
        r = svig.classify_row(self.provision_row(), self.synthetic(tmp_path))
        o2 = codes(r)["O2"]
        assert o2.strength == "STR"
        assert "skipping of exon 2/5" in o2.reason and "preserves the reading frame" in o2.reason

    def test_benign_computational_codes_are_not_applied(self, tmp_path):
        r = svig.classify_row(self.provision_row(), self.synthetic(tmp_path))
        assert "B3" not in codes(r) and "B4" not in codes(r)

    def test_spliceai_route_to_o6_is_not_counted_twice(self, tmp_path):
        """O2's strength already rests on an assumed splicing effect, as for canonical sites."""
        r = svig.classify_row(self.provision_row(), self.synthetic(tmp_path))
        assert "O6" not in codes(r)

    def test_review_flag_and_note_are_raised(self, tmp_path):
        r = svig.classify_row(self.provision_row(), self.synthetic(tmp_path))
        assert any("gene-specific provision" in n for n in r.notes)
        assert any("driver" in n for n in r.notes)
        assert any("gene-specific guideline provision" in f for f in svig.review_flags(r))

    def test_lowest_strength_without_a_transcript_model(self, tmp_path):
        res = self.synthetic(tmp_path)
        res.transcripts = {}
        r = svig.classify_row(self.provision_row(), res)
        o2 = codes(r)["O2"]
        assert o2.strength == "MOD" and "lowest strength applied" in o2.reason

    def test_other_positions_in_the_same_gene_are_untouched(self, tmp_path):
        r = svig.classify_row(self.provision_row(cds_pos=150), self.synthetic(tmp_path))
        assert "O2" not in codes(r)
        assert not any("gene-specific provision" in n for n in r.notes)

    def test_provision_requires_the_gene_to_be_a_tumour_suppressor(self, tmp_path):
        """O2 is only evaluated for genes with a tumour-suppressor role [S T1 O2]."""
        res = self.synthetic(tmp_path)
        res.gene_roles = {"GENE1": {"ONC"}}
        r = svig.classify_row(self.provision_row(), res)
        assert "O2" not in codes(r)


# ===========================================================================
# O10 / B6 from functional assays [S T1 O10 p.13; S T2 B6 p.17]
#   "Variant is considered to have functionally abnormal effect if protein activity assay or
#    functional impact is <25% of wildtype level"; an intermediate result may be downgraded
#    "by one or more levels"; the strength itself follows ClinGen SVI / CanVIG-UK calibration.
# ===========================================================================

ASSAY_PARAMS = os.path.join(SCRIPTS_DIR, "..", "reference_lists", "svig_uk_functional_assays.tsv")
ASSAY_SCORES = ASSAY_PARAMS.replace("_assays.tsv", "_assay_scores.tsv")


def assay_resources(tmp_path, scores, strength="MOD", abnormal="25", intermediate="",
                    normal="", applies_to="missense", transcript="NM_000001.1", roles=None):
    """Resources with one synthetic assay covering GENE1. scores: {protein: pct_wt}."""
    res = make_resources(roles=roles or TSG)
    params = tmp_path / "assays.tsv"
    params.write_text(
        "assay_id\tgene\ttranscript\tapplies_to\tmetric\tabnormal_max_pct\t"
        "intermediate_max_pct\tnormal_min_pct\tstrength\tsource\tpmid\tlicence\turn\n"
        f"A1\tGENE1\t{transcript}\t{applies_to}\ttest activity\t{abnormal}\t{intermediate}\t"
        f"{normal}\t{strength}\ttest\t12345\tCC0\turn:test\n")
    scores_file = tmp_path / "assay_scores.tsv"
    scores_file.write_text("assay_id\tprotein\tpct_wt\n"
                           + "".join(f"A1\t{p}\t{v}\n" for p, v in scores.items()))
    res.load_functional_assays(str(params), str(scores_file))
    return res


def assay_row(**kw):
    base = {"NM_Transcript": "NM_000001.1", "HGVSp": "p.Arg100His"}
    base.update(kw)
    return row(**base)


class TestFunctionalAssayO10:
    def test_below_25_percent_is_functionally_abnormal(self, tmp_path):
        res = assay_resources(tmp_path, {"R100H": "14.4"})
        o10 = codes(svig.classify_row(assay_row(), res))["O10"]
        assert o10.strength == "MOD" and o10.points == 2
        assert "14.4% of wild-type activity" in o10.reason and "functionally abnormal" in o10.reason

    def test_at_the_threshold_is_not_abnormal(self, tmp_path):
        """<25%, so exactly 25% does not qualify."""
        res = assay_resources(tmp_path, {"R100H": "25.0"})
        assert "O10" not in codes(svig.classify_row(assay_row(), res))

    def test_strength_comes_from_the_assay_table(self, tmp_path):
        for strength, points in (("VSTR", 8), ("STR", 4), ("MOD", 2), ("SUPP", 1)):
            res = assay_resources(tmp_path, {"R100H": "5"}, strength=strength)
            o10 = codes(svig.classify_row(assay_row(), res))["O10"]
            assert (o10.strength, o10.points) == (strength, points)

    def test_intermediate_result_is_downgraded_one_level(self, tmp_path):
        res = assay_resources(tmp_path, {"R100H": "40"}, strength="STR", intermediate="50")
        o10 = codes(svig.classify_row(assay_row(), res))["O10"]
        assert o10.strength == "MOD" and "intermediate" in o10.reason
        assert "downgraded from STR" in o10.reason

    def test_no_intermediate_boundary_means_no_code(self, tmp_path):
        """The guideline fixes no numeric boundary for 'intermediate', so none is assumed."""
        res = assay_resources(tmp_path, {"R100H": "40"}, intermediate="")
        assert "O10" not in codes(svig.classify_row(assay_row(), res))

    def test_supporting_cannot_be_downgraded_further(self, tmp_path):
        res = assay_resources(tmp_path, {"R100H": "40"}, strength="SUPP", intermediate="50")
        assert "O10" not in codes(svig.classify_row(assay_row(), res))

    def test_loss_of_function_needs_a_tumour_suppressor_role(self, tmp_path):
        """A null result supports oncogenicity only where loss of function is the mechanism."""
        res = assay_resources(tmp_path, {"R100H": "5"}, roles=ONCOGENE)
        assert "O10" not in codes(svig.classify_row(assay_row(), res))

    def test_b6_only_when_a_normal_boundary_is_recorded(self, tmp_path):
        res = assay_resources(tmp_path, {"R100H": "95"}, normal="")
        assert "B6" not in codes(svig.classify_row(assay_row(), res))
        res = assay_resources(tmp_path, {"R100H": "95"}, normal="75", strength="STR")
        b6 = codes(svig.classify_row(assay_row(), res))["B6"]
        assert b6.strength == "STR" and b6.points == -4 and "no damaging effect" in b6.reason

    def test_other_transcripts_and_classes_are_not_matched(self, tmp_path):
        """"Evidence should refer to biologically relevant transcripts"; applies_to is respected."""
        res = assay_resources(tmp_path, {"R100H": "5"})
        assert "O10" not in codes(svig.classify_row(assay_row(NM_Transcript="NM_999999.9"), res))
        res = assay_resources(tmp_path, {"R100H": "5"}, applies_to="inframe")
        assert "O10" not in codes(svig.classify_row(assay_row(), res))

    def test_variant_not_in_the_assay_falls_back_to_the_oncokb_proxy(self, tmp_path):
        res = assay_resources(tmp_path, {"R200H": "5"})       # a different residue
        r = svig.classify_row(assay_row(VARIANT_IN_ONCOKB="True",
                                        MUTATION_EFFECT="Loss-of-function"), res)
        o10 = codes(r)["O10"]
        assert o10.strength == "SUPP" and "proxy" in o10.reason

    def test_assay_takes_precedence_over_the_proxy(self, tmp_path):
        res = assay_resources(tmp_path, {"R100H": "5"}, strength="STR")
        r = svig.classify_row(assay_row(VARIANT_IN_ONCOKB="True",
                                        MUTATION_EFFECT="Loss-of-function"), res)
        o10 = codes(r)["O10"]
        assert o10.strength == "STR" and "proxy" not in o10.reason
        assert any("functional assay" in n for n in r.notes)
        assert any("functional assay" in f for f in svig.review_flags(r))


class TestShippedFunctionalAssayTable:
    def test_tp53_transactivation_assay_is_present_and_sourced(self):
        res = svig.Resources()
        res.load_functional_assays(ASSAY_PARAMS, ASSAY_SCORES)
        assays = res.functional_assays["TP53"]
        assert len(assays) == 1
        a = assays[0]
        assert a["abnormal_max_pct"] == "25"        # the guideline's own threshold
        assert a["strength"] in ("VSTR", "STR", "MOD", "SUPP")
        assert a["pmid"] and a["licence"] == "CC0" and a["transcript"].startswith("NM_000546")
        # the boundaries the guideline does not fix are left empty on purpose
        assert a["intermediate_max_pct"] == "" and a["normal_min_pct"] == ""

    def test_scores_are_percentages_of_wild_type(self):
        res = svig.Resources()
        res.load_functional_assays(ASSAY_PARAMS, ASSAY_SCORES)
        values = [v for (assay, _), v in res.functional_scores.items()
                  if assay == "TP53_transactivation_median"]
        assert len(values) > 2000
        assert all(0 <= v <= 250 for v in values)

    def test_a_known_null_hotspot_is_called_functionally_abnormal(self):
        """R175H and R273H are archetypal TP53 null variants; both must be below 25% of wild type."""
        res = svig.Resources()
        res.load_functional_assays(ASSAY_PARAMS, ASSAY_SCORES)
        for protein in ("R175H", "R273H", "R248Q"):
            pct = res.functional_scores[("TP53_transactivation_median", protein)]
            assert pct < 25, f"{protein} measured at {pct}% of wild type"
