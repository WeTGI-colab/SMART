#!/usr/bin/env python3

# =============================================================================
# svig_uk_classifier.py
# =============================================================================
#
# DESCRIPTION:
#   Automated oncogenicity classification of somatic small variants following
#   the ACGS / SVIG-UK 2025 guidelines ("Association for Clinical Genomic Science
#   guidelines for the classification of oncogenicity of somatic variants in
#   cancer", ratified 22/07/2025; version 1.1, 26/05/2026) and their Supplementary
#   Material. Figure and combination numbers follow v1.1 (Supplementary Table 6).
#
#   The scorer reads a SMART MAF (Final_result_tier1.maf) and appends four
#   columns: SVIG_UK_CLASS, SVIG_UK_SCORE, SVIG_UK_CODES and SVIG_UK_NOTES.
#
#   DESIGN PRINCIPLES
#     1. Every rule follows the guideline text. Citations use the form
#        [G §x] for the main guideline and [S Tn] / [S Fn] for Supplementary
#        Tables / Figures.
#     2. Only the official strengths are used: +8/+4/+2/+1 and -8/-4/-2/-1
#        [G Fig 1B], plus the stand-alone outcomes of O1, B1 and B2.
#     3. Where the guideline requires a judgement that cannot be taken from
#        annotation data, the rule is either omitted or replaced by a clearly
#        labelled "SMART proxy" — an automatable approximation chosen to be
#        conservative (never stronger than the guideline would allow).
#     4. Each code is applied once, the code-combination restrictions of
#        [S F1] and [S F2] are enforced, and the minimum-two-codes rule of
#        [G §6.3] is applied.
#     5. The weights were fixed from the guideline alone, before any evaluation
#        against expert-classified variants.
#
#   CODES
#     Implemented as written : O1, O3, O6, O7 (Cancer Hotspots route), O8,
#                              B1, B3, B4 (in silico route)
#     Implemented with proxy : O2 (PVS1 tree from RefSeq exon structure; region
#                              criticality not judged), O4 (GENIE; COSMIC
#                              all-tissue counts as fallback), O9 and B5 (VEP domains),
#                              B2 (null variants in oncogenes only),
#                              O10 / B6 (OncoKB variant-level mutation effect)
#                              O5 missense route (canonical list + catalogue
#                              references, REVEL from the local REVEL table)
#     Not automatable        : O5 splice / initiation-codon routes, O11 / B7 (tumour phenotype),
#                              O2(RNA) / B4(RNA) (splicing assays)
#
# INPUT:
#   --maf        SMART Final_result_tier1.maf
#   --refs       reference directory (defaults below are resolved against it)
#   --tumour-type used when the MAF has no TUMOUR_TYPE column
#                (myeloid | lymphoid | haem | solid | unknown)
#
# OUTPUT:
#   --out        the input MAF with the four SVIG_UK_* columns appended
# =============================================================================

import argparse
import csv
import gzip
import json
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field

# =============================================================================
# Guideline constants
# =============================================================================

POINTS = {"VSTR": 8, "STR": 4, "MOD": 2, "SUPP": 1}  # [G Fig 1B]

CLASS_THRESHOLDS = [  # [G Fig 1C]
    (10, "Oncogenic"),
    (6, "Likely Oncogenic"),
    (0, "VUS"),
    (-6, "Likely Benign"),
]
ONC_CLASSES = {"Oncogenic", "Likely Oncogenic"}
BEN_CLASSES = {"Likely Benign", "Benign"}

# gnomAD v4 sample sizes per genetic-ancestry group. Genomes = gnomAD v3 column
# and exomes = v4 total minus v3 genomes, from the gnomAD v4.0 release table
# (gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0). The FIN/NFE split and
# the genome AMI/remaining split are approximate. These sizes are only used to
# (a) identify sub-populations with >=100,000 alleles and (b) estimate allele
# counts from allele frequencies for the O3 ">1 allele" rule [S T1 O3].
GNOMAD_SAMPLES = {
    ("gnomADe", "AFR"): 16801, ("gnomADe", "AMR"): 22372, ("gnomADe", "ASJ"): 13068,
    ("gnomADe", "EAS"): 19844, ("gnomADe", "FIN"): 26710, ("gnomADe", "NFE"): 556002,
    ("gnomADe", "MID"): 2873, ("gnomADe", "SAS"): 43127, ("gnomADe", "REMAINING"): 29669,
    ("gnomADg", "AFR"): 20744, ("gnomADg", "AMR"): 7647, ("gnomADg", "ASJ"): 1736,
    ("gnomADg", "EAS"): 2604, ("gnomADg", "FIN"): 5316, ("gnomADg", "NFE"): 34029,
    ("gnomADg", "MID"): 158, ("gnomADg", "SAS"): 2419, ("gnomADg", "AMI"): 456,
    ("gnomADg", "REMAINING"): 1047,
}
B1_POPULATIONS = ("AFR", "EAS", "NFE", "AMR", "SAS")  # [S T2 B1] five continental groups

# VEP DOMAINS sources treated as "known functional domain" (O9) — curated domain
# and motif databases. Structural mappings (PDB, AFDB), whole-family (PANTHER),
# structural classification (Gene3D, Superfamily) and sequence-feature
# predictors (low complexity, coiled coils, disorder, signal peptides,
# transmembrane) are excluded. SMART proxy.
FUNCTIONAL_DOMAIN_DBS = {"Pfam", "SMART", "PROSITE_profiles", "PROSITE_patterns",
                         "CDD", "Prints", "HAMAP", "NCBIFAM", "PIRSF"}
REPEAT_DBS = {"Low_complexity_(Seg)"}

# Genes whose mode of action depends on tumour type and is stated explicitly in
# the guideline: "EZH2 acts as a tumour suppressor gene in myeloid disorders
# (e.g. MDS/AML), whereas in lymphoid disorders (e.g. B-cell lymphoma), EZH2
# variants are activating" [G §7.2].
CONTEXT_ROLE_OVERRIDES = {
    ("EZH2", "myeloid"): {"TSG"},
    ("EZH2", "lymphoid"): {"ONC"},
}

# O4 case-counting thresholds [S T1 O4]. The guideline sets them for GENIE
# (missense >10 strong, 5-10 moderate; truncating and in-frame >50 / 20-50 / 10-19)
# and requires "much higher thresholds" when COSMIC is used, giving ">50 solid
# tumour entries for missense variants to achieve a strong strength" as the
# example. SMART proxy: COSMIC Cancer Mutation Census sample counts across all
# tissues with every GENIE threshold x5, which reproduces that example (>10 -> >50).
# On-target (tumour-type) counting is not possible: the CMC tissue breakdown is
# empty for ~99.5% of mutations. Entries from other tumour types may be combined
# where the mechanism of action is consistent [S T1 O4 guidance].
O4_THRESHOLDS = {
    # class: (strong_gt, moderate_ge, supporting_ge or None)
    "missense": (50, 25, None),
    "truncating": (250, 100, 50),
    "inframe": (250, 100, 50),
}

# O4 from AACR Project GENIE, thresholds as written in [S T1 O4] for an international
# curated database: missense and splice >10 strong, 5-10 moderate; frameshift/nonsense
# (same position or downstream) and in-frame (exact insertion; deletion equal to or
# within the deleted region) >50 strong, 20-50 moderate, 10-19 supporting.
GENIE_O4_THRESHOLDS = {
    "missense": (10, 5, None),
    "splice": (10, 5, None),
    "truncating": (50, 20, 10),
    "inframe": (50, 20, 10),
}
# 'On-target' tumour types [S T1 O4]: "when assessing a variant detected in a patient with
# acute myeloid leukaemia (AML), 'on-target' entries may be specific to AML, to myeloid
# disorders or to all haematopoietic disorders (myeloid and lymphoid). The definition of
# an 'on-target' tumour type is currently at the discretion of individual laboratories".
# The tumour type is an OncoTree code (e.g. AML) or a broad group (myeloid, lymphoid,
# haem, solid). --o4-on-target chooses the level counted:
#   code    the OncoTree code and its subtypes (AML -> AML, AMLMRC, AMLNOS, ...)
#   tissue  the OncoTree tissue (AML -> all Myeloid)
#   lineage haematopoietic (myeloid + lymphoid) or solid tumours (default)
#   all     every tumour type ("off-target" entries included)
# A level finer than the tumour type given falls back to the finest one available; an
# unknown tumour type counts all entries (the criterion "has been designed to be tumour
# agnostic"). Genes whose mode of action differs by lineage (CONTEXT_ROLE_OVERRIDES) are
# never counted beyond their own tissue.
O4_ON_TARGET_LEVELS = ("code", "tissue", "lineage", "all")
HAEM_TISSUES = {"Myeloid", "Lymphoid"}

# O5 missense route [S T1 O5]: "Novel missense change at an amino acid residue where a
# different amino acid change determined to be (likely) oncogenic has been seen before":
# moderate when the variant's REVEL is equivalent (difference <0.02) or higher than the
# reference variant's, or both are >0.773; supporting when the variant's REVEL is >0.7
# and <0.773 and lower than the reference's.
O5_REVEL_EQUIVALENT = 0.02
O5_REVEL_BOTH_HIGH = 0.773
O5_REVEL_SUPPORTING_MIN = 0.7

# OncoKB variant-level mutation effect -> O10 / B6 (SMART proxy). SVIG-UK O10/B6
# require well-established functional studies assessed with the Brnich et al.
# framework [S T1 O10], which cannot be done automatically. OncoKB's curated
# mutation effect summarises published functional evidence, so a definitive
# effect consistent with the gene mechanism is used at SUPPORTING strength only
# — the lowest level — reflecting that "the stringency of current guidance ...
# frequently prevents the application of O10 above a supporting level" [G §7.14].
ONCOKB_GOF = {"Gain-of-function", "Switch-of-function"}
ONCOKB_LOF = {"Loss-of-function"}
ONCOKB_NEUTRAL = {"Neutral"}

# O3 CHIP correction [S T1 O3]: "Somatic mosaicism of variants in genes such as (but
# not limited to) DNMT3A, ASXL1, TET2 and TP53 during hematopoietic clonal expansion
# (CHIP) ... If the variant under investigation does not meet any population rule codes
# (O3 or B1) AND has a total allele frequency >0.00001 (0.001%) then low variant allele
# frequency (VAF) alleles ... can be excluded and the total allele frequency based on the
# number of alleles with VAF between 0.3 and 0.65 can be re-calculated."
CHIP_GENES = {"DNMT3A", "ASXL1", "TET2", "TP53"}

# Gene-specific population thresholds for hereditary cancer genes [S T1 O3; S T2 B1]: "if the
# somatic variant is in a gene known to cause predisposition to hereditary cancer, ClinGen VCEP
# gene-specific guidelines must be consulted (if they exist) to establish a cut-off"; B1 applies at
# ">1% (or as specified by respective expert panel group)". PM2 -> O3, BA1 -> B1 stand-alone,
# BS1 -> B1 strong (reference_lists/svig_uk_vcep_population_thresholds.tsv, from ClinGen CSpec).
# SMART proxy: VCEPs mostly specify the GrpMax filtering allele frequency (FAF95); the annotation
# carries allele frequencies only, so allele counts are estimated from the gnomAD v4 group sizes and
# the FAF95 is computed from them. Founder-effect groups are ignored as the VCEPs specify.
CONTINENTAL_GROUPS = ("AFR", "AMR", "EAS", "NFE", "SAS")
VCEP_MIN_GROUP_ALLELES = 2000
CHIP_MIN_TOTAL_AF = 1e-5
# gnomAD ab_hist_alt_bin_freq bin edges 0.00|0.05|...|1.00: bins 6-12 = VAF 0.30-0.65
CHIP_AB_BINS = slice(6, 13)

# O2 / PVS1 decision tree constants [S T1 O2; Abou Tayoun et al. 2018].
# "if the premature termination codon (PTC) occurs in the last 50 bp's of the
# penultimate exon or in the final exon, then NMD would generally not be expected"
NMD_PENULTIMATE_EXON_BP = 50
# "is the missing region critical to normal function/does it account for >10% of
# the entire protein?" — region criticality cannot be judged automatically, so
# only the >10% branch is used.
PVS1_PROTEIN_FRACTION = 0.10
# "look for nearby (suggest up to approximately +/- 20 bp) cryptic splice sites or
# newly generated splice sites using splice predictive tools (e.g. SpliceAI)".
# A gain is counted at the SpliceAI delta >= 0.2 used for O6 [S T1 O6].
CRYPTIC_SITE_WINDOW_BP = 20
CRYPTIC_SITE_MIN_DS = 0.2
STOP_CODONS = {"TAA", "TAG", "TGA"}
GENETIC_CODE = {a + b + c: aa for (a, b, c), aa in zip(
    ((a, b, c) for a in "TCAG" for b in "TCAG" for c in "TCAG"),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}

AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
    "Sec": "U", "Pyl": "O", "Xaa": "X",
}

# =============================================================================
# Parsing helpers
# =============================================================================


def fnum(v):
    if v is None:
        return None
    s = str(v).strip()
    if s in ("", ".", "-", "nan", "NaN", "None", "NA"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def first_int(v):
    m = re.search(r"\d+", str(v or ""))
    return int(m.group()) if m else None


def normalise_protein(hgvsp):
    """Return a canonical one-letter protein change without prefix or
    frameshift tail, e.g. 'ENSP..:p.Trp288CysfsTer12' -> 'W288Cfs',
    'p.(Arg201Ter)' -> 'R201*', 'p.Pro95_Arg102del' -> 'P95_R102del',
    'p.Thr477%3D' -> 'T477='. Returns '' when no protein change is available."""
    s = str(hgvsp or "").strip()
    if s in ("", ".", "nan", "None"):
        return ""
    s = s.split(":")[-1].replace("%3D", "=").replace("(", "").replace(")", "")
    if s.startswith("p."):
        s = s[2:]
    if s in ("", "?", "0", "="):
        return ""
    for three, one in AA3.items():
        s = s.replace(three, one)
    m = re.match(r"^(.*?fs)", s)
    return m.group(1) if m else s


PROT_RE = re.compile(r"^([A-Z*])(\d+)(?:_([A-Z*])(\d+))?(.*)$")


def protein_parts(protein):
    """'R882H' -> ('R', 882, 882, 'H'); 'P95_R102del' -> ('P', 95, 102, 'del')."""
    m = PROT_RE.match(protein or "")
    if not m:
        return None
    start = int(m.group(2))
    end = int(m.group(4)) if m.group(4) else start
    return m.group(1), start, end, m.group(5)


def terms(consequence):
    return set(str(consequence or "").split("&"))


def variant_class(consequence):
    t = terms(consequence)
    if t & {"stop_gained", "frameshift_variant"}:
        return "truncating"
    if t & {"splice_donor_variant", "splice_acceptor_variant"}:
        return "canonical_splice"
    if "start_lost" in t:
        return "start_lost"
    if "stop_lost" in t:
        return "stop_lost"
    if t & {"inframe_insertion", "inframe_deletion", "protein_altering_variant"}:
        return "inframe"
    if "missense_variant" in t:
        return "missense"
    if t & {"synonymous_variant", "stop_retained_variant", "start_retained_variant"}:
        return "synonymous"
    if t & {"5_prime_UTR_variant", "3_prime_UTR_variant"}:
        return "utr"
    if "intron_variant" in t:
        return "intronic"
    return "other"


NULL_CLASSES = {"truncating", "canonical_splice", "start_lost", "donor_last_base"}
# classes whose strength already rests on an assumed splicing effect: the SpliceAI
# route to O6 would count the same evidence twice
SPLICE_ASSUMED_CLASSES = {"canonical_splice", "donor_last_base"}
SPLICE_REGION_TERMS = {"splice_region_variant", "splice_donor_5th_base_variant",
                       "splice_donor_region_variant", "splice_polypyrimidine_tract_variant"}


def exon_fraction(v):
    m = re.match(r"^\s*(\d+)(?:-\d+)?/(\d+)", str(v or ""))
    return (int(m.group(1)), int(m.group(2))) if m else (None, None)


def intronic_offset(hgvsc):
    m = re.search(r"c\.[*-]?\d+([+-]\d+)", str(hgvsc or ""))
    return int(m.group(1)) if m else None


def touches_canonical_splice_site(hgvsc):
    """True when the HGVS c. description (3'-normalised) includes a +/-1 or +/-2
    intronic position or spans an exon-intron boundary. VEP's splice_donor/acceptor
    consequence is computed on the left-aligned allele, so a repeat-shifted deletion
    such as c.100+5_100+8del in a GTTT repeat can be flagged although the donor dinucleotide is intact."""
    m = re.search(r"c\.([^A-Za-z>]*)", str(hgvsc or ""))
    if not m:
        return None
    tokens = re.findall(r"[*-]?\d+([+-]\d+)?", m.group(1))
    if not tokens:
        return None
    offsets = [int(t) for t in tokens if t]
    if not offsets:
        return False
    if len(offsets) < len(tokens):
        return True  # one end exonic, the other intronic: the boundary is included
    lo, hi = min(offsets), max(offsets)
    return any(lo <= o <= hi for o in (1, 2, -1, -2))


def spliceai_scores(row):
    vals = [fnum(row.get(c)) for c in ("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
                                       "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")]
    vals = [v for v in vals if v is not None]
    return vals or None


def ptc_residues(hgvsp):
    """Return (first altered residue, residue of the new stop codon) for a truncating
    protein change; the stop is None when HGVS gives no in-frame stop ('fs*?').
    'p.Arg100Ter' -> (100, 100); 'p.Phe126LeufsTer5' -> (126, 130) because the
    frameshift stop is counted from the first changed residue (= position 1)."""
    s = str(hgvsp or "").split(":")[-1].replace("(", "").replace(")", "")
    if s.startswith("p."):
        s = s[2:]
    for three, one in AA3.items():
        s = s.replace(three, one)
    m = re.match(r"^[A-Z](\d+)(?:_[A-Z]\d+)?[A-Z]*fs\*(\d+|\?)$", s)
    if m:
        first = int(m.group(1))
        return first, (None if m.group(2) == "?" else first + int(m.group(2)) - 1)
    m = re.match(r"^[A-Z](\d+)\*$", s)
    if m:
        return int(m.group(1)), int(m.group(1))
    m = re.match(r"^[A-Z](\d+)(?:_[A-Z](\d+))?(delins|ins|dup)([A-Z]*)\*", s)
    if m:
        start = int(m.group(1))
        end = int(m.group(2)) if m.group(2) else start
        stop = (end + 1 if m.group(3) == "ins" else start) + len(m.group(4))
        return start, stop
    return None


def indel_length_change(row):
    ref = str(row.get("Reference_Allele") or "").replace("-", "")
    alt = str(row.get("Tumor_Seq_Allele2") or "").replace("-", "")
    return len(alt) - len(ref)


def domain_sources(row):
    return {d.split(":")[0] for d in str(row.get("DOMAINS") or "").split("&") if d}


def tumour_group(tumour_type):
    t = str(tumour_type or "").strip().lower()
    if not t:
        return "unknown"
    if "myelo" in t or t in {"aml", "mds", "mpn", "cmml", "cml"}:
        return "myeloid"
    if "lympho" in t or t in {"all", "cll", "dlbcl"}:
        return "lymphoid"
    if t.startswith("haem") or t.startswith("hem") or "leuk" in t:
        return "haem"
    if "solid" in t:
        return "solid"
    return "unknown"


# =============================================================================
# Resources
# =============================================================================


def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _strip_version(tx):
    return str(tx or "").split(".")[0]


def _revcomp(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


class TranscriptModel:
    """Exon structure of one RefSeq transcript (UCSC ncbiRefSeq table: 0-based,
    half-open genomic coordinates). Exons are stored in transcript order."""

    def __init__(self, name, chrom, strand, starts, ends, cds_start, cds_end):
        exons = sorted(zip(starts, ends))
        self.name, self.chrom, self.strand = name, chrom, strand
        self.exons = exons[::-1] if strand == "-" else exons
        self.cds_start, self.cds_end = cds_start, cds_end
        self.tx_start, self.tx_end = exons[0][0], exons[-1][1]
        self.exon_lengths = [e - s for s, e in self.exons]
        self.coding = [max(0, min(e, cds_end) - max(s, cds_start)) for s, e in self.exons]
        self.cds_length = sum(self.coding)
        self.protein_length = self.cds_length // 3 - 1  # excluding the stop codon
        utr5 = 0
        for (s, e), c in zip(self.exons, self.coding):
            if c:
                utr5 += (e - min(e, cds_end)) if strand == "-" else (max(s, cds_start) - s)
                break
            utr5 += e - s
        self.utr5 = utr5

    @property
    def n_exons(self):
        return len(self.exons)

    def exon_has_start(self, k):
        s, e = self.exons[k - 1]
        pos = self.cds_end - 1 if self.strand == "-" else self.cds_start
        return s <= pos < e

    def exon_has_stop(self, k):
        s, e = self.exons[k - 1]
        pos = self.cds_start if self.strand == "-" else self.cds_end - 1
        return s <= pos < e

    def genomic_to_cds(self, pos):
        """1-based genomic position -> 1-based CDS nucleotide, or None if not coding."""
        p, cum = pos - 1, 0
        for s, e in self.exons:
            cs, ce = max(s, self.cds_start), min(e, self.cds_end)
            if cs >= ce:
                continue
            if cs <= p < ce:
                return cum + (p - cs if self.strand == "+" else ce - 1 - p) + 1
            cum += ce - cs
        return None

    def cds_to_genomic(self, c):
        """1-based CDS nucleotide -> 1-based genomic position, or None."""
        cum = 0
        for s, e in self.exons:
            cs, ce = max(s, self.cds_start), min(e, self.cds_end)
            if cs >= ce:
                continue
            if c <= cum + (ce - cs):
                k = c - cum - 1
                return (cs + k if self.strand == "+" else ce - 1 - k) + 1
            cum += ce - cs
        return None

    def residue_at(self, pos):
        c = self.genomic_to_cds(pos)
        return (c - 1) // 3 + 1 if c else None

    def exon_sequence(self, genome, k):
        s, e = self.exons[k - 1]
        seq = genome.fetch(self.chrom, s, e)
        return _revcomp(seq) if self.strand == "-" else seq


def locate_in_exons(exon_lengths, mrna_pos):
    """1-based mRNA position -> (exon number, nucleotides to the exon's 3' end)."""
    cum = 0
    for i, length in enumerate(exon_lengths, 1):
        if mrna_pos <= cum + length:
            return i, cum + length - mrna_pos
        cum += length
    return None


def nmd_predicted(exon_lengths, mrna_pos):
    """NMD rule quoted at NMD_PENULTIMATE_EXON_BP. Returns (predicted, where)."""
    loc = locate_in_exons(exon_lengths, mrna_pos)
    n = len(exon_lengths)
    if loc is None:
        return False, "beyond the transcript end"
    exon, to_end = loc
    if exon == n:
        return False, f"in the final exon ({exon}/{n})"
    if exon == n - 1 and to_end < NMD_PENULTIMATE_EXON_BP:
        return False, f"in the last {NMD_PENULTIMATE_EXON_BP} bp of the penultimate exon ({exon}/{n})"
    return True, f"in exon {exon}/{n}"


class Genome:
    """Minimal random-access reader for an uncompressed, faidx-indexed FASTA
    (the hg38.fa already required by SMART's VEP step)."""

    def __init__(self, fasta):
        self.index = {}
        with open(fasta + ".fai") as fh:
            for line in fh:
                f = line.split("\t")
                self.index[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
        self.fh = open(fasta, "rb")

    def fetch(self, chrom, start, end):
        length, offset, line_bases, line_width = self.index[chrom]
        start, end = max(0, start), min(end, length)

        def off(p):
            return offset + (p // line_bases) * line_width + p % line_bases

        self.fh.seek(off(start))
        raw = self.fh.read(off(end) - off(start))
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode().upper()


def left_align_maf(genome, chrom, start, end, ref, alt):
    """Left-align a MAF-style pure insertion (ref '-', start/end = flanking bases) or
    deletion (alt '-', start..end = deleted bases) on the reference genome."""
    if genome is None or (ref != "-" and alt != "-"):
        return start, end, ref, alt
    seq = (alt if ref == "-" else ref).upper()
    anchor = start if ref == "-" else start - 1  # last reference base before the event
    while anchor >= 1:
        base = genome.fetch(chrom, anchor - 1, anchor)
        if base != seq[-1]:
            break
        seq = base + seq[:-1]
        start, end, anchor = start - 1, end - 1, anchor - 1
    return (start, end, "-", seq) if ref == "-" else (start, end, seq, "-")


def right_align_maf(genome, chrom, start, end, ref, alt):
    """Shift a MAF-style pure insertion or deletion to its right-most (3' on the plus strand)
    equivalent position, as HGVS does."""
    if genome is None or (ref != "-" and alt != "-"):
        return start, end, ref, alt
    seq = (alt if ref == "-" else ref).upper()
    nxt = end if ref == "-" else end + 1  # first reference base after the event
    while True:
        base = genome.fetch(chrom, nxt - 1, nxt)
        if not base or base != seq[0]:
            break
        seq = seq[1:] + base
        start, end, nxt = start + 1, end + 1, nxt + 1
    return (start, end, "-", seq) if ref == "-" else (start, end, seq, "-")


def entry_residue(tx, e, genome=None):
    """Residue of the first transcript base a GENIE entry alters. Indels are shifted 3' in
    transcript direction first (HGVS convention): GENIE entries are left-aligned, which on the plus
    strand places a repeat-shifted insertion or deletion one residue early (NPM1 c.860_863dup)."""
    start, end, ref, alt = e.start, e.end, e.ref, e.alt
    if tx.strand == "+" and (ref == "-" or alt == "-"):
        start, end, ref, alt = right_align_maf(genome, e.chrom, start, end, ref, alt)
    if tx.strand == "+":
        return tx.residue_at(end if ref == "-" else start)
    return tx.residue_at(start if ref == "-" else end)


def maf_allele(row, genome=None):
    """The row's variant as (chrom, start, end, ref, alt) in left-aligned MAF style."""
    chrom = str(row.get("Chromosome") or "")
    chrom = chrom if chrom.startswith("chr") else "chr" + chrom
    pos = first_int(row.get("POS") or row.get("Start_Position"))
    ref = str(row.get("REF") or row.get("Reference_Allele") or "").replace("-", "").upper()
    alt = str(row.get("ALT") or row.get("Tumor_Seq_Allele2") or "").replace("-", "").upper()
    if pos is None or (not ref and not alt):
        return None
    while ref and alt and ref[0] == alt[0] and (len(ref) > 1 or len(alt) > 1):
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    if not ref:
        start, end, ref = pos - 1, pos, "-"
    elif not alt:
        start, end, alt = pos, pos + len(ref) - 1, "-"
    else:
        start, end = pos, pos + len(ref) - 1
    start, end, ref, alt = left_align_maf(genome, chrom, start, end, ref, alt)
    return chrom, start, end, ref, alt


@dataclass
class GenieEntry:
    vclass: str
    chrom: str
    start: int
    end: int
    ref: str
    alt: str
    ref_aa: str
    alt_aa: str
    lineage: str
    n: int
    code: str = ""
    tissue: str = ""

    def __post_init__(self):
        # one patient may carry the variant in samples of several tumour types ('|'-separated)
        lineages = self.lineage.split("|")
        codes = self.code.split("|") if self.code else [""] * len(lineages)
        tissues = self.tissue.split("|") if self.code else [""] * len(lineages)
        self.types = list(zip(codes, tissues, lineages))


class GenieO4:
    """Unique-patient GENIE counts built by utils/build_genie_o4.py (GRCh38, left-aligned)."""

    def __init__(self, path, genes=None):
        self.by_gene = defaultdict(list)
        with _open_text(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if genes is not None and r["gene"] not in genes:
                    continue
                self.by_gene[r["gene"]].append(GenieEntry(
                    r["variant_class"], r["chrom"], int(r["start"]), int(r["end"]), r["ref"], r["alt"],
                    r["ref_aa"], r["alt_aa"], r["lineage"], int(r["n_patients"]),
                    r.get("oncotree_code", ""), r.get("tissue", "")))


class OncoTree:
    """OncoTree release (JSON with tumorTypes: code, tissue, parent)."""

    def __init__(self, path):
        with open(path) as fh:
            types = json.load(fh)["tumorTypes"]
        self.tissue = {t["code"]: t["tissue"] for t in types}
        self.children = defaultdict(set)
        for t in types:
            self.children[t.get("parent")].add(t["code"])

    def subtree(self, code):
        out, stack = set(), [code]
        while stack:
            c = stack.pop()
            if c not in out:
                out.add(c)
                stack.extend(self.children.get(c, ()))
        return out


def resolve_tumour_type(value, oncotree=None):
    """Tumour type -> (OncoTree code or None, tissue or None, lineage set or None)."""
    v = str(value or "").strip()
    if oncotree is not None and v.upper() in oncotree.tissue:
        code = v.upper()
        tissue = oncotree.tissue[code]
        return code, tissue, (set(HAEM_TISSUES) if tissue in HAEM_TISSUES else {"Solid"})
    group = tumour_group(v)
    if group in ("myeloid", "lymphoid"):
        return None, group.capitalize(), set(HAEM_TISSUES)
    if group == "haem":
        return None, None, set(HAEM_TISSUES)
    if group == "solid":
        return None, None, {"Solid"}
    return None, None, None


def genie_target(value, res, gene):
    """Filter for on-target GENIE entries and a label describing it."""
    code, tissue, lineages = resolve_tumour_type(value, res.oncotree)
    level = res.o4_level
    if gene in {g for g, _ in CONTEXT_ROLE_OVERRIDES} and tissue in HAEM_TISSUES \
            and level in ("lineage", "all"):
        level = "tissue"  # lineage-dependent mode of action: never count the other lineage
    # a patient counts once if any of their samples with the variant is on-target
    if level == "code" and code and res.oncotree is not None:
        codes = res.oncotree.subtree(code)
        return (lambda e: any(c in codes for c, _, _ in e.types)), f"{code} and subtypes"
    if level in ("code", "tissue") and tissue:
        return (lambda e: any(t == tissue or (not t and l == tissue) for _, t, l in e.types)), f"tissue {tissue}"
    if level in ("code", "tissue", "lineage") and lineages:
        return (lambda e: any(l in lineages for _, _, l in e.types)), "+".join(sorted(lineages))
    return (lambda e: True), "all tumour types"


class RevelTabix:
    """REVEL scores by GRCh38 position from the tabix-indexed REVEL table
    (columns chr, hg19_pos, grch38_pos, ref, alt, aaref, aaalt, REVEL, ...)."""

    def __init__(self, path):
        self.path = path
        self.cache = {}

    def at(self, chrom, pos):
        key = (chrom, pos)
        if key not in self.cache:
            import subprocess
            region = f"{str(chrom).replace('chr', '')}:{pos}-{pos}"
            try:
                out = subprocess.run(["tabix", self.path, region], capture_output=True,
                                     text=True, check=True).stdout
            except (OSError, subprocess.CalledProcessError):
                out = ""
            rows = []
            for line in out.splitlines():
                f = line.split("\t")
                if len(f) >= 8 and f[2] == str(pos) and fnum(f[7]) is not None:
                    rows.append((f[3], f[4], f[5], f[6], float(f[7])))
            self.cache[key] = rows
        return self.cache[key]


class Resources:
    """Lookup tables. Every component is optional; a missing component means the
    dependent codes are reported as not evaluated rather than silently absent."""

    def __init__(self):
        self.canonical = []            # list of dicts from the O1 list
        self.gene_roles = {}           # gene -> set of {"ONC","TSG"}
        self.gene_roles_source = ""
        self.hotspots_single = {}      # (gene, aa_pos) -> entry
        self.hotspots_indel = defaultdict(list)
        self.mis_z = {}                # gene -> {transcript: z, "_mane": z}
        self.cosmic = None             # CosmicO4 or None
        self.genie = None              # GenieO4 or None
        self.o5_references = defaultdict(list)  # gene -> [(transcript, (ref, pos, alt), source)]
        self.revel = None              # RevelTabix or None
        self.o4_level = "lineage"      # GENIE on-target level, see O4_ON_TARGET_LEVELS
        self.oncotree = None           # OncoTree or None
        self.transcripts = {}          # unversioned accession -> [TranscriptModel]
        self.genome = None             # Genome or None
        self.chip_ab = None            # (chrom, pos, ref, alt) -> {dataset: (AN, VAF 0.3-0.65 alleles)}
        self.vcep = {}                 # gene -> VcepThresholds (hereditary cancer genes)
        self.gene_specific_splice = {}  # (gene, transcript) -> gene-specific splice provisions
        self.functional_assays = {}   # gene -> [assay parameters]
        self.functional_scores = {}   # (assay_id, protein) -> % of wild-type activity

    # --------------------------------------------------------------- loaders
    def load_canonical(self, path):
        with open(path) as fh:
            self.canonical = list(csv.DictReader(fh, delimiter="\t"))
        for c in self.canonical:
            c["_protein"] = normalise_protein(c.get("protein"))

    def load_gene_specific_splice(self, path):
        self.gene_specific_splice = load_gene_specific_splice(path)

    def load_functional_assays(self, params_path, scores_path):
        """Assay parameters and the per-variant measurements they apply to [S T1 O10]."""
        with open(params_path, newline="") as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
                self.functional_assays.setdefault(r["gene"], []).append(r)
        if not os.path.exists(scores_path):
            return
        with open(scores_path, newline="") as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
                self.functional_scores[(r["assay_id"], r["protein"])] = fnum(r["pct_wt"])

    def add_o5_reference(self, gene, transcript, protein, source):
        parts = protein_parts(normalise_protein(protein))
        if parts and parts[1] == parts[2] and re.fullmatch(r"[A-Z]", parts[3] or "") \
                and parts[0] != parts[3]:
            self.o5_references[gene].append((transcript or "", (parts[0], parts[1], parts[3]), source))

    def load_o5_references(self, path=None):
        """(Likely) oncogenic missense reference variants for O5: the SVIG-UK canonical
        list, plus an optional table (gene, transcript, protein, assessment, source)."""
        for c in self.canonical:
            if c.get("assessment") in ONC_CLASSES:
                self.add_o5_reference(c["gene"], c.get("transcript"), c.get("protein"),
                                      "SVIG-UK canonical list")
        if path:
            with _open_text(path) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r.get("assessment") in ONC_CLASSES:
                        self.add_o5_reference(r["gene"], r.get("transcript"), r["protein"], r.get("source", path))

    def load_gene_roles(self, cgc_path=None, oncokb_path=None):
        roles, sources = {}, []
        if oncokb_path and os.path.exists(oncokb_path):
            mapping = {"ONCOGENE": {"ONC"}, "TSG": {"TSG"}, "ONCOGENE_AND_TSG": {"ONC", "TSG"}}
            with open(oncokb_path) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r.get("geneType") in mapping:
                        roles[r["hugoSymbol"]] = set(mapping[r["geneType"]])
            sources.append("OncoKB (fallback)")
        if cgc_path and os.path.exists(cgc_path):
            with open(cgc_path) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    tokens = {t.strip().lower() for t in r.get("role_in_cancer", "").split(",")}
                    s = set()
                    if "oncogene" in tokens:
                        s.add("ONC")
                    if "tsg" in tokens:
                        s.add("TSG")
                    if s:
                        roles[r["gene"]] = s  # CGC takes precedence [S T2 B2]
            sources.insert(0, "COSMIC Cancer Gene Census")
        self.gene_roles = roles
        self.gene_roles_source = " + ".join(sources)

    def load_transcript_models(self, path, accessions=None):
        with _open_text(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 11 or not f[1].startswith("NM_") or "_" in f[2]:
                    continue  # coding RefSeq transcripts on primary assembly contigs only
                base = _strip_version(f[1])
                if accessions is not None and base not in accessions:
                    continue
                cds_start, cds_end = int(f[6]), int(f[7])
                if cds_end <= cds_start:
                    continue
                starts = [int(x) for x in f[9].rstrip(",").split(",")]
                ends = [int(x) for x in f[10].rstrip(",").split(",")]
                self.transcripts.setdefault(base, []).append(
                    TranscriptModel(f[1], f[2], f[3], starts, ends, cds_start, cds_end))

    def load_vcep_thresholds(self, path):
        def num(v):
            return float(v) if str(v).strip() != "" else None
        with open(path) as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
                self.vcep[r["gene"]] = VcepThresholds(
                    r["gene"], r["vcep_cspec_id"], r["vcep_version"], num(r["o3_max_af"]),
                    r["o3_inclusive"] == "yes", num(r["o3_max_af_multiallelic_group"]),
                    num(r["o3_max_af_single_allele"]), r["o3_absent_only"] == "yes",
                    float(r["b1_standalone_af"]), r["b1_standalone_inclusive"] == "yes",
                    float(r["b1_strong_af"]), r["b1_strong_inclusive"] == "yes", int(r["b1_min_alleles"]))

    def load_chip_allele_balance(self, path):
        self.chip_ab = {}
        with _open_text(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                hist = [int(x) for x in r["ab_hist_alt_bin_freq"].split("|") if x != ""]
                if len(hist) != 20:
                    continue
                key = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
                self.chip_ab.setdefault(key, {})[r["dataset"]] = (int(r["AN"]), sum(hist[CHIP_AB_BINS]))

    def transcript_for(self, row):
        """Transcript model for the row's annotated transcript. A versioned accession
        must match exactly; the model must contain the variant position."""
        chrom = str(row.get("Chromosome") or "")
        if chrom and not chrom.startswith("chr"):
            chrom = "chr" + chrom
        pos = first_int(row.get("Start_Position"))
        for tx in (row.get("Feature"), row.get("MANE_SELECT"), row.get("NM_Transcript")):
            tx = str(tx or "").strip()
            if not tx.startswith("NM_"):
                continue
            models = [m for m in self.transcripts.get(_strip_version(tx), [])
                      if (not chrom or m.chrom == chrom)
                      and (pos is None or m.tx_start - 10 <= pos <= m.tx_end + 10)
                      and ("." not in tx or m.name == tx)]
            if models:
                return max(models, key=lambda m: int(m.name.split(".")[1]) if "." in m.name else 0)
        return None

    def load_hotspots(self, path):
        with open(path) as fh:
            for e in json.load(fh):
                gene = e.get("hugoSymbol")
                if e.get("type") == "single residue":
                    pos = (e.get("aminoAcidPosition") or {}).get("start")
                    self.hotspots_single[(gene, pos)] = e
                else:
                    self.hotspots_indel[gene].append(e)

    def load_constraint(self, path):
        with _open_text(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for r in reader:
                z = fnum(r.get("mis.z_score"))
                if z is None:
                    continue
                gene = r.get("#gene") or r.get("gene")
                d = self.mis_z.setdefault(gene, {})
                d[_strip_version(r.get("transcript"))] = z
                if str(r.get("mane_select", "")).lower() == "true":
                    d["_mane"] = z

    # --------------------------------------------------------------- queries
    def roles(self, gene, group):
        override = CONTEXT_ROLE_OVERRIDES.get((gene, group))
        if override is not None:
            return set(override), "SVIG-UK tumour-type rule [G §7.2]"
        r = self.gene_roles.get(gene)
        return (set(r), self.gene_roles_source) if r else (set(), "")

    def mis_z_for(self, gene, transcripts):
        d = self.mis_z.get(gene)
        if not d:
            return None
        for tx in transcripts:
            key = _strip_version(tx)
            if key and key in d:
                return d[key]
        return d.get("_mane")


class CosmicO4:
    """Protein-level COSMIC sample counts built by utils/build_svig_uk_resources.py."""

    def __init__(self, path, genes=None):
        self.missense = defaultdict(dict)      # (gene, protein) -> {tx: n}
        self.inframe_ins = defaultdict(dict)
        self.truncating = defaultdict(lambda: defaultdict(list))  # gene -> tx -> [(start, n)]
        self.inframe_del = defaultdict(lambda: defaultdict(list))  # gene -> tx -> [(s, e, n)]
        with _open_text(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                gene = r["gene"]
                if genes is not None and gene not in genes:
                    continue
                n = int(r["n_samples"])
                tx, cls, prot = r["transcript"], r["mut_class"], r["protein"]
                s, e = first_int(r["aa_start"]), first_int(r["aa_end"])
                if cls == "missense":
                    self.missense[(gene, prot)][tx] = n
                elif cls == "inframe_ins":
                    self.inframe_ins[(gene, prot)][tx] = n
                elif cls == "truncating" and s is not None:
                    self.truncating[gene][tx].append((s, n))
                elif cls == "inframe_del" and s is not None:
                    self.inframe_del[gene][tx].append((s, e if e is not None else s, n))

    # Counts are taken from the single COSMIC transcript with the most entries so
    # that the same mutation reported on several transcripts is not double counted.
    def missense_count(self, gene, protein):
        return max(self.missense.get((gene, protein), {}).values(), default=0)

    def insertion_count(self, gene, protein):
        return max(self.inframe_ins.get((gene, protein), {}).values(), default=0)

    def truncating_count(self, gene, start):
        """Entries 'at the same position or further downstream' [S T1 O4]."""
        return max((sum(n for s, n in lst if s >= start)
                    for lst in self.truncating.get(gene, {}).values()), default=0)

    def deletion_count(self, gene, start, end):
        """Deletions that 'exactly match the VUA or are smaller than and within
        the deleted region' [S T1 O4]."""
        return max((sum(n for s, e, n in lst if s >= start and e <= end)
                    for lst in self.inframe_del.get(gene, {}).values()), default=0)


# =============================================================================
# Evidence
# =============================================================================


@dataclass
class Evidence:
    code: str
    strength: str          # VSTR/STR/MOD/SUPP, or STANDALONE
    points: int
    reason: str

    def label(self):
        if self.strength == "STANDALONE":
            return f"{self.code}_standalone [{self.reason}]"
        return f"{self.code}_{self.strength}({self.points:+d}) [{self.reason}]"


def onc(code, strength, reason):
    return Evidence(code, strength, POINTS[strength], reason)


def ben(code, strength, reason):
    return Evidence(code, strength, -POINTS[strength], reason)


@dataclass
class Result:
    classification: str
    score: object
    codes: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    tumour_type: str = ""
    tumour_group: str = ""
    gene_role: str = ""

    def codes_str(self):
        return "; ".join(e.label() for e in self.codes)

    def notes_str(self):
        return "; ".join(self.notes)


def class_from_points(points):
    for threshold, name in CLASS_THRESHOLDS:
        if points >= threshold:
            return name
    return "Benign"


# =============================================================================
# Rules
# =============================================================================


def rule_o1(row, res, protein):
    """O1 — SVIG-UK canonical variants list: stand-alone classification [S T3]."""
    gene = row.get("Hugo_Symbol", "")
    chrom = str(row.get("Chromosome", "")).replace("chr", "")
    for c in res.canonical:
        if c.get("gene") != gene:
            continue
        rule = c.get("match_rule", "protein")
        hit = False
        if rule == "protein":
            hit = bool(protein) and protein == c["_protein"]
        elif rule.startswith("genomic:"):
            _, cc, pos, ref, alt = rule.split(":")
            hit = (chrom == cc.replace("chr", "") and str(row.get("Start_Position")) == pos
                   and row.get("Reference_Allele") == ref and row.get("Tumor_Seq_Allele2") == alt)
        elif rule == "flt3_itd":
            exon, _ = exon_fraction(row.get("EXON"))
            hit = (terms(row.get("Consequence")) & {"inframe_insertion"}
                   and exon == 14) or ("dup" in protein and exon == 14)
        if hit:
            return c
    return None


@dataclass
class VcepThresholds:
    gene: str
    cspec_id: str
    version: str
    o3_max_af: object
    o3_inclusive: bool
    o3_max_af_multiallelic_group: object
    o3_max_af_single_allele: object
    o3_absent_only: bool
    b1_standalone_af: float
    b1_standalone_inclusive: bool
    b1_strong_af: float
    b1_strong_inclusive: bool
    b1_min_alleles: int

    @property
    def label(self):
        return f"ClinGen {self.gene} VCEP {self.cspec_id} v{self.version}"


def faf95(allele_count, allele_number):
    """gnomAD-style filtering allele frequency: the highest true AF for which the observed allele
    count is still within the 95% Poisson confidence bound. Normal approximation above 1,000 alleles."""
    import math
    ac, an = allele_count, allele_number
    if ac <= 0 or an <= 0:
        return 0.0
    if ac > 1000:
        return max(0.0, (ac - 1.645 * math.sqrt(ac)) / an)

    def upper_tail(lam):  # P(X >= ac) for X ~ Poisson(lam)
        term, cdf = math.exp(-lam), 0.0
        for k in range(int(ac)):
            cdf += term
            term *= lam / (k + 1)
        return 1 - cdf

    lo, hi = 0.0, ac / an
    for _ in range(60):
        mid = (lo + hi) / 2
        if upper_tail(mid * an) < 0.05:
            lo = mid
        else:
            hi = mid
    return lo


def continental_group_stats(row):
    """[(dataset, group, af, allele_number, estimated allele count, faf95)] for groups with data."""
    stats = []
    for ds in ("gnomADe", "gnomADg"):
        for grp_name in CONTINENTAL_GROUPS:
            af = fnum(row.get(f"{ds}_{grp_name}_AF"))
            n = GNOMAD_SAMPLES.get((ds, grp_name))
            if af is None or not n:
                continue
            an = 2 * n
            ac = round(af * an)
            stats.append((ds, grp_name, af, an, ac, faf95(ac, an)))
    return stats


def _above(value, threshold, inclusive):
    return value >= threshold if inclusive else value > threshold


def _below(value, threshold, inclusive):
    return value <= threshold if inclusive else value < threshold


def rule_b1_vcep(row, t):
    groups = [g for g in continental_group_stats(row) if g[3] >= VCEP_MIN_GROUP_ALLELES and g[4] >= t.b1_min_alleles]
    if not groups:
        return None
    ds, grp_name, af, an, ac, faf = max(groups, key=lambda g: g[5])
    detail = f"GrpMax FAF95 {faf:.3g} ({ds} {grp_name}, ~{ac} alleles)"
    if _above(faf, t.b1_standalone_af, t.b1_standalone_inclusive):
        return Evidence("B1", "STANDALONE", 0, f"{detail} meets BA1 of {t.label} ({t.b1_standalone_af:g})")
    if _above(faf, t.b1_strong_af, t.b1_strong_inclusive):
        return ben("B1", "STR", f"{detail} meets BS1 of {t.label} ({t.b1_strong_af:g})")
    return None


def rule_b1(row, res=None):
    """B1 — high population frequency: >1% stand-alone benign, >0.1% strong [S T2 B1].
    Overall AF or the highest of the five continental groups, exomes or genomes. Hereditary
    cancer genes with ClinGen VCEP rules use the panel's BA1/BS1 cut-offs (see VcepThresholds)."""
    t = res.vcep.get(row.get("Hugo_Symbol", "")) if res is not None else None
    if t is not None:
        return rule_b1_vcep(row, t)
    vals = [fnum(row.get("gnomADe_AF")), fnum(row.get("gnomADg_AF"))]
    for ds in ("gnomADe", "gnomADg"):
        vals += [fnum(row.get(f"{ds}_{pop}_AF")) for pop in B1_POPULATIONS]
    vals = [v for v in vals if v is not None]
    if not vals:
        return None
    top = max(vals)
    if top > 0.01:
        return Evidence("B1", "STANDALONE", 0, f"gnomAD v4.1 max AF {top:.4g} > 1%")
    if top > 0.001:
        return ben("B1", "STR", f"gnomAD v4.1 max AF {top:.4g} > 0.1%")
    return None


def o3_threshold(gene, res=None):
    """Allele frequency at or below which O3 applies at supporting strength."""
    t = res.vcep.get(gene) if res is not None else None
    if t is None:
        return 1e-5, "0.001%"
    if t.o3_absent_only:
        return 0.0, f"absent only ({t.label})"
    return t.o3_max_af, f"{t.o3_max_af:g} ({t.label})"


def rule_o3_vcep(row, t):
    if t.o3_absent_only:
        return None
    overall = max([v for v in (fnum(row.get("gnomADe_AF")), fnum(row.get("gnomADg_AF"))) if v is not None], default=0)
    total_alleles = 2 * max(sum(n for (ds, _), n in GNOMAD_SAMPLES.items() if ds == "gnomADe"),
                            sum(n for (ds, _), n in GNOMAD_SAMPLES.items() if ds == "gnomADg"))
    limit = t.o3_max_af
    if t.o3_max_af_single_allele is not None and round(overall * total_alleles) <= 1:
        limit = t.o3_max_af_single_allele
    if not _below(overall, limit, t.o3_inclusive):
        return None
    groups = continental_group_stats(row)
    for ds, grp_name, af, an, ac, faf in groups:
        if faf and not _below(faf, t.o3_max_af, t.o3_inclusive):
            return None
        if ac >= 2 and t.o3_max_af_multiallelic_group is not None \
                and not _below(af, t.o3_max_af_multiallelic_group, t.o3_inclusive or t.o3_max_af_multiallelic_group == 0):
            return None
    return onc("O3", "SUPP", f"gnomAD v4.1 AF {overall:.3g} within PM2_Supporting of {t.label}")


def rule_o3(row, res=None):
    """O3 — absent (+2) or very rare (+1) in gnomAD v4.1 [S T1 O3]. Hereditary cancer genes with
    ClinGen VCEP rules use the panel's PM2 cut-off for the very-rare branch."""
    afs = {k: fnum(row.get(f"{k[0]}_{k[1]}_AF")) for k in GNOMAD_SAMPLES}
    overall = [fnum(row.get("gnomADe_AF")), fnum(row.get("gnomADg_AF"))]
    if not any(v for v in list(afs.values()) + overall):
        return onc("O3", "MOD", "absent from gnomAD v4.1 exomes and genomes")
    t = res.vcep.get(row.get("Hugo_Symbol", "")) if res is not None else None
    if t is not None:
        return rule_o3_vcep(row, t)
    for k, af in afs.items():
        alleles = 2 * GNOMAD_SAMPLES[k]
        if not af:
            continue
        if alleles >= 100000 and af > 1e-5:
            return None
        if alleles < 100000 and af * alleles > 1.5:
            # "present in >1 allele in any sub-population consisting of <100,000 alleles"
            return None
    return onc("O3", "SUPP", "present in gnomAD v4.1 at <=0.001% and in <=1 allele in smaller sub-populations")


def rule_o3_chip(row, res):
    """O3 re-assessed after excluding likely CHIP alleles (see CHIP_GENES) [S T1 O3].
    Called only when neither O3 nor B1 applies. Returns (evidence or None, note or None).
    SMART proxy: the recalculated frequency is judged with the O3 thresholds on the
    combined exome + genome total, because gnomAD gives the allele-balance histogram
    only for all populations together."""
    gene = row.get("Hugo_Symbol", "")
    total = max([v for v in (fnum(row.get("gnomADe_AF")), fnum(row.get("gnomADg_AF"))) if v is not None],
                default=0)
    threshold, threshold_label = o3_threshold(gene, res)
    if gene not in CHIP_GENES or total <= threshold:
        return None, None
    if res.chip_ab is None:
        return None, "O3 CHIP correction not evaluated: gnomAD allele-balance table not provided"
    chrom = str(row.get("Chromosome") or "")
    chrom = chrom if chrom.startswith("chr") else "chr" + chrom
    pos = first_int(row.get("POS") or row.get("Start_Position"))
    entry = res.chip_ab.get((chrom, pos, row.get("REF") or row.get("Reference_Allele"),
                             row.get("ALT") or row.get("Tumor_Seq_Allele2")))
    if not entry:
        return None, "O3 CHIP correction not evaluated: variant not found in the gnomAD allele-balance table"
    an = sum(v[0] for v in entry.values())
    mid = sum(v[1] for v in entry.values())
    af = mid / an if an else 0
    detail = (f"{gene} CHIP correction: {mid} gnomAD v4.1 alleles with VAF 0.30-0.65 "
              f"(total AF {total:.3g}; recalculated AF {af:.3g})")
    if mid == 0:
        return onc("O3", "MOD", detail + ", absent after excluding low-VAF alleles"), None
    if af <= threshold:
        return onc("O3", "SUPP", f"{detail} <= {threshold_label}"), None
    return None, f"{detail} > {threshold_label}: O3 still not applicable"


def _protein_loss(tx, first_residue, context, residues=None):
    """Abou Tayoun 2018 branch for variants not predicted to undergo NMD: STR when
    the altered region is >10% of the protein, otherwise MOD. The altered region runs
    from first_residue to the C-terminus unless an in-frame length is given."""
    if residues is None:
        residues = tx.protein_length - first_residue + 1
    lost = residues / tx.protein_length
    detail = f"{context}; alters {lost:.0%} of the {tx.protein_length}-aa protein ({tx.name})"
    if lost > PVS1_PROTEIN_FRACTION:
        return onc("O2", "STR", detail + " > 10%")
    return onc("O2", "MOD", detail + " <= 10%")


def _vep_exon_total(row):
    _, exon_total = exon_fraction(row.get("EXON"))
    if exon_total is not None:
        return exon_total
    _, intron_total = exon_fraction(row.get("INTRON"))
    return intron_total + 1 if intron_total is not None else None


def _cryptic_site_nearby(row):
    for ds, dp in (("SpliceAI_pred_DS_AG", "SpliceAI_pred_DP_AG"),
                   ("SpliceAI_pred_DS_DG", "SpliceAI_pred_DP_DG")):
        score, pos = fnum(row.get(ds)), fnum(row.get(dp))
        if score is not None and score >= CRYPTIC_SITE_MIN_DS and pos is not None \
                and abs(pos) <= CRYPTIC_SITE_WINDOW_BP:
            return f"SpliceAI {ds[-2:]} {score:.2f} at {int(pos):+d} bp"
    return None


def translate_exonic_indel(hgvsc, tx, genome):
    """Apply a simple exonic c. deletion/duplication/insertion to the reference mRNA and
    return (first altered residue, new stop residue or None) for frameshifts, or None.
    Used when VEP reports no protein change (e.g. an indel flagged as a splice acceptor
    variant before 3' shifting, such as a 1-bp deletion of the first base of an exon)."""
    m = re.match(r"^c\.(\d+)(?:_(\d+))?(del|dup|ins)([ACGT]*)$", str(hgvsc or "").split(":")[-1])
    if not m or genome is None:
        return None
    start = int(m.group(1))
    end = int(m.group(2)) if m.group(2) else start
    mrna = "".join(tx.exon_sequence(genome, i) for i in range(1, tx.n_exons + 1))
    u = tx.utr5
    if m.group(3) == "del":
        mut = mrna[:u + start - 1] + mrna[u + end:]
    elif m.group(3) == "dup":
        mut = mrna[:u + end] + mrna[u + start - 1:u + end] + mrna[u + end:]
    else:
        if not m.group(2) or not m.group(4):
            return None
        mut = mrna[:u + start] + m.group(4) + mrna[u + start:]
    if (len(mut) - len(mrna)) % 3 == 0 or end > tx.cds_length:
        return None
    first = None
    for i in range(u + 3 * ((start - 1) // 3), len(mut) - 2, 3):
        codon = mut[i:i + 3]
        if first is None and GENETIC_CODE.get(codon) != GENETIC_CODE.get(mrna[i:i + 3]):
            first = (i - u) // 3 + 1  # HGVS: first residue whose amino acid changes
        if codon in STOP_CODONS:
            return first, (i - u) // 3 + 1
    return first, None


def _o2_truncating(row, tx, genome=None):
    parts = ptc_residues(row.get("HGVSp") or row.get("HGVSp_Short"))
    if not parts:
        parts = translate_exonic_indel(row.get("HGVSc"), tx, genome)
    if not parts or not 1 <= parts[0] <= tx.protein_length:
        return None
    first, stop = parts
    if stop is None:
        return _protein_loss(tx, first, "no new in-frame stop codon reported (fs*?); NMD not predicted")
    # stop residue is numbered on the variant protein; shift back to reference mRNA coordinates
    mrna_pos = tx.utr5 + 3 * (stop - 1) + 1 - indel_length_change(row)
    nmd, where = nmd_predicted(tx.exon_lengths, mrna_pos)
    if nmd:
        return onc("O2", "VSTR", f"PTC at residue {stop} {where}, >{NMD_PENULTIMATE_EXON_BP} bp before "
                                 f"the last exon junction: NMD predicted ({tx.name})")
    return _protein_loss(tx, first, f"PTC at residue {stop} {where}: NMD escape")


def _o2_skipped_exon(row, tx, genome):
    """Canonical +/-1,2 variant: the exon nearest the affected intron is assumed skipped."""
    intron_n, _ = exon_fraction(row.get("INTRON"))
    offset = intronic_offset(row.get("HGVSc"))
    if intron_n is None or offset is None:
        return None
    # "it is assumed that the exon it is closest to is skipped" [S T1 O2]
    k = intron_n if offset > 0 else intron_n + 1
    return _o2_exon_skip_effect(row, tx, genome, k)


def _o2_exon_skip_effect(row, tx, genome, k):
    """Strength for skipping exon k, following the PVS1 splice branch."""
    if not 1 <= k <= tx.n_exons:
        return None
    coding = tx.coding[k - 1]
    label = f"skipping of exon {k}/{tx.n_exons}"
    if coding == 0:
        return "not_applicable"
    if tx.exon_has_start(k):
        return onc("O2", "SUPP", f"{label} removes the initiation codon; start-loss branch, "
                                 f"alternative start sites not assessable ({tx.name})")
    if k == tx.n_exons or k == 1:
        return onc("O2", "MOD", f"{label}: first/final intron, impact difficult to predict, "
                                f"lowest strength applied [S T1 O2; Walker 2023] ({tx.name})")
    cryptic = _cryptic_site_nearby(row)
    if cryptic:
        return onc("O2", "MOD", f"possible cryptic splice site within {CRYPTIC_SITE_WINDOW_BP} bp ({cryptic}): "
                                f"in-frame and out-of-frame scenarios possible, lowest strength applied")
    first_residue = sum(tx.coding[:k - 1]) // 3 + 1
    if tx.exon_has_stop(k):
        return _protein_loss(tx, first_residue, f"{label} removes the stop codon; NMD not predicted")
    if coding % 3 == 0:
        return _protein_loss(tx, first_residue, f"{label} ({coding} coding bp) preserves the reading frame",
                             residues=coding // 3)
    if genome is None:
        return _protein_loss(tx, first_residue, f"{label} ({coding} coding bp) disrupts the reading frame; "
                                                f"NMD not assessable without the reference genome")
    exons = [tx.exon_sequence(genome, i) for i in range(1, tx.n_exons + 1) if i != k]
    lengths = [len(e) for e in exons]
    mrna = "".join(exons)
    stop_pos = None
    for i in range(tx.utr5, len(mrna) - 2, 3):
        if mrna[i:i + 3] in STOP_CODONS:
            stop_pos = i + 1
            break
    if stop_pos is None:
        return _protein_loss(tx, first_residue, f"{label} disrupts the reading frame with no downstream stop codon")
    nmd, where = nmd_predicted(lengths, stop_pos)
    new_stop = (stop_pos - tx.utr5 - 1) // 3 + 1
    if nmd:
        return onc("O2", "VSTR", f"{label} ({coding} coding bp) disrupts the reading frame; new stop at residue "
                                 f"{new_stop} {where} of the skipped transcript: NMD predicted ({tx.name})")
    return _protein_loss(tx, first_residue, f"{label} disrupts the reading frame; new stop {where}: NMD escape")


def exon_cds_span(tx, k):
    """1-based CDS range (first, last) of exon k, or None when the exon is untranslated."""
    if not 1 <= k <= tx.n_exons or tx.coding[k - 1] == 0:
        return None
    first = sum(tx.coding[:k - 1]) + 1
    return first, first + tx.coding[k - 1] - 1


def exon_edge_position(tx, cds_pos):
    """(exon, label) when a coding position is the first base or one of the last three
    bases of its exon, where synonymous changes are most likely to disrupt the donor or
    acceptor motif [S T2 B4, p.16]. Otherwise (exon, None)."""
    for k in range(1, tx.n_exons + 1):
        span = exon_cds_span(tx, k)
        if span is None or not span[0] <= cds_pos <= span[1]:
            continue
        if cds_pos == span[0]:
            return k, "first base of the exon"
        if span[1] - cds_pos <= 2:
            nth = {0: "last", 1: "second to last", 2: "third to last"}[span[1] - cds_pos]
            return k, f"{nth} base of the exon"
        return k, None
    return None, None


def load_gene_specific_splice(path):
    """Gene-specific splice provisions: (gene, transcript without version) -> entries."""
    table = {}
    if not path or not os.path.exists(path):
        return table
    with open(path, newline="") as fh:
        for row in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            key = (row["gene"], row["transcript"].split(".")[0])
            table.setdefault(key, []).append({
                "exon": int(row["exon"]), "cds_position": int(row["cds_position"]),
                "applies_to": row.get("applies_to", "synonymous"),
                "assessment": row.get("assessment", ""), "source": row.get("source", "")})
    return table


def gene_specific_splice(row, res, vclass):
    """The entry matching this variant, when a gene-specific splice provision applies."""
    table = getattr(res, "gene_specific_splice", None)
    if not table:
        return None
    entries = table.get((row.get("Hugo_Symbol", ""), _strip_version(row.get("NM_Transcript") or row.get("Feature"))))
    if not entries:
        return None
    cds = first_int(row.get("CDS_position"))
    if cds is None:
        m = re.match(r"^c\.(\d+)[ACGT]>[ACGT]$", str(row.get("HGVSc") or ""))
        cds = int(m.group(1)) if m else None
    if cds is None:
        return None
    for e in entries:
        if e["cds_position"] == cds and vclass in e["applies_to"].split("|"):
            return e
    return None


def rule_o2(row, vclass, res=None):
    """O2 — null variant in a TSG, PVS1-based strength [S T1 O2; S F3; Abou Tayoun 2018].
    With a transcript model (RefSeq exon structure) the tree is followed: NMD is
    predicted from the PTC position, and variants escaping NMD are weighted by the
    fraction of protein altered (>10% STR, otherwise MOD). Canonical splice variants
    are assumed to skip the nearest exon; the new reading frame is translated from the
    reference genome when available. SMART proxy: whether a region is 'critical to
    protein function' cannot be judged, so only the >10% branch is used. Without a
    transcript model (or when its exon count disagrees with VEP) the lowest strength
    reachable in the tree is used."""
    exon_n, exon_total = exon_fraction(row.get("EXON"))
    tx = res.transcript_for(row) if res is not None and res.transcripts else None
    if tx is not None and _vep_exon_total(row) != tx.n_exons:
        tx = None
    if vclass == "truncating":
        cds = first_int(row.get("CDS_position"))
        if exon_n == 1 and cds is not None and cds <= 100:
            return onc("O2", "MOD", "stop within the first 100 bp of exon 1; re-initiation not assessable")
        if tx is not None:
            e = _o2_truncating(row, tx, res.genome)
            if e:
                return e
        if exon_n is None or exon_total is None:
            return onc("O2", "MOD", "exon position unknown; NMD not assessable")
        if exon_n == exon_total:
            return onc("O2", "MOD", "premature stop in the final exon (NMD escape; protein fraction not assessable)")
        if exon_n == exon_total - 1:
            return onc("O2", "MOD", "premature stop in the penultimate exon; NMD boundary not assessable")
        return onc("O2", "VSTR", "premature stop predicted to trigger NMD")
    if vclass == "canonical_splice":
        hgvsc = str(row.get("HGVSc") or "")
        sai = spliceai_scores(row)
        if re.search(r"\+2T>C", hgvsc) and (sai is None or max(sai) < 0.8):
            return None  # "+2T>C ... PVS1 should not be applied ... SpliceAI delta score <0.8"
        if re.search(r"\+2C>T", hgvsc):
            return None  # "+2C>T ... O2 is not applicable"
        if touches_canonical_splice_site(hgvsc) is False:
            # not a +/-1,2 variant once described on the 3'-shifted allele; an exonic
            # frameshift is assessed as a truncating variant instead
            e = _o2_truncating(row, tx, res.genome) if tx is not None else None
            return e
        if tx is not None:
            e = _o2_skipped_exon(row, tx, res.genome)
            if e == "not_applicable":
                return None  # skipped exon is untranslated: no effect on the protein
            if e:
                return e
        return onc("O2", "MOD", "canonical +/-1,2 splice site; in-frame vs out-of-frame not assessable, lowest strength applied")
    if vclass == "donor_last_base":
        # A gene-specific provision says this exonic position disrupts the splice motif, so the
        # exon it sits in is the one assumed skipped [G 7.11; S T2 B4 p.16].
        if tx is not None and exon_n is not None:
            e = _o2_exon_skip_effect(row, tx, res.genome, exon_n)
            if e == "not_applicable":
                return None
            if e:
                return e
        return onc("O2", "MOD", "exonic position disrupting the splice motif (gene-specific provision); "
                                "exon skipping assumed, in-frame vs out-of-frame not assessable, "
                                "lowest strength applied")
    if vclass == "start_lost":
        return onc("O2", "SUPP", "initiation codon loss; alternative start / upstream oncogenic variants not assessable")
    return None


def rule_o9_b5(row, vclass, protein, roles):
    """O9 — in-frame length change in a functional domain, or truncation in the final
    exon of an oncogene predicted to be gain-of-function [S T1 O9].
    B5 — in-frame deletion/insertion in a repetitive region of unknown function [S T2 B5].
    SMART proxy: VEP DOMAINS defines 'known functional domain' and low-complexity
    sequence stands in for 'repetitive region'."""
    sources = domain_sources(row)
    functional = bool(sources & FUNCTIONAL_DOMAIN_DBS)
    if vclass == "inframe":
        if sources & REPEAT_DBS and not functional:
            return ben("B5", "SUPP", "in-frame indel in low-complexity sequence outside any functional domain")
        if functional:
            parts = protein_parts(protein)
            size = (parts[2] - parts[1] + 1) if parts else 1
            ins = re.search(r"ins([A-Z*]+)$", protein or "")
            if ins:
                size = max(size, len(ins.group(1)))
            if size > 1:
                return onc("O9", "MOD", f"in-frame change of {size} residues within a functional domain")
            return onc("O9", "SUPP", "single-residue in-frame change within a functional domain")
    if vclass == "truncating" and "ONC" in roles and "TSG" not in roles:
        exon_n, exon_total = exon_fraction(row.get("EXON"))
        if exon_n is not None and exon_n == exon_total:
            return onc("O9", "MOD", "truncating variant in the final exon of an oncogene (predicted gain-of-function)")
    return None


def rule_o5(row, res, vclass, protein):
    """O5 — different (likely) oncogenic amino acid change at the same residue [S T1 O5],
    missense route. Reference residues given on another transcript are mapped through
    the genome to the transcript SMART annotated; references without a transcript must
    agree on the reference amino acid. The reference variant's REVEL is the highest score
    among nucleotide changes producing that amino acid change. SMART proxy: the reference
    lists are curated catalogues rather than accredited-laboratory classifications, and
    the mechanism of action is assumed to be shared within a gene. The splice and
    initiation-codon routes are not automated."""
    if vclass != "missense" or res.revel is None:
        return None
    gene = row.get("Hugo_Symbol", "")
    parts = protein_parts(protein)
    vua_revel = fnum(row.get("REVEL"))
    tx = res.transcript_for(row) if res.transcripts else None
    if not parts or vua_revel is None or tx is None or not res.o5_references.get(gene):
        return None
    ref_aa, pos, _, alt_aa = parts
    codon = [tx.cds_to_genomic(3 * pos - k) for k in (2, 1, 0)]
    if None in codon:
        return None
    best = None
    for transcript, (r_ref, r_pos, r_alt), source in res.o5_references[gene]:
        if r_ref != ref_aa or r_alt == alt_aa:
            continue
        acc = _strip_version(transcript)
        if acc and acc != _strip_version(tx.name):
            models = res.transcripts.get(acc, [])
            if not models:
                continue
            g = models[0].cds_to_genomic(3 * r_pos - 2)
            if g is None or tx.residue_at(g) != pos:
                continue
        elif r_pos != pos:
            continue
        scores = [sc for p in codon for (_, _, aaref, aaalt, sc) in res.revel.at(tx.chrom, p)
                  if aaref == ref_aa and aaalt == r_alt]
        if not scores:
            continue
        ref_revel = max(scores)
        detail = (f"{source}: {gene} {r_ref}{pos}{r_alt} (REVEL {ref_revel:.3f}) at the same residue; "
                  f"variant REVEL {vua_revel:.3f}")
        if vua_revel > ref_revel - O5_REVEL_EQUIVALENT or \
                (vua_revel > O5_REVEL_BOTH_HIGH and ref_revel > O5_REVEL_BOTH_HIGH):
            e = onc("O5", "MOD", detail)
        elif O5_REVEL_SUPPORTING_MIN < vua_revel < O5_REVEL_BOTH_HIGH:
            e = onc("O5", "SUPP", detail)
        else:
            continue
        if best is None or e.points > best.points:
            best = e
    return best


def rule_o6_b3(row, vclass, roles, functional_region):
    """O6 / B3 — computational evidence, supporting only [G §7.6; S T1 O6; S T2 B3]."""
    revel = fnum(row.get("REVEL"))
    sai = spliceai_scores(row)
    sai_max = max(sai) if sai else None
    oncogene_only = "ONC" in roles and "TSG" not in roles
    o6 = None
    if vclass == "missense" and revel is not None and revel >= 0.7:
        # "For GOF genes, O6 ... should be limited to variants within functional regions"
        if not oncogene_only or functional_region:
            o6 = onc("O6", "SUPP", f"REVEL {revel:.3f} >= 0.7")
    if o6 is None and vclass not in SPLICE_ASSUMED_CLASSES and sai_max is not None and sai_max >= 0.2:
        o6 = onc("O6", "SUPP", f"SpliceAI max delta {sai_max:.2f} >= 0.2")
    if o6:
        return o6
    if vclass == "missense":
        # REVEL "should NOT be used as evidence of benignity for non-LOF (e.g. activating) variants"
        if revel is not None and revel <= 0.4 and "TSG" in roles and (sai_max is None or sai_max <= 0.1):
            return ben("B3", "SUPP", f"REVEL {revel:.3f} <= 0.4 in a loss-of-function gene")
    elif vclass in {"synonymous", "intronic", "utr"} and sai and max(sai) <= 0.1:
        return ben("B3", "SUPP", f"all SpliceAI deltas <= 0.1 (max {max(sai):.2f})")
    return None


def rule_b4(row, vclass, b3, synonymous_ok=True):
    """B4 — synonymous / deep intronic / UTR with no predicted splicing impact,
    supporting strength without an RNA assay [S T2 B4]. B4 "should be used with caution for
    synonymous variants ... located at the first base or last three bases of the exon"
    [S T2 B4, p.16]; synonymous_ok is False for those, so it is not applied."""
    if not b3 or not b3.reason.startswith("all SpliceAI"):
        return None
    if vclass == "synonymous" and synonymous_ok and not (terms(row.get("Consequence")) & SPLICE_REGION_TERMS):
        return ben("B4", "SUPP", "synonymous, outside the splice region, no predicted splicing impact")
    if vclass == "intronic":
        off = intronic_offset(row.get("HGVSc"))
        if off is not None and (off >= 7 or off <= -21):
            return ben("B4", "SUPP", f"intronic at {off:+d} (at or beyond +7/-21), no predicted splicing impact")
    if vclass == "utr":
        return ben("B4", "SUPP", "untranslated region, no predicted splicing impact")
    return None


def rule_o7(row, res, vclass, protein):
    """O7 — Cancer Hotspots route [S T1 O7]: >=50 entries at the residue with >=10
    for the same change (strong); <50 with >=10 (moderate); 2-9 same change (supporting).
    The functional-domain routes need judgement on local benign variation and are not used."""
    gene = row.get("Hugo_Symbol", "")
    parts = protein_parts(protein)
    if not parts:
        return None
    total = same = 0
    if vclass == "missense":
        e = res.hotspots_single.get((gene, parts[1]))
        if not e or not str(e.get("residue", "")).startswith(parts[0]):
            return None
        alt = parts[3][:1]
        total, same = e.get("tumorCount") or 0, (e.get("variantAminoAcid") or {}).get(alt, 0)
    elif vclass == "inframe":
        for e in res.hotspots_indel.get(gene, []):
            if protein in (e.get("variantAminoAcid") or {}):
                total, same = e.get("tumorCount") or 0, e["variantAminoAcid"][protein]
                break
    else:
        return None
    detail = f"Cancer Hotspots {gene} {protein}: {same} same change / {total} at residue"
    if total >= 50 and same >= 10:
        return onc("O7", "STR", detail)
    if same >= 10:
        return onc("O7", "MOD", detail)
    if 2 <= same <= 9:
        return onc("O7", "SUPP", detail)
    return None


def rule_o8(row, res, vclass):
    """O8 — missense constraint, gnomAD missense Z >= 3.09 [S T1 O8]."""
    if vclass != "missense":
        return None
    z = res.mis_z_for(row.get("Hugo_Symbol", ""),
                      [row.get("MANE_SELECT"), row.get("NM_Transcript"), row.get("Feature")])
    if z is not None and z >= 3.09:
        return onc("O8", "SUPP", f"gnomAD missense Z {z:.2f} >= 3.09")
    return None


def rule_o4(row, res, vclass, protein):
    """O4 — enrichment in a somatic database [S T1 O4] (SMART proxy: COSMIC, see O4_THRESHOLDS)."""
    if res.cosmic is None:
        return None
    gene = row.get("Hugo_Symbol", "")
    parts = protein_parts(protein)
    if not parts:
        return None
    if vclass == "missense":
        cls, n = "missense", res.cosmic.missense_count(gene, protein)
    elif vclass == "truncating":
        cls, n = "truncating", res.cosmic.truncating_count(gene, parts[1])
    elif vclass == "inframe":
        cls = "inframe"
        if "del" in protein and "ins" not in protein:
            n = res.cosmic.deletion_count(gene, parts[1], parts[2])
        else:
            n = res.cosmic.insertion_count(gene, protein)
    else:
        return None
    strong, moderate, supporting = O4_THRESHOLDS[cls]
    detail = f"COSMIC samples = {n} ({cls})"
    if n > strong:
        return onc("O4", "STR", detail)
    if n >= moderate:
        return onc("O4", "MOD", detail)
    if supporting is not None and n >= supporting:
        return onc("O4", "SUPP", detail)
    return None


def genie_matches(row, res, vclass, protein, entries):
    """GENIE entries that count towards O4 for this variant [S T1 O4]: same amino acid change
    (missense), exact nucleotide change (canonical splice), frameshift/nonsense at the same residue or
    downstream (truncating), exact insertion or deletions equal to or within the deleted region
    (in-frame). Entries are matched by genomic change and residues are computed on the transcript SMART
    annotated. Returns (O4 class, description, matched entries, note); class is None when not countable."""
    gene = row.get("Hugo_Symbol", "")
    allele = maf_allele(row, res.genome)
    tx = res.transcript_for(row) if res.transcripts else None
    parts = protein_parts(protein)

    def residue(e):
        return entry_residue(tx, e, res.genome)

    if vclass == "missense":
        if tx is None or not parts:
            return None, "", [], "O4 (GENIE) not evaluated: transcript model or protein change unavailable"
        ref_aa, pos, _, alt_aa = parts
        return "missense", f"{gene} {protein}", [
            e for e in entries if e.vclass == "missense" and e.ref_aa == ref_aa and e.alt_aa == alt_aa
            and residue(e) == pos], None
    if vclass == "canonical_splice":
        if allele is None:
            return None, "", [], None
        return "splice", "exact nucleotide change", [
            e for e in entries if (e.chrom, e.start, e.end, e.ref, e.alt) == allele], None
    if vclass == "truncating":
        if tx is None:
            return None, "", [], "O4 (GENIE) not evaluated: transcript model unavailable"
        ptc = ptc_residues(row.get("HGVSp") or row.get("HGVSp_Short")) or \
            translate_exonic_indel(row.get("HGVSc"), tx, res.genome)
        if not ptc:
            return None, "", [], None
        first = ptc[0]
        return "truncating", f"frameshift/nonsense at residue >= {first}", [
            e for e in entries if e.vclass == "truncating"
            and (r := residue(e)) is not None and first <= r <= tx.protein_length], None
    if vclass == "inframe" and allele is not None:
        chrom, start, end, ref, alt = allele
        if alt == "-":
            return "inframe", "in-frame deletions equal to or within the deleted region", [
                e for e in entries if e.vclass == "inframe_del" and e.chrom == chrom
                and start <= e.start and e.end <= end], None
        return "inframe", "exact in-frame insertion", [
            e for e in entries if e.vclass in ("inframe_ins", "inframe_del")
            and (e.chrom, e.start, e.end, e.ref, e.alt) == allele], None
    return None, "", [], None


def rule_o4_genie(row, res, vclass, protein, tumour_type):
    """O4 — case counting in GENIE [S T1 O4]; see genie_matches(). On-target entries follow
    --o4-on-target (see O4_ON_TARGET_LEVELS). Returns (evidence or None, note or None)."""
    gene = row.get("Hugo_Symbol", "")
    on_target, target_label = genie_target(tumour_type, res, gene)
    entries = [e for e in res.genie.by_gene.get(gene, []) if on_target(e)]
    cls, what, matched, note = genie_matches(row, res, vclass, protein, entries)
    if cls is None:
        return None, note
    n = sum(e.n for e in matched)
    strong, moderate, supporting = GENIE_O4_THRESHOLDS[cls]
    detail = f"GENIE {n} patients ({what}; {target_label})"
    if n > strong:
        return onc("O4", "STR", detail), None
    if n >= moderate:
        return onc("O4", "MOD", detail), None
    if supporting is not None and n >= supporting:
        return onc("O4", "SUPP", detail), None
    return None, None


DOWNGRADE = {"VSTR": "STR", "STR": "MOD", "MOD": "SUPP", "SUPP": None}


def rule_o10_functional(row, res, roles, vclass, protein):
    """O10 / B6 from a functional assay [S T1 O10 p.13; S T2 B6 p.17].

    "Variant is considered to have functionally abnormal effect if protein activity assay or
    functional impact is <25% of wildtype level"; an intermediate result may be downgraded "by
    one or more levels". The boundary for abnormal, the boundary for intermediate and the
    strength are all read from the assay table, because the guideline fixes only the first and
    leaves the weighting to ClinGen SVI / CanVIG-UK calibration of the individual assay.
    Applied at variant level only, and only where the direction agrees with the gene's role.
    """
    if not protein:
        return None
    for assay in res.functional_assays.get(row.get("Hugo_Symbol", ""), []):
        if vclass not in (assay.get("applies_to") or vclass).split("|"):
            continue
        if assay.get("transcript") and _strip_version(assay["transcript"]) != \
                _strip_version(row.get("NM_Transcript") or row.get("Feature")):
            continue
        pct = res.functional_scores.get((assay["assay_id"], protein))
        if pct is None:
            continue
        abnormal = fnum(assay.get("abnormal_max_pct"))
        intermediate = fnum(assay.get("intermediate_max_pct"))
        normal = fnum(assay.get("normal_min_pct"))
        strength = assay.get("strength") or "SUPP"
        measure = f"{pct:.1f}% of wild-type activity ({assay['assay_id']}, PMID {assay['pmid']})"
        if abnormal is not None and pct < abnormal:
            # a loss-of-function result supports oncogenicity only in a tumour suppressor
            if "TSG" in roles:
                return onc("O10", strength, f"{measure}: below {abnormal:g}%, functionally abnormal")
            continue
        if intermediate is not None and pct < intermediate:
            lower = DOWNGRADE.get(strength)
            if lower and "TSG" in roles:
                return onc("O10", lower, f"{measure}: intermediate (below {intermediate:g}%), "
                                         f"downgraded from {strength}")
            continue
        if normal is not None and pct >= normal:
            return ben("B6", strength, f"{measure}: at or above {normal:g}%, no damaging effect")
    return None


def rule_o10_b6(row, roles, vclass):
    """O10 / B6 — SMART proxy from OncoKB variant-level mutation effect, used where no functional
    assay covers the variant (see ONCOKB_* notes)."""
    if str(row.get("VARIANT_IN_ONCOKB", "")).lower() != "true":
        return None
    effect = str(row.get("MUTATION_EFFECT") or "").strip()
    if effect in ONCOKB_GOF and "ONC" in roles:
        return onc("O10", "SUPP", f"OncoKB mutation effect '{effect}' in an oncogene (proxy)")
    if effect in ONCOKB_LOF and "TSG" in roles:
        return onc("O10", "SUPP", f"OncoKB mutation effect '{effect}' in a tumour suppressor (proxy)")
    if effect in ONCOKB_NEUTRAL:
        return ben("B6", "SUPP", f"OncoKB mutation effect '{effect}' (proxy)")
    return None


# =============================================================================
# Classification
# =============================================================================


def classify_row(row, res, tumour_type=None):
    """Classify one MAF row; the Result also records the tumour type and gene role used."""
    result = _classify_row(row, res, tumour_type)
    tumour_value = row.get("TUMOUR_TYPE") or tumour_type or "unknown"
    group = tumour_group(tumour_value)
    code, tissue, _ = resolve_tumour_type(tumour_value, res.oncotree)
    if code:
        group = tissue.lower() if tissue in HAEM_TISSUES else "solid"
    roles, source = res.roles(row.get("Hugo_Symbol", ""), group)
    names = {"ONC": "oncogene", "TSG": "tumour suppressor"}
    result.tumour_type = str(tumour_value)
    result.tumour_group = group
    result.gene_role = (" + ".join(names[r] for r in sorted(roles)) + f" ({source})") if roles else "unknown"
    return result


def _classify_row(row, res, tumour_type=None):
    tumour_value = row.get("TUMOUR_TYPE") or tumour_type
    group = tumour_group(tumour_value)
    code, tissue, _ = resolve_tumour_type(tumour_value, res.oncotree)
    if code:  # OncoTree code: group from its tissue (e.g. DLBCLNOS -> lymphoid)
        group = tissue.lower() if tissue in HAEM_TISSUES else "solid"
    gene = row.get("Hugo_Symbol", "")
    protein = normalise_protein(row.get("HGVSp") or row.get("HGVSp_Short"))
    vclass = variant_class(row.get("Consequence"))
    roles, roles_src = res.roles(gene, group)
    notes = []

    # Synonymous variants at an exon edge: B4 "should be used with caution" there, and a few
    # positions carry a gene-specific provision that they disrupt splicing [G 7.11; S T2 B4 p.16].
    synonymous_b4_ok = True
    if vclass == "synonymous":
        provision = gene_specific_splice(row, res, vclass)
        tx_edge = res.transcript_for(row) if res.transcripts else None
        cds_pos = first_int(row.get("CDS_position"))
        edge = exon_edge_position(tx_edge, cds_pos)[1] if (tx_edge and cds_pos) else None
        if provision:
            vclass = "donor_last_base"
            notes.append(f"gene-specific provision: {provision['source']}; the guideline classes these as "
                         f"{provision['assessment']} variants, which is not a code or a strength - "
                         f"scientist review required")
        elif edge:
            synonymous_b4_ok = False
            notes.append(f"synonymous at the {edge}: B4 not applied, higher likelihood of splice-motif "
                         f"disruption [S T2 B4 p.16]")

    # ---- stand-alone codes -------------------------------------------------
    canon = rule_o1(row, res, protein)
    if canon:
        note = f"SVIG-UK canonical list: {canon['gene']} {canon.get('protein') or canon.get('hgvsc')}"
        return Result(canon["assessment"], "", [Evidence("O1", "STANDALONE", 0, note)],
                      ["stand-alone O1 [G §6.3]"])

    b1 = rule_b1(row, res)
    if b1 and b1.strength == "STANDALONE":
        return Result("Benign", "", [b1], ["stand-alone B1 [G §6.3]"])

    null_in_oncogene = (vclass in NULL_CLASSES and "ONC" in roles and "TSG" not in roles)
    if null_in_oncogene:
        exon_n, exon_total = exon_fraction(row.get("EXON"))
        final_exon_truncation = vclass == "truncating" and exon_n is not None and exon_n == exon_total
        if not final_exon_truncation:
            b2 = Evidence("B2", "STANDALONE", 0,
                          f"null variant in {gene}, an oncogene without a tumour-suppressor role ({roles_src})")
            return Result("VUS", "", [b2], ["stand-alone B2 [S T2 B2]"])

    if not roles:
        notes.append(f"mode of action unknown for {gene}: O2/B2 not evaluable")

    # ---- collect evidence --------------------------------------------------
    codes = {}

    def add(e):
        if e:
            codes[e.code] = e

    add(b1)
    add(rule_o3(row, res))
    if "O3" not in codes and b1 is None:
        o3_chip, chip_note = rule_o3_chip(row, res)
        add(o3_chip)
        if chip_note:
            notes.append(chip_note)
    if vclass in NULL_CLASSES and "TSG" in roles:
        add(rule_o2(row, vclass, res))
        if not res.transcripts:
            notes.append("O2 without transcript structure: transcript models not provided, lowest strength used")
    add(rule_o9_b5(row, vclass, protein, roles))
    functional_region = bool(domain_sources(row) & FUNCTIONAL_DOMAIN_DBS)
    o7 = rule_o7(row, res, vclass, protein) if res.hotspots_single else None
    if res.hotspots_single:
        add(o7)
    else:
        notes.append("O7 not evaluated: Cancer Hotspots counts not provided")
    b3_or_o6 = rule_o6_b3(row, vclass, roles, functional_region or o7 is not None)
    add(b3_or_o6)
    add(rule_b4(row, vclass, b3_or_o6 if b3_or_o6 and b3_or_o6.code == "B3" else None,
                synonymous_ok=synonymous_b4_ok))
    if res.mis_z:
        add(rule_o8(row, res, vclass))
        # A gene absent from the table is not a gene without missense constraint, and the
        # guard above tests only whether the *table* loaded. It always has, so this note
        # never appeared: loeuf_dataset_grch38.tsv.gz covers 22 autosomes and carries no
        # chrX or chrY rows at all, so mis_z_for returns None for every X-linked tumour
        # suppressor -- ATRX, KDM6A, BCOR, STAG2, DDX3X, PHF6, AMER1, MED12, KDM5C, ZRSR2,
        # GATA1, BTK, AR -- and O8 was recorded as "not met" when it had never been looked
        # up. gnomAD v4 does include chrX; this is a defect of our build of the table.
        #
        # A note only: the score and the codes are untouched, so no classification changes.
        if vclass == "missense" and res.mis_z_for(
                row.get("Hugo_Symbol", ""),
                [row.get("MANE_SELECT"), row.get("NM_Transcript"),
                 row.get("Feature")]) is None:
            notes.append(
                f"O8 not evaluated: {row.get('Hugo_Symbol', 'this gene')} is absent from the "
                f"gnomAD constraint table, so missense constraint was not looked up rather "
                f"than found to be below the threshold")
    else:
        notes.append("O8 not evaluated: gnomAD constraint not provided")
    if res.genie is not None:
        o4, o4_note = rule_o4_genie(row, res, vclass, protein, tumour_value)
        add(o4)
        if o4_note:
            notes.append(o4_note)
    elif res.cosmic is not None:
        add(rule_o4(row, res, vclass, protein))
    else:
        notes.append("O4 not evaluated: GENIE or COSMIC counts not provided")
    functional = rule_o10_functional(row, res, roles, vclass, protein)
    add(functional)
    if functional is None:
        add(rule_o10_b6(row, roles, vclass))
    elif res.functional_assays:
        notes.append(f"O10/B6 from a functional assay rather than the OncoKB proxy; the strength "
                     f"reflects the calibration recorded in the assay table")
    add(rule_o5(row, res, vclass, protein))

    # ---- combination rules [S F1, S F2; S T1/T2 guidance] -------------------
    def drop(code, why):
        if code in codes:
            notes.append(f"{codes.pop(code).code} removed: {why}")

    if "B1" in codes:
        drop("O4", "O4 cannot be used when B1 applies [S T1 O4]")
        drop("O7", "O7 requires consistency with mode of action, i.e. B1 not applied [S T1 O7]")
    if "O4" in codes and "O3" not in codes:
        drop("O4", "O3 is a prerequisite for O4 [S T1 O4]")
    if "O2" in codes:
        drop("O6", "not permitted with O2 [S F1 #2]")
        drop("O7", "not permitted with O2 [S F1 #4]")
        drop("O9", "not permitted with O2 [S F1 #10]")
        drop("O10", "O10 with O2 only for missense variants [S F1 #15]")
    if "O9" in codes:
        drop("O6", "not permitted with O9 [S F1 #12]")
    if "O5" in codes and "O7" in codes:
        # "O5 and O7 cannot be applied together; use evidence code which achieves the highest strength"
        weaker = "O5" if codes["O5"].points <= codes["O7"].points else "O7"
        drop(weaker, "O5 and O7 cannot be applied together, the stronger code is kept [S F1 #5]")
    if "O5" in codes:
        drop("O8", "O8 should not be used in combination with O5 [S F1 #8]")
        drop("O10", "O10 with O5 only if the reference variant was classified without functional "
                    "evidence, which cannot be verified [S F1 #16]")
    if "O7" in codes and "O8" in codes:
        drop("O8", "not permitted with O7 [S F1 #9]")
    if "O7" in codes and "O6" in codes and codes["O7"].points + codes["O6"].points > 4:
        drop("O6", "O6 + O7 limited to 4 points [S F1 #6]")
    sai = spliceai_scores(row)
    if "B6" in codes and sai and max(sai) >= 0.2:
        drop("B6", "B6 not applied when in silico tools predict a splicing impact [S F2 #4]")

    # O3 "should not be used if all other evidence suggests that a variant is likely
    # benign and application of O3 would move it into the VUS category" [S T1 O3]
    if "O3" in codes:
        without = sum(e.points for c, e in codes.items() if c != "O3")
        if class_from_points(without) in BEN_CLASSES and class_from_points(without + codes["O3"].points) == "VUS":
            drop("O3", "would move otherwise likely benign evidence into VUS [S T1 O3]")

    score = sum(e.points for e in codes.values())
    classification = class_from_points(score)

    # minimum of two evidence items in the same direction [G §6.3]
    n_onc = sum(1 for e in codes.values() if e.points > 0)
    n_ben = sum(1 for e in codes.values() if e.points < 0)
    if classification in ONC_CLASSES and n_onc < 2:
        notes.append(f"{classification} requires >=2 oncogenic codes, only {n_onc} -> VUS [G §6.3]")
        classification = "VUS"
    elif classification in BEN_CLASSES and n_ben < 2:
        notes.append(f"{classification} requires >=2 benign codes, only {n_ben} -> VUS [G §6.3]")
        classification = "VUS"

    ordered = sorted(codes.values(), key=lambda e: -e.points)
    return Result(classification, score, ordered, notes)


# =============================================================================
# CLI
# =============================================================================

OUT_COLUMNS = ["SVIG_UK_CLASS", "SVIG_UK_SCORE", "SVIG_UK_CODES", "SVIG_UK_NOTES"]

# =============================================================================
# Tier 4 — classification output with a two-row header (field names; then
# "description | source | version", as in Tier 2 and Tier 3)
# =============================================================================

CLASSIFIER_VERSION = "2.3"
GUIDELINE = "ACGS/SVIG-UK somatic oncogenicity guidelines v1.1 (26/05/2026)"
CODE_DEFINITIONS = {
    "O1": "Stand-alone: variant on the SVIG-UK canonical (likely) oncogenic variants list [S T3]",
    "O2": "Null variant in a tumour suppressor; PVS1 decision tree: NMD predicted from the stop position, NMD escape "
          "weighted by protein fraction (>10%), canonical splice sites assessed by exon skipping [S T1 O2]",
    "O3": "Absent (+2) or <=0.001% (+1) in gnomAD v4.1; CHIP allele-balance correction for DNMT3A, ASXL1, TET2, TP53 [S T1 O3]",
    "O4": "Enriched in a somatic database: GENIE unique patients in on-target tumour types (COSMIC fallback); needs O3 [S T1 O4]",
    "O5": "Different (likely) oncogenic amino acid change at the same residue with equivalent or higher REVEL [S T1 O5]",
    "O6": "Computational evidence of impact: REVEL >=0.7 or SpliceAI >=0.2 [S T1 O6]",
    "O7": "Cancer Hotspots recurrence at the residue / same amino acid change [S T1 O7]",
    "O8": "Missense variant in a missense-constrained gene (gnomAD missense Z >=3.09) [S T1 O8]",
    "O9": "In-frame length change in a functional domain, or final-exon truncation of an oncogene [S T1 O9]",
    "O10": "Functional studies: a calibrated assay where one covers the variant (<25% of wild-type "
           "activity is functionally abnormal), otherwise the SMART proxy of OncoKB mutation effect "
           "consistent with the gene role at supporting [S T1 O10]",
    "O11": "Tumour phenotype highly specific for the gene (not automated) [S T1 O11]",
    "B1": "Stand-alone >1% or strong >0.1% population frequency in gnomAD v4.1 [S T2 B1]",
    "B2": "Stand-alone: null variant in an oncogene without a tumour-suppressor role (classified VUS) [S T2 B2]",
    "B3": "Computational evidence of no impact: REVEL <=0.4 and SpliceAI <=0.1 [S T2 B3]",
    "B4": "Synonymous or intronic variant with no predicted splicing impact [S T2 B4]",
    "B5": "In-frame deletion/insertion in a repetitive region of unknown function [S T2 B5]",
    "B6": "Functional studies showing no damaging effect: a calibrated assay where one covers the variant "
           "and a normal boundary is recorded, otherwise the SMART proxy of OncoKB mutation effect Neutral at supporting [S T2 B6]",
    "B7": "Variant found with a different, established driver explaining the phenotype (not automated) [S T2 B7]",
}
NOT_AUTOMATED = {"O11", "B7"}
TIER4_ID_COLUMNS = ["Tumor_Sample_Barcode", "Hugo_Symbol", "NM_Transcript", "HGVSc", "HGVSp_Short",
                    "Chromosome", "Start_Position", "REF", "ALT", "Consequence", "EXON", "INTRON"]
TIER4_EVIDENCE_COLUMNS = ["gnomADe_AF", "gnomADg_AF", "REVEL", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
                          "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL", "VARIANT_IN_ONCOKB", "MUTATION_EFFECT",
                          "DOMAINS"]


REVIEW_STATUS = "Automated suggestion - review and confirmation by a clinical scientist required"
CLASS_BOUNDARIES = (10, 6, 0, -1, -7)  # lowest score of Oncogenic, LO, VUS; highest of LB (-1), Benign (-7)


def review_flags(result):
    """Reasons a scientist should look closely at this automated classification."""
    flags = []
    notes = " ".join(result.notes)
    if "not evaluated" in notes or "not provided" in notes:
        flags.append("evidence not evaluated (missing resource)")
    if "requires >=2" in notes:
        flags.append("minimum-two-codes rule changed the class")
    points = [e.points for e in result.codes]
    if any(p > 0 for p in points) and any(p < 0 for p in points):
        flags.append("conflicting oncogenic and benign evidence")
    if isinstance(result.score, int):
        if any(abs(result.score - b) <= 1 for b in CLASS_BOUNDARIES):
            flags.append("score within 1 point of a class boundary")
    if result.gene_role == "unknown":
        flags.append("gene mode of action unknown")
    if "O2 without transcript structure" in notes:
        flags.append("O2 assessed without transcript structure")
    if "functional assay" in notes:
        flags.append("O10/B6 from a functional assay: check the assay strength is the one agreed")
    if "gene-specific provision" in notes:
        flags.append("gene-specific guideline provision applies: confirm the class manually")
    if "B4 not applied" in notes:
        flags.append("synonymous variant at an exon edge: possible splice-motif disruption")
    if "CHIP" in notes or any("CHIP" in e.reason for e in result.codes):
        flags.append("CHIP gene: population frequency corrected")
    if "O4" in notes and "not evaluated" in notes:
        flags.append("O4 not evaluated")
    flags.append("O11 and B7 (tumour phenotype, other driver) not assessed")
    return flags


def tier4_header(smart_version=""):
    src = f"SMART SVIG-UK classifier v{CLASSIFIER_VERSION}"
    ver = f"{GUIDELINE}; SMART {smart_version}".strip("; ")
    meta = {c: f"Copied from Tier 1 | SMART annotation | {smart_version or 'see Tier 1'}"
            for c in TIER4_ID_COLUMNS + TIER4_EVIDENCE_COLUMNS}
    meta.update({
        "TUMOUR_TYPE_USED": f"Tumour type used for the classification (TUMOUR_TYPE column or --tumour-type; OncoTree "
                            f"code or broad group) | {src} | {ver}",
        "TUMOUR_GROUP": f"Tumour group derived from the tumour type (myeloid, lymphoid, haem, solid, unknown) | {src} | {ver}",
        "GENE_ROLE": f"Mode of action used for O2/B2/O10 (oncogene, tumour suppressor) and its source | {src} | {ver}",
        "SVIG_UK_CLASS": f"Final class: Oncogenic (>=10), Likely Oncogenic (6-9), VUS (0-5), Likely Benign (-1 to -6), "
                         f"Benign (<=-7); stand-alone O1/B1/B2 and the minimum-two-codes rule applied | {src} | {ver}",
        "SVIG_UK_SCORE": f"Sum of points of the codes applied (empty for stand-alone O1, B1, B2) | {src} | {ver}",
        "SVIG_UK_REVIEW_STATUS": f"The classification is an automated suggestion; SVIG-UK requires professional "
                                 f"judgement, so a clinical scientist must review and confirm it | {src} | {ver}",
        "SVIG_UK_REVIEW_FLAGS": f"Reasons to review this variant closely: missing resources, minimum-two-codes "
                                f"downgrades, conflicting evidence, score within 1 point of a class boundary, unknown "
                                f"gene role, CHIP correction, codes not automated | {src} | {ver}",
    })
    for code, definition in CODE_DEFINITIONS.items():
        value = "not automated in SMART" if code in NOT_AUTOMATED else "strength and points applied, empty if not met"
        meta[code] = f"{definition}; {value} | {src} | {ver}"
    meta["SVIG_UK_CODES"] = f"Every code applied with its strength, points and the evidence behind it | {src} | {ver}"
    meta["SVIG_UK_NOTES"] = (f"Codes removed by SVIG-UK combination rules, minimum-two-codes downgrades and resources "
                             f"not evaluated | {src} | {ver}")
    columns = (TIER4_ID_COLUMNS + ["TUMOUR_TYPE_USED", "TUMOUR_GROUP", "GENE_ROLE", "SVIG_UK_CLASS", "SVIG_UK_SCORE",
                                   "SVIG_UK_REVIEW_STATUS", "SVIG_UK_REVIEW_FLAGS"]
               + list(CODE_DEFINITIONS) + ["SVIG_UK_CODES", "SVIG_UK_NOTES"] + TIER4_EVIDENCE_COLUMNS)
    return columns, [meta[c] for c in columns]


def tier4_json(rows, results, smart_version, args):
    """Structured export of the Tier 4 classification, for import into review tools."""
    variants = []
    for row, result in zip(rows, results):
        variants.append({
            "variant": {c: row.get(c, "") for c in TIER4_ID_COLUMNS},
            "tumour_type": result.tumour_type, "tumour_group": result.tumour_group, "gene_role": result.gene_role,
            "classification": result.classification,
            "score": result.score if result.score != "" else None,
            "codes": [{"code": e.code, "strength": e.strength, "points": e.points, "evidence": e.reason}
                      for e in result.codes],
            "codes_not_met": [c for c in CODE_DEFINITIONS if c not in {e.code for e in result.codes}
                              and c not in NOT_AUTOMATED],
            "codes_not_automated": sorted(NOT_AUTOMATED),
            "notes": result.notes,
            "review_status": REVIEW_STATUS,
            "review_flags": review_flags(result),
        })
    return {
        "classification_system": GUIDELINE,
        "classifier": f"SMART SVIG-UK classifier v{CLASSIFIER_VERSION}",
        "smart_version": smart_version,
        "settings": {"tumour_type": args.tumour_type, "o4_on_target": getattr(args, "o4_on_target", "lineage")},
        "code_definitions": CODE_DEFINITIONS,
        "variants": variants,
    }


def tier4_row(row, result):
    out = {c: row.get(c, "") for c in TIER4_ID_COLUMNS + TIER4_EVIDENCE_COLUMNS}
    out.update(TUMOUR_TYPE_USED=result.tumour_type, TUMOUR_GROUP=result.tumour_group, GENE_ROLE=result.gene_role,
               SVIG_UK_CLASS=result.classification, SVIG_UK_SCORE=result.score,
               SVIG_UK_CODES=result.codes_str(), SVIG_UK_NOTES=result.notes_str(),
               SVIG_UK_REVIEW_STATUS=REVIEW_STATUS, SVIG_UK_REVIEW_FLAGS="; ".join(review_flags(result)))
    for code in CODE_DEFINITIONS:
        out[code] = ""
    for e in result.codes:
        out[e.code] = "stand-alone" if e.strength == "STANDALONE" else f"{e.strength} {e.points:+d}"
    return out


def load_resources(args, genes=None, accessions=None):
    res = Resources()
    here = os.path.dirname(os.path.abspath(__file__))
    canonical = args.canonical or os.path.join(here, "..", "reference_lists", "svig_uk_canonical_variants.tsv")
    res.load_canonical(canonical)
    vcep = getattr(args, "vcep_thresholds", None) or os.path.join(here, "..", "reference_lists",
                                                                   "svig_uk_vcep_population_thresholds.tsv")
    if os.path.exists(vcep):
        res.load_vcep_thresholds(vcep)
    gss = getattr(args, "gene_specific_splice", None) or os.path.join(
        here, "..", "reference_lists", "svig_uk_gene_specific_splice.tsv")
    if os.path.exists(gss):
        res.load_gene_specific_splice(gss)
    fa = getattr(args, "functional_assays", None) or os.path.join(
        here, "..", "reference_lists", "svig_uk_functional_assays.tsv")
    if os.path.exists(fa):
        res.load_functional_assays(fa, fa.replace("_assays.tsv", "_assay_scores.tsv"))
    refs = args.refs or ""

    def ref(path, default):
        p = path or os.path.join(refs, default)
        return p if p and os.path.exists(p) else None

    res.load_gene_roles(ref(args.cgc, "SVIG_UK/cgc_gene_roles.tsv"),
                        ref(args.oncokb_roles, "oncokb_gene_roles.tsv"))
    hs = ref(args.hotspots, "CancerHotSpots/cancerhotspots_counts.json")
    if hs:
        res.load_hotspots(hs)
    cons = ref(args.constraint, "gnomAD_constraints/loeuf_dataset_grch38.tsv.gz")
    if cons:
        res.load_constraint(cons)
    cosmic = ref(args.cosmic_o4, "SVIG_UK/cosmic_o4_counts.tsv.gz")
    if cosmic:
        res.cosmic = CosmicO4(cosmic, genes=genes)
    models = ref(getattr(args, "transcript_models", None), "TranscriptModels/ncbiRefSeq.hg38.txt.gz")
    if models:
        if accessions is not None:  # O5 maps canonical-list residues through their own transcripts
            accessions = set(accessions) | {_strip_version(c.get("transcript")) for c in res.canonical}
        res.load_transcript_models(models, accessions=accessions)
    genie = ref(getattr(args, "genie_o4", None), "GENIE/genie_o4_counts.tsv.gz")
    if genie:
        res.genie = GenieO4(genie, genes=genes)
        res.o4_level = getattr(args, "o4_on_target", None) or "lineage"
    oncotree = ref(getattr(args, "oncotree", None), "OncoTree/oncotree_2025_10_03.json")
    if oncotree:
        res.oncotree = OncoTree(oncotree)
    revel = ref(getattr(args, "revel", None), "REVEL/new_tabbed_revel_grch38.tsv.gz")
    if revel and os.path.exists(revel + ".tbi"):
        res.revel = RevelTabix(revel)
        res.load_o5_references(ref(getattr(args, "o5_references", None), "SVIG_UK/o5_reference_variants.tsv"))
    chip = ref(getattr(args, "chip_allele_balance", None), "SVIG_UK/gnomad_v41_chip_allele_balance.tsv.gz")
    if chip:
        res.load_chip_allele_balance(chip)
    fasta = ref(getattr(args, "genome_fasta", None), "liftover/hg38.fa")
    if fasta and os.path.exists(fasta + ".fai"):
        res.genome = Genome(fasta)
    return res


def read_maf(path):
    with open(path) as fh:
        lines = fh.readlines()
    preamble = [l for l in lines if l.startswith("#")]
    body = [l for l in lines if not l.startswith("#")]
    return preamble, list(csv.DictReader(body, delimiter="\t"))


def main(argv=None):
    ap = argparse.ArgumentParser(description="SVIG-UK oncogenicity classification of a SMART MAF")
    ap.add_argument("--maf", required=True)
    ap.add_argument("--out", help="copy of the input MAF with the SVIG_UK_* columns appended")
    ap.add_argument("--tier4", help="Tier 4 classification table with a two-row header")
    ap.add_argument("--tier4-json", help="structured JSON export of the Tier 4 classification")
    ap.add_argument("--refs", help="reference directory used to resolve default resource paths")
    ap.add_argument("--tumour-type", default="unknown",
                    help="OncoTree code (e.g. AML, MDS, DLBCLNOS) or myeloid | lymphoid | haem | solid | unknown "
                         "(used if the MAF has no TUMOUR_TYPE)")
    ap.add_argument("--canonical")
    ap.add_argument("--vcep-thresholds", help="ClinGen VCEP population thresholds for hereditary cancer genes "
                    "(default: reference_lists/svig_uk_vcep_population_thresholds.tsv)")
    ap.add_argument("--gene-specific-splice", help="gene-specific splice provisions "
                    "(default: reference_lists/svig_uk_gene_specific_splice.tsv)")
    ap.add_argument("--functional-assays", help="functional-assay parameters for O10/B6 "
                    "(default: reference_lists/svig_uk_functional_assays.tsv; the scores file is "
                    "found beside it)")
    ap.add_argument("--cgc")
    ap.add_argument("--oncokb-roles")
    ap.add_argument("--hotspots")
    ap.add_argument("--constraint")
    ap.add_argument("--cosmic-o4")
    ap.add_argument("--transcript-models", help="UCSC ncbiRefSeq table (exon structure for the O2/PVS1 tree)")
    ap.add_argument("--genie-o4", help="GENIE unique-patient counts (utils/build_genie_o4.py); "
                    "used for O4 in preference to COSMIC")
    ap.add_argument("--revel", help="tabix-indexed GRCh38 REVEL table (reference variant scores for O5)")
    ap.add_argument("--o5-references", help="extra (likely) oncogenic missense references for O5 "
                    "(utils/build_o5_references.py); the SVIG-UK canonical list is always used")
    ap.add_argument("--o4-on-target", choices=O4_ON_TARGET_LEVELS, default="lineage",
                    help="GENIE O4 on-target entries: code (OncoTree code and subtypes), tissue, "
                         "lineage (haematopoietic or solid; default) or all")
    ap.add_argument("--oncotree", help="OncoTree release JSON (resolves OncoTree codes given as tumour type)")
    ap.add_argument("--chip-allele-balance", help="gnomAD v4.1 allele-balance table for the O3 CHIP "
                    "correction (utils/build_gnomad_chip_ab.py)")
    ap.add_argument("--genome-fasta", help="faidx-indexed GRCh38 FASTA (reading frame after exon skipping)")
    args = ap.parse_args(argv)
    if not args.out and not args.tier4 and not args.tier4_json:
        ap.error("give --out, --tier4 and/or --tier4-json")

    preamble, rows = read_maf(args.maf)
    genes = {r.get("Hugo_Symbol", "") for r in rows}
    accessions = {_strip_version(r.get(c)) for r in rows
                  for c in ("Feature", "MANE_SELECT", "NM_Transcript") if r.get(c)}
    res = load_resources(args, genes=genes, accessions=accessions)

    results = [classify_row(r, res, args.tumour_type) for r in rows]

    if args.out:
        fieldnames = list(rows[0].keys()) if rows else []
        fieldnames += [c for c in OUT_COLUMNS if c not in fieldnames]
        with open(args.out, "w", newline="") as fh:
            fh.writelines(preamble)
            w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n",
                               extrasaction="ignore")
            w.writeheader()
            for r, result in zip(rows, results):
                w.writerow(dict(r, SVIG_UK_CLASS=result.classification, SVIG_UK_SCORE=result.score,
                                SVIG_UK_CODES=result.codes_str(), SVIG_UK_NOTES=result.notes_str()))
        print(f"classified {len(rows)} variants -> {args.out}", file=sys.stderr)

    smart_version = next((l.split()[-1] for l in preamble if l.startswith("#SMART_VERSION")), "")
    if args.tier4_json:
        with open(args.tier4_json, "w") as fh:
            json.dump(tier4_json(rows, results, smart_version, args), fh, indent=1)
        print(f"Tier 4 JSON: {len(rows)} variants -> {args.tier4_json}", file=sys.stderr)

    if args.tier4:
        columns, metadata = tier4_header(smart_version)
        with open(args.tier4, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(columns)
            w.writerow(metadata)
            for r, result in zip(rows, results):
                out = tier4_row(r, result)
                w.writerow([out[c] for c in columns])
        print(f"Tier 4: {len(rows)} variants -> {args.tier4}", file=sys.stderr)


if __name__ == "__main__":
    main()
