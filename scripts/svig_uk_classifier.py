#!/usr/bin/env python3

# =============================================================================
# svig_uk_classifier.py
# =============================================================================
#
# DESCRIPTION:
#   Automated oncogenicity classification of somatic small variants following
#   the ACGS / SVIG-UK 2025 guidelines ("Association for Clinical Genomic Science
#   guidelines for the classification of oncogenicity of somatic variants in
#   cancer", ratified 22/07/2025) and their Supplementary Material.
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
#     Implemented with proxy : O2 (PVS1 tree without exon-boundary and
#                              protein-length data), O4 (COSMIC all-tissue
#                              counts instead of GENIE), O9 and B5 (VEP domains),
#                              B2 (null variants in oncogenes only),
#                              O10 / B6 (OncoKB variant-level mutation effect)
#     Not automatable        : O5 (needs REVEL of the (likely) oncogenic
#                              reference variant), O11 / B7 (tumour phenotype),
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


NULL_CLASSES = {"truncating", "canonical_splice", "start_lost"}
SPLICE_REGION_TERMS = {"splice_region_variant", "splice_donor_5th_base_variant",
                       "splice_donor_region_variant", "splice_polypyrimidine_tract_variant"}


def exon_fraction(v):
    m = re.match(r"^\s*(\d+)(?:-\d+)?/(\d+)", str(v or ""))
    return (int(m.group(1)), int(m.group(2))) if m else (None, None)


def intronic_offset(hgvsc):
    m = re.search(r"c\.[*-]?\d+([+-]\d+)", str(hgvsc or ""))
    return int(m.group(1)) if m else None


def spliceai_scores(row):
    vals = [fnum(row.get(c)) for c in ("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
                                       "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")]
    vals = [v for v in vals if v is not None]
    return vals or None


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

    # --------------------------------------------------------------- loaders
    def load_canonical(self, path):
        with open(path) as fh:
            self.canonical = list(csv.DictReader(fh, delimiter="\t"))
        for c in self.canonical:
            c["_protein"] = normalise_protein(c.get("protein"))

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


def rule_b1(row):
    """B1 — high population frequency: >1% stand-alone benign, >0.1% strong [S T2 B1].
    Overall AF or the highest of the five continental groups, exomes or genomes."""
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


def rule_o3(row):
    """O3 — absent (+2) or very rare (+1) in gnomAD v4.1 [S T1 O3]."""
    afs = {k: fnum(row.get(f"{k[0]}_{k[1]}_AF")) for k in GNOMAD_SAMPLES}
    overall = [fnum(row.get("gnomADe_AF")), fnum(row.get("gnomADg_AF"))]
    if not any(v for v in list(afs.values()) + overall):
        return onc("O3", "MOD", "absent from gnomAD v4.1 exomes and genomes")
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


def rule_o2(row, vclass):
    """O2 — null variant in a TSG, PVS1-based strength [S T1 O2; S F3].
    SMART proxy: exon boundaries and protein length are not in the MAF, so the
    '>10% of protein' and 'last 50 bp of penultimate exon' branches cannot be
    evaluated; where the tree needs them, the lowest reachable strength is used."""
    exon_n, exon_total = exon_fraction(row.get("EXON"))
    if vclass == "truncating":
        cds = first_int(row.get("CDS_position"))
        if exon_n == 1 and cds is not None and cds <= 100:
            return onc("O2", "MOD", "stop within the first 100 bp of exon 1; re-initiation not assessable")
        if exon_n is None or exon_total is None:
            return onc("O2", "MOD", "exon position unknown; NMD not assessable")
        if exon_n == exon_total:
            return onc("O2", "MOD", "premature stop in the final exon (NMD escape; region criticality not assessable)")
        return onc("O2", "VSTR", "premature stop predicted to trigger NMD")
    if vclass == "canonical_splice":
        hgvsc = str(row.get("HGVSc") or "")
        sai = spliceai_scores(row)
        if re.search(r"\+2T>C", hgvsc) and (sai is None or max(sai) < 0.8):
            return None  # "+2T>C ... PVS1 should not be applied ... SpliceAI delta score <0.8"
        if re.search(r"\+2C>T", hgvsc):
            return None  # "+2C>T ... O2 is not applicable"
        return onc("O2", "MOD", "canonical +/-1,2 splice site; in-frame vs out-of-frame not assessable, lowest strength applied")
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
    if o6 is None and vclass != "canonical_splice" and sai_max is not None and sai_max >= 0.2:
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


def rule_b4(row, vclass, b3):
    """B4 — synonymous / deep intronic / UTR with no predicted splicing impact,
    supporting strength without an RNA assay [S T2 B4]."""
    if not b3 or not b3.reason.startswith("all SpliceAI"):
        return None
    if vclass == "synonymous" and not (terms(row.get("Consequence")) & SPLICE_REGION_TERMS):
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


def rule_o10_b6(row, roles, vclass):
    """O10 / B6 — SMART proxy from OncoKB variant-level mutation effect (see ONCOKB_* notes)."""
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
    group = tumour_group(row.get("TUMOUR_TYPE") or tumour_type)
    gene = row.get("Hugo_Symbol", "")
    protein = normalise_protein(row.get("HGVSp") or row.get("HGVSp_Short"))
    vclass = variant_class(row.get("Consequence"))
    roles, roles_src = res.roles(gene, group)
    notes = []

    # ---- stand-alone codes -------------------------------------------------
    canon = rule_o1(row, res, protein)
    if canon:
        note = f"SVIG-UK canonical list: {canon['gene']} {canon.get('protein') or canon.get('hgvsc')}"
        return Result(canon["assessment"], "", [Evidence("O1", "STANDALONE", 0, note)],
                      ["stand-alone O1 [G §6.3]"])

    b1 = rule_b1(row)
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
    add(rule_o3(row))
    if vclass in NULL_CLASSES and "TSG" in roles:
        add(rule_o2(row, vclass))
    add(rule_o9_b5(row, vclass, protein, roles))
    functional_region = bool(domain_sources(row) & FUNCTIONAL_DOMAIN_DBS)
    o7 = rule_o7(row, res, vclass, protein) if res.hotspots_single else None
    if res.hotspots_single:
        add(o7)
    else:
        notes.append("O7 not evaluated: Cancer Hotspots counts not provided")
    b3_or_o6 = rule_o6_b3(row, vclass, roles, functional_region or o7 is not None)
    add(b3_or_o6)
    add(rule_b4(row, vclass, b3_or_o6 if b3_or_o6 and b3_or_o6.code == "B3" else None))
    if res.mis_z:
        add(rule_o8(row, res, vclass))
    else:
        notes.append("O8 not evaluated: gnomAD constraint not provided")
    if res.cosmic is not None:
        add(rule_o4(row, res, vclass, protein))
    else:
        notes.append("O4 not evaluated: COSMIC counts not provided")
    add(rule_o10_b6(row, roles, vclass))

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
        drop("O9", "not permitted with O2 [S F1 #9]")
        drop("O10", "O10 with O2 only for missense variants [S F1 #13]")
    if "O9" in codes:
        drop("O6", "not permitted with O9 [S F1 #11]")
    if "O7" in codes and "O8" in codes:
        drop("O8", "not permitted with O7 [S F1 #8]")
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


def load_resources(args, genes=None):
    res = Resources()
    here = os.path.dirname(os.path.abspath(__file__))
    canonical = args.canonical or os.path.join(here, "..", "reference_lists", "svig_uk_canonical_variants.tsv")
    res.load_canonical(canonical)
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
    ap.add_argument("--out", required=True)
    ap.add_argument("--refs", help="reference directory used to resolve default resource paths")
    ap.add_argument("--tumour-type", default="unknown",
                    help="myeloid | lymphoid | haem | solid | unknown (used if the MAF has no TUMOUR_TYPE)")
    ap.add_argument("--canonical")
    ap.add_argument("--cgc")
    ap.add_argument("--oncokb-roles")
    ap.add_argument("--hotspots")
    ap.add_argument("--constraint")
    ap.add_argument("--cosmic-o4")
    args = ap.parse_args(argv)

    preamble, rows = read_maf(args.maf)
    genes = {r.get("Hugo_Symbol", "") for r in rows}
    res = load_resources(args, genes=genes)

    fieldnames = list(rows[0].keys()) if rows else []
    fieldnames += [c for c in OUT_COLUMNS if c not in fieldnames]
    with open(args.out, "w", newline="") as fh:
        fh.writelines(preamble)
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            result = classify_row(r, res, args.tumour_type)
            r["SVIG_UK_CLASS"] = result.classification
            r["SVIG_UK_SCORE"] = result.score
            r["SVIG_UK_CODES"] = result.codes_str()
            r["SVIG_UK_NOTES"] = result.notes_str()
            w.writerow(r)
    print(f"classified {len(rows)} variants -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
