#!/usr/bin/env python3
"""
Recompute the SMART multi-evidence (SVIG-UK v1) score IN PLACE on the master table.

This is the audit script: it reads `smart_multievidence_full_annotated.tsv`, recomputes
EVIDENCE_CODES / SMART_SCORE / SMART_MULTI_VERDICT (and TP-TN_FN-FP_MULTI) from the raw
annotation columns already in that table, and writes them back into the SAME file.

Design for auditability:
  * Every fired evidence code carries its point contribution inline, e.g.
    "O10:OncoKB_Oncogenic(+10)". SMART_SCORE is the plain sum of those numbers, so any
    score can be checked by eye from EVIDENCE_CODES. The script asserts this equality
    for every row.
  * Before overwriting, it prints how many rows' SMART_SCORE / SMART_MULTI_VERDICT
    differ from whatever was previously stored, so you can see if the old values were off.
  * The file is rewritten atomically (temp file + os.replace); no extra tables are created.

Scoring scheme (the values we decided — simplified, auditable SVIG-UK / ACGS 2025):
  Implemented: O2(conservative LOF, no TSG check), O3, O6, O7, O8, O10(OncoKB+ClinVar),
               B1, B3, B4, B6.
  NOT implemented (data not in this table): O4 (GENIE/COSMIC somatic counts — the
               strongest discriminator), O5, O9, O11, B2, B5, B7.
  Thresholds follow the documented framework; calibration is a v1 approximation, not the
  certified WGLS pipeline.

Usage:
    python3 score_multievidence_inplace.py            # audit + rewrite in place
    python3 score_multievidence_inplace.py --dry-run  # audit only, do not write
"""
import csv
import os
import re
import sys
import tempfile
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
TABLE = os.path.join(HERE, "smart_multievidence_full_annotated.tsv")

# columns this script (re)writes; everything else is passed through untouched
COMPUTED = ["EVIDENCE_CODES", "SMART_SCORE", "SMART_MULTI_VERDICT",
            "TP-TN_FN-FP", "TP-TN_FN-FP_MULTI"]

LOF_CONS = ("frameshift_variant", "stop_gained", "splice_donor_variant",
            "splice_acceptor_variant", "start_lost")

ONC_POS = ("Oncogenic", "Likely Oncogenic")  # positive class for the binary task


def fnum(v):
    if v in ("", ".", "-", None):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def spliceai_max(r):
    vals = [fnum(r.get(c, "")) for c in
            ("SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
             "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")]
    vals = [v for v in vals if v is not None]
    return max(vals) if vals else None


def svig_score(r):
    """Return (points, [codes]). Each code string ends in '(+N)' or '(-N)'."""
    codes = []

    onc = r.get("ONCOGENIC", "")
    cln = r.get("ClinVar_CLNSIG", "")
    onc_clinvar = r.get("ClinVar_ONC", "")
    cons = r.get("Consequence", "")
    impact = r.get("IMPACT", "")
    revel = fnum(r.get("REVEL", ""))
    sai = spliceai_max(r)
    max_af = fnum(r.get("MAX_AF", ""))
    loeuf = fnum(r.get("LOEUF", ""))
    okb_hot = r.get("ONCOKB_HOTSPOT", "") == "TRUE"
    ch_hot = r.get("CancerHotspots_HOTSPOT", "") in ("1", "TRUE", "True")

    # --- O10 / B6 via OncoKB functional/clinical call ----------------------
    if onc == "Oncogenic":
        codes.append("O10:OncoKB_Oncogenic(+10)")
    elif onc in ("Likely Oncogenic", "Resistance"):
        codes.append("O10:OncoKB_LikelyOnc(+6)")
    elif onc == "Likely Neutral":
        codes.append("B6:OncoKB_LikelyNeutral(-4)")
    elif onc == "Neutral":
        codes.append("B6:OncoKB_Neutral(-7)")

    # --- O10 / B6 via ClinVar ----------------------------------------------
    if cln in ("Pathogenic", "Pathogenic/Likely_pathogenic"):
        codes.append("O10:ClinVar_Path(+5)")
    elif cln == "Likely_pathogenic":
        codes.append("O10:ClinVar_LikelyPath(+3)")
    elif cln == "Benign":
        codes.append("B6:ClinVar_Benign(-5)")
    elif cln in ("Likely_benign", "Benign/Likely_benign"):
        codes.append("B6:ClinVar_LikelyBenign(-3)")
    if onc_clinvar and "Oncogenic" in onc_clinvar:
        codes.append("O10:ClinVar_ONC(+4)")

    # --- O2 conservative LOF (no TSG mode-of-action check available) --------
    if impact == "HIGH" and any(c in cons for c in LOF_CONS):
        codes.append("O2:LOF_HIGH(+4,no-TSG-check)")

    # --- O7 hotspot ---------------------------------------------------------
    if ch_hot or okb_hot:
        codes.append("O7:hotspot(+2)")

    # --- O3 rarity / B1 common (population AF) ------------------------------
    if max_af is not None:
        if max_af >= 0.01:
            codes.append("B1:common_AF>=1%(-8)")
        elif max_af >= 0.001:
            codes.append("B1:common_AF>=0.1%(-4)")
        elif max_af == 0:
            codes.append("O3:absent_gnomAD(+2)")
        elif max_af <= 1e-5:
            codes.append("O3:ultrarare(+1)")

    # --- O6 / B3 computational ---------------------------------------------
    comp_onc = (revel is not None and revel >= 0.7) or (sai is not None and sai >= 0.2)
    comp_ben = (revel is not None and revel < 0.7) and (sai is None or sai < 0.1)
    if comp_onc:
        codes.append("O6:insilico_damaging(+1)")
    elif comp_ben:
        codes.append("B3:insilico_benign(-1)")

    # --- O8 constraint ------------------------------------------------------
    if loeuf is not None and loeuf < 0.6:
        codes.append("O8:constrained_LOEUF(+1)")

    # --- B4 synonymous / low impact ----------------------------------------
    if "synonymous_variant" in cons and "missense" not in cons:
        codes.append("B4:synonymous(-4)")

    pts = sum(int(m) for c in codes for m in re.findall(r"\(([+-]\d+)", c))
    return pts, codes


def score_to_5class(pts):
    if pts >= 10:  return "Oncogenic"
    if pts >= 6:   return "Likely Oncogenic"
    if pts >= 0:   return "VUS"
    if pts >= -6:  return "Likely Benign"
    return "Benign"


def tp_label(verdict, wgls):
    """Confusion label for one verdict vs the geneticist (WGLS).

    Positive = (likely) oncogenic on BOTH sides. Works for the OncoKB-only verdict
    (MY_VERDICT, whose only oncogenic value is 'Oncogenic') and for the multi-evidence
    verdict (SMART_MULTI_VERDICT, which can also be 'Likely Oncogenic').
    """
    sp = verdict in ONC_POS
    wp = wgls in ONC_POS
    if sp and wp:           return "TP"
    if (not sp) and (not wp): return "TN"
    if sp and not wp:       return "FP"
    return "FN"


def okb_3class(v):
    """OncoKB-only verdict (MY_VERDICT) collapsed to ONC / VUS / BEN."""
    return {"Oncogenic": "ONC", "VUS": "VUS", "Benign": "BEN"}.get(v, "VUS")


def wgls_3class(v):
    """Geneticist verdict (WGLS) collapsed to ONC / VUS / BEN."""
    return {"Oncogenic": "ONC", "Likely Oncogenic": "ONC", "VUS": "VUS",
            "Likely Benign": "BEN", "Benign": "BEN"}.get(v, "VUS")


def breakdown_3class(title, pred_name, c3):
    """Print an ONC/VUS/BEN breakdown of a classifier vs the geneticist.

    c3 is a Counter keyed by (predicted_3class, geneticist_3class). For each class it
    shows how many the geneticist assigned, how many the classifier assigned, and on how
    many they agree (= the classifier 'got right', taking the geneticist as truth).
    """
    n = sum(c3.values())
    print(f"{title} — 3-class breakdown vs geneticist:")
    print(f"  class   geneticist   {pred_name:6}   agree ({pred_name} correct)")
    for cls in ("ONC", "VUS", "BEN"):
        gen = sum(v for (p, w), v in c3.items() if w == cls)   # geneticist total
        pred = sum(v for (p, w), v in c3.items() if p == cls)  # classifier total
        agree = c3[(cls, cls)]                                  # both call this class
        pct = f"{100*agree/gen:.1f}%" if gen else "n/a"
        print(f"  {cls:5}   {gen:9}   {pred:6}   {agree:5}  ({pct} of geneticist {cls})")
    diag = sum(c3[(c, c)] for c in ("ONC", "VUS", "BEN"))
    overall = f"{100*diag/n:.1f}%" if n else "n/a"
    print(f"  total   {n:9}   {n:6}   {diag:5}  ({overall} overall agreement)")
    print()


def _review_strength(r):
    """Heuristic 'how worth a second look' score for a SMART-oncogenic / geneticist-not
    variant. Higher = more independent, harder evidence the geneticist should re-review."""
    s = 0
    if r.get("ONCOGENIC") == "Oncogenic":               s += 3   # OncoKB definitive
    elif r.get("ONCOGENIC") == "Likely Oncogenic":      s += 1
    cln = r.get("ClinVar_CLNSIG", "")
    if cln in ("Pathogenic", "Pathogenic/Likely_pathogenic"): s += 3
    elif cln == "Likely_pathogenic":                    s += 2
    if r.get("ONCOKB_HOTSPOT") == "TRUE" or \
       r.get("CancerHotspots_HOTSPOT") in ("1", "TRUE", "True"): s += 2
    if r.get("MUTATION_EFFECT", "") in ("Gain-of-function", "Loss-of-function"): s += 2
    return s


def review_candidates(rows, top=20):
    """Print variants SMART called oncogenic but the geneticist did not (OncoKB-only FP),
    deduplicated and ranked by evidence strength. These are the candidates the geneticist
    may want to re-review — strongest first. NOT claims of error: SMART/OncoKB is
    tumour-agnostic while WGLS is tumour-specific, so a 'general' oncogenic call can be a
    legitimate VUS in this tumour. ClinVar (Likely) Pathogenic is the most independent flag.
    """
    seen, fp = set(), []
    for r in rows:
        if r.get("TP-TN_FN-FP") != "FP":
            continue
        k = (r.get("Hugo_Symbol", ""), r.get("HGVSp_Short", ""), r.get("Start_Position", ""))
        if k in seen:
            continue
        seen.add(k)
        fp.append(r)
    fp.sort(key=_review_strength, reverse=True)

    print(f"Candidates for geneticist re-review (SMART oncogenic, geneticist not) — "
          f"{len(fp)} unique, strongest first:")
    print(f"  {'gene':7} {'variant':14} {'WGLS':14} {'OncoKB':16} "
          f"{'ClinVar':24} hot effect")
    for r in fp[:top]:
        hot = "Y" if (r.get("ONCOKB_HOTSPOT") == "TRUE" or
                      r.get("CancerHotspots_HOTSPOT") in ("1", "TRUE", "True")) else "-"
        print(f"  {r.get('Hugo_Symbol',''):7} {r.get('HGVSp_Short',''):14} "
              f"{r.get('WGLS',''):14} {r.get('ONCOGENIC',''):16} "
              f"{(r.get('ClinVar_CLNSIG','') or '-'):24} {hot}   "
              f"{r.get('MUTATION_EFFECT','')}")
    if len(fp) > top:
        print(f"  ... ({len(fp)-top} more)")
    path = [r for r in fp if r.get("ClinVar_CLNSIG", "") in
            ("Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic")]
    print(f"  -> {len(path)} have independent ClinVar (Likely) Pathogenic support "
          f"(highest priority): " + ", ".join(
              f"{r.get('Hugo_Symbol','')} {r.get('HGVSp_Short','') or '(splice)'}" for r in path))
    print()


def confusion_columns(r):
    """Return (okb_label, multi_label): the two TP/TN/FP/FN columns for this row.

    - TP-TN_FN-FP       : OncoKB-only call  (MY_VERDICT)         vs WGLS
    - TP-TN_FN-FP_MULTI : SMART multi-evidence (SMART_MULTI_VERDICT) vs WGLS
    """
    wgls = r.get("WGLS", "")
    return tp_label(r.get("MY_VERDICT", ""), wgls), tp_label(r.get("SMART_MULTI_VERDICT", ""), wgls)


def main():
    dry = "--dry-run" in sys.argv[1:]

    with open(TABLE, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames
        rows = list(reader)

    missing = [c for c in COMPUTED if c not in fields]
    if missing:
        sys.exit(f"ERROR: table is missing computed columns {missing}; aborting.")

    changed_score = changed_verdict = changed_okb = changed_multi = 0  # noqa: E741
    okb_counts = {"TP": 0, "TN": 0, "FP": 0, "FN": 0}
    okb_eval_counts = {"TP": 0, "TN": 0, "FP": 0, "FN": 0}  # OncoKB rows that have data
    multi_counts = {"TP": 0, "TN": 0, "FP": 0, "FN": 0}
    okb_3x3 = Counter()    # (okb_3class, wgls_3class) over ALL rows
    eval_3x3 = Counter()   # (okb_3class, wgls_3class) over evaluable rows only
    multi_3x3 = Counter()  # (multi_3class, wgls_3class) over ALL rows
    for r in rows:
        old_score = r.get("SMART_SCORE", "")
        old_verdict = r.get("SMART_MULTI_VERDICT", "")
        old_okb = r.get("TP-TN_FN-FP", "")
        old_multi = r.get("TP-TN_FN-FP_MULTI", "")

        pts, codes = svig_score(r)
        verdict = score_to_5class(pts)
        codes_str = ";".join(codes)

        # SMART_MULTI_VERDICT feeds the multi-evidence confusion label, so set it on the
        # row first, then derive both TP/TN/FP/FN columns from the (updated) row.
        r["SMART_MULTI_VERDICT"] = verdict
        okb_label, multi_label = confusion_columns(r)

        # audit: sum of code points must equal the score (guards the scheme itself)
        check = sum(int(m) for c in codes for m in re.findall(r"\(([+-]\d+)", c))
        assert check == pts, f"code/score mismatch: {codes_str} -> {check} != {pts}"

        if old_score != str(pts):
            changed_score += 1
        if old_verdict != verdict:
            changed_verdict += 1
        if old_okb != okb_label:
            changed_okb += 1
        if old_multi != multi_label:
            changed_multi += 1

        okb_counts[okb_label] += 1
        multi_counts[multi_label] += 1
        w3 = wgls_3class(r.get("WGLS", ""))
        okb_3x3[(okb_3class(r.get("MY_VERDICT", "")), w3)] += 1
        multi_3x3[(wgls_3class(verdict), w3)] += 1  # multi verdict shares WGLS label space
        if r.get("AGREEMENT", "") != "ONCOKB_NO_DATA":  # OncoKB had a call here
            okb_eval_counts[okb_label] += 1
            eval_3x3[(okb_3class(r.get("MY_VERDICT", "")), w3)] += 1

        r["EVIDENCE_CODES"] = codes_str
        r["SMART_SCORE"] = str(pts)
        r["TP-TN_FN-FP"] = okb_label
        r["TP-TN_FN-FP_MULTI"] = multi_label

    def line(c):
        n = sum(c.values())
        sens = c["TP"] / (c["TP"] + c["FN"]) if (c["TP"] + c["FN"]) else float("nan")
        spec = c["TN"] / (c["TN"] + c["FP"]) if (c["TN"] + c["FP"]) else float("nan")
        return (f"TP={c['TP']} TN={c['TN']} FP={c['FP']} FN={c['FN']} | "
                f"sens={sens:.3f} spec={spec:.3f} (N={n})")

    print(f"Rows audited: {len(rows)}")
    print(f"  SMART_SCORE differing from stored value : {changed_score}")
    print(f"  SMART_MULTI_VERDICT differing from stored: {changed_verdict}")
    print(f"  TP-TN_FN-FP differing from stored        : {changed_okb}")
    print(f"  TP-TN_FN-FP_MULTI differing from stored  : {changed_multi}")
    print("  EVIDENCE_CODES sum == SMART_SCORE for every row: OK")
    print()
    print("Confusion vs geneticist (WGLS), positive = (Likely) Oncogenic:")
    print(f"  OncoKB-only, all variants        : {line(okb_counts)}")
    print(f"  OncoKB-only, evaluable (has data): {line(okb_eval_counts)}")
    print(f"  Multi-evidence, all variants     : {line(multi_counts)}")
    print()
    breakdown_3class("OncoKB-only, all variants", "OncoKB", okb_3x3)
    breakdown_3class("OncoKB-only, evaluable (has data)", "OncoKB", eval_3x3)
    breakdown_3class("Multi-evidence, all variants", "SMART", multi_3x3)
    review_candidates(rows)

    if dry:
        print("\n--dry-run: table NOT modified.")
        return

    # atomic in-place rewrite (temp file in same dir, then replace)
    fd, tmp = tempfile.mkstemp(dir=HERE, suffix=".tmp")
    with os.fdopen(fd, "w", newline="") as g:
        w = csv.DictWriter(g, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, TABLE)
    print(f"\nRewrote in place: {os.path.basename(TABLE)} "
          f"({len(rows)} rows, {len(fields)} cols unchanged)")


if __name__ == "__main__":
    main()
