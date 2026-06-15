#!/usr/bin/env python3
"""
SMART multi-evidence (SVIG-UK) classifier vs OncoKB-only vs geneticist (WGLS).

Goal: show how much SMART adds *over raw OncoKB* by combining the other evidence
columns (ClinVar, REVEL, SpliceAI, gnomAD, hotspots, LOEUF, consequence) using a
simplified, auditable implementation of the SVIG-UK / ACGS 2025 point system.

Input : output/output/Final_result_tier1.maf
Output: validation_report/
          - smart_multievidence_full.tsv   per-variant score, codes, both verdicts, WGLS
          - rescued_by_multievidence.tsv    OncoKB no-data cases that multi-evidence reclassifies
          - multievidence_summary.md        3-way comparison

NOTE ON SCOPE / LIMITATIONS (stated honestly in the report too):
  * Codes implemented from available columns: O3, O6, O7, O8, O10(proxy via OncoKB
    & ClinVar), B1, B3, B4, and a conservative O2 (LOF consequence).
  * NOT implemented (data not in this MAF): O4 (GENIE/COSMIC counts), O5, O9 repeat
    context, O11/B7 (tumour phenotype), B2/B5 (gene mode-of-action, repeat region).
  * O2 here does not check TSG status (no gene mode-of-action table), so it is applied
    conservatively and flagged. Thresholds follow the documented framework but the
    calibration is a v1 approximation, not the certified WGLS pipeline.
"""
import csv
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
INPUT = os.path.join(HERE, "..", "output", "output", "Final_result_tier1.maf")


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


LOF_CONS = ("frameshift_variant", "stop_gained", "splice_donor_variant",
            "splice_acceptor_variant", "start_lost")


def svig_score(r):
    """Return (points, list_of_codes) from automatable SVIG-UK evidence."""
    pts = 0
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

    # --- O10 / O1 proxy: OncoKB functional/clinical call -------------------
    if onc == "Oncogenic":
        pts += 10; codes.append("O10:OncoKB_Oncogenic(+10)")
    elif onc in ("Likely Oncogenic", "Resistance"):
        pts += 6; codes.append("O10:OncoKB_LikelyOnc(+6)")
    elif onc == "Likely Neutral":
        pts -= 4; codes.append("B6:OncoKB_LikelyNeutral(-4)")
    elif onc == "Neutral":
        pts -= 7; codes.append("B6:OncoKB_Neutral(-7)")

    # --- O10 / B6 via ClinVar ---------------------------------------------
    if cln in ("Pathogenic", "Pathogenic/Likely_pathogenic"):
        pts += 5; codes.append("O10:ClinVar_Path(+5)")
    elif cln == "Likely_pathogenic":
        pts += 3; codes.append("O10:ClinVar_LikelyPath(+3)")
    elif cln == "Benign":
        pts -= 5; codes.append("B6:ClinVar_Benign(-5)")
    elif cln in ("Likely_benign", "Benign/Likely_benign"):
        pts -= 3; codes.append("B6:ClinVar_LikelyBenign(-3)")
    if onc_clinvar and "Oncogenic" in onc_clinvar:
        pts += 4; codes.append("O10:ClinVar_ONC(+4)")

    # --- O2 conservative LOF (no TSG check available) ----------------------
    if impact == "HIGH" and any(c in cons for c in LOF_CONS):
        pts += 4; codes.append("O2:LOF_HIGH(+4,no-TSG-check)")

    # --- O7 hotspot --------------------------------------------------------
    if ch_hot or okb_hot:
        pts += 2; codes.append("O7:hotspot(+2)")

    # --- O3 rarity / B1 common (MAX_AF) ------------------------------------
    if max_af is not None:
        if max_af >= 0.01:
            pts -= 8; codes.append("B1:common_AF>=1%(-8)")
        elif max_af >= 0.001:
            pts -= 4; codes.append("B1:common_AF>=0.1%(-4)")
        elif max_af == 0:
            pts += 2; codes.append("O3:absent_gnomAD(+2)")
        elif max_af <= 1e-5:
            pts += 1; codes.append("O3:ultrarare(+1)")

    # --- O6 / B3 computational --------------------------------------------
    comp_onc = (revel is not None and revel >= 0.7) or (sai is not None and sai >= 0.2)
    comp_ben = (revel is not None and revel < 0.7) and (sai is None or sai < 0.1)
    if comp_onc:
        pts += 1; codes.append("O6:insilico_damaging(+1)")
    elif comp_ben:
        pts -= 1; codes.append("B3:insilico_benign(-1)")

    # --- O8 constraint -----------------------------------------------------
    if loeuf is not None and loeuf < 0.6:
        pts += 1; codes.append("O8:constrained_LOEUF(+1)")

    # --- B4 synonymous / low impact ---------------------------------------
    if "synonymous_variant" in cons and "missense" not in cons:
        pts -= 4; codes.append("B4:synonymous(-4)")

    return pts, codes


def score_to_5class(pts):
    if pts >= 10:  return "Oncogenic"
    if pts >= 6:   return "Likely Oncogenic"
    if pts >= 0:   return "VUS"
    if pts >= -6:  return "Likely Benign"
    return "Benign"


def to3(label):
    return {"Benign": 0, "Likely Benign": 0, "VUS": 1,
            "Likely Oncogenic": 2, "Oncogenic": 2}[label]


def wgls_pos(v): return v in ("Oncogenic", "Likely Oncogenic")
def label3_to_pos(c): return c == 2


def confusion(pairs):
    """pairs: list of (smart_positive_bool, wgls_positive_bool)."""
    tp = tn = fp = fn = 0
    for sp, wp in pairs:
        if sp and wp:           tp += 1
        elif not sp and not wp: tn += 1
        elif sp and not wp:     fp += 1
        else:                   fn += 1
    return tp, tn, fp, fn


def metrics(tp, tn, fp, fn):
    n = tp + tn + fp + fn
    sens = tp / (tp + fn) if (tp + fn) else float("nan")
    spec = tn / (tn + fp) if (tn + fp) else float("nan")
    ppv = tp / (tp + fp) if (tp + fp) else float("nan")
    acc = (tp + tn) / n if n else float("nan")
    f1 = 2 * ppv * sens / (ppv + sens) if (ppv + sens) else float("nan")
    return acc, sens, spec, ppv, f1


def weighted_kappa(triples, weight):
    """triples: list of (smart3, wgls3)."""
    cats = [0, 1, 2]
    n = len(triples)
    obs = Counter(triples)
    sm = Counter(a for a, _ in triples)
    wm = Counter(b for _, b in triples)
    num = den = 0.0
    for i in cats:
        for j in cats:
            w = weight(i, j)
            num += w * obs[(i, j)] / n
            den += w * (sm[i] / n) * (wm[j] / n)
    return 1 - num / den if den else float("nan")


def main():
    with open(INPUT) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    out_full = []
    for r in rows:
        pts, codes = svig_score(r)
        multi = score_to_5class(pts)
        out_full.append({
            "Hugo_Symbol": r.get("Hugo_Symbol", ""),
            "HGVSp_Short": r.get("HGVSp_Short", ""),
            "Consequence": r.get("Consequence", ""),
            "AGREEMENT": r.get("AGREEMENT", ""),
            "ONCOKB_VERDICT": r.get("MY_VERDICT", ""),
            "SMART_MULTI_VERDICT": multi,
            "SMART_SCORE": pts,
            "WGLS": r.get("WGLS", ""),
            "EVIDENCE_CODES": ";".join(codes),
        })

    # write full audit table (summary columns only)
    with open(os.path.join(HERE, "smart_multievidence_full.tsv"), "w", newline="") as g:
        w = csv.DictWriter(g, fieldnames=list(out_full[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(out_full)

    # write fully-annotated table: relevant comparison columns first, then all
    # the rest of the original MAF annotation.
    with open(INPUT) as f:
        orig_fields = csv.DictReader(f, delimiter="\t").fieldnames

    def multi_label(o):
        sp = label3_to_pos(to3(o["SMART_MULTI_VERDICT"]))
        wp = wgls_pos(o["WGLS"])
        if sp and wp:           return "TP"
        if (not sp) and (not wp): return "TN"
        if sp and not wp:       return "FP"
        return "FN"

    # leading columns (relevant for the comparison); the rest follow as-is
    lead = ["Hugo_Symbol", "HGVSp_Short", "Consequence", "AGREEMENT",
            "WGLS", "MY_VERDICT", "TP-TN_FN-FP",
            "SMART_MULTI_VERDICT", "SMART_SCORE", "TP-TN_FN-FP_MULTI",
            "EVIDENCE_CODES"]
    rest = [c for c in orig_fields if c not in lead]
    out_fields = lead + rest

    with open(os.path.join(HERE, "smart_multievidence_full_annotated.tsv"), "w", newline="") as g:
        w = csv.DictWriter(g, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        for r, o in zip(rows, out_full):
            r = dict(r)
            r["SMART_MULTI_VERDICT"] = o["SMART_MULTI_VERDICT"]
            r["SMART_SCORE"] = o["SMART_SCORE"]
            r["EVIDENCE_CODES"] = o["EVIDENCE_CODES"]
            r["TP-TN_FN-FP_MULTI"] = multi_label(o)
            w.writerow(r)

    # --- 3-way comparison ---------------------------------------------------
    okb_pairs = [(o["ONCOKB_VERDICT"] == "Oncogenic", wgls_pos(o["WGLS"])) for o in out_full]
    multi_pairs = [(label3_to_pos(to3(o["SMART_MULTI_VERDICT"])), wgls_pos(o["WGLS"])) for o in out_full]
    okb_trip = [(to3({"Oncogenic": "Oncogenic", "VUS": "VUS", "Benign": "Benign"}[o["ONCOKB_VERDICT"]]),
                 to3(o["WGLS"])) for o in out_full]
    multi_trip = [(to3(o["SMART_MULTI_VERDICT"]), to3(o["WGLS"])) for o in out_full]

    okb_cm = confusion(okb_pairs)
    multi_cm = confusion(multi_pairs)
    okb_m = metrics(*okb_cm)
    multi_m = metrics(*multi_cm)
    okb_wk = weighted_kappa(okb_trip, lambda i, j: abs(i - j))
    multi_wk = weighted_kappa(multi_trip, lambda i, j: abs(i - j))

    # --- what happens to the OncoKB no-data cases --------------------------
    nodata = [o for o in out_full if o["AGREEMENT"] == "ONCOKB_NO_DATA"]
    reclassified = [o for o in nodata if o["SMART_MULTI_VERDICT"] not in ("VUS",)]
    to_onc = [o for o in reclassified if o["SMART_MULTI_VERDICT"] in ("Oncogenic", "Likely Oncogenic")]
    to_ben = [o for o in reclassified if o["SMART_MULTI_VERDICT"] in ("Benign", "Likely Benign")]
    onc_correct = sum(1 for o in to_onc if wgls_pos(o["WGLS"]))
    ben_match = sum(1 for o in to_ben if to3(o["WGLS"]) == 0)
    ben_vs_vus = sum(1 for o in to_ben if o["WGLS"] == "VUS")
    ben_dangerous = sum(1 for o in to_ben if wgls_pos(o["WGLS"]))

    with open(os.path.join(HERE, "rescued_by_multievidence.tsv"), "w", newline="") as g:
        w = csv.DictWriter(g, fieldnames=list(out_full[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(reclassified)

    def fm(cm, m):
        return (f"TP={cm[0]} TN={cm[1]} FP={cm[2]} FN={cm[3]} | "
                f"acc={m[0]:.3f} sens={m[1]:.3f} spec={m[2]:.3f} ppv={m[3]:.3f} f1={m[4]:.3f}")

    md = f"""# 3-way comparison: OncoKB-only vs SMART multi-evidence vs Geneticist

{len(rows)} tier-1 variants. Binary task = oncogenic vs not (positive = geneticist
called Oncogenic/Likely Oncogenic).

## Headline

| Classifier | Acc | Sens | Spec | PPV | F1 | Weighted κ (3-class) |
|---|---|---|---|---|---|---|
| OncoKB-only (current `MY_VERDICT`) | {okb_m[0]:.3f} | {okb_m[1]:.3f} | {okb_m[2]:.3f} | {okb_m[3]:.3f} | {okb_m[4]:.3f} | {okb_wk:.3f} |
| SMART multi-evidence (SVIG-UK v1) | {multi_m[0]:.3f} | {multi_m[1]:.3f} | {multi_m[2]:.3f} | {multi_m[3]:.3f} | {multi_m[4]:.3f} | {multi_wk:.3f} |

```
OncoKB-only        : {fm(okb_cm, okb_m)}
SMART multi-evid.  : {fm(multi_cm, multi_m)}
```

## Where the multi-evidence value actually comes from

Of the **{len(nodata)} variants OncoKB had no data on** (OncoKB forces VUS), the
multi-evidence classifier moves **{len(reclassified)}** off VUS. The direction matters:

**Toward Oncogenic: {len(to_onc)}** — only {onc_correct} confirmed by the geneticist.
→ Multi-evidence does **NOT** recover the missed actionable variants. This is faithful
  to SVIG-UK: with OncoKB silent, high REVEL alone is only supporting (+1), and the
  strong discriminator O4 (GENIE/COSMIC somatic counts) is not in this MAF. **Adding O4
  is the single change most likely to recover oncogenic misses.**

**Toward Benign: {len(to_ben)}** — {ben_match} match the geneticist (benign/LB),
  {ben_vs_vus} the geneticist left as VUS (SMART more confident; defensible),
  and **{ben_dangerous} the geneticist called oncogenic (dangerous under-calls)**.

So the value SMART adds over raw OncoKB here is on the **benign side** — it correctly
de-prioritises ~{ben_match} non-actionable variants OncoKB couldn't judge, which is what
lifts specificity ({okb_m[2]:.3f}→{multi_m[2]:.3f}), PPV ({okb_m[3]:.3f}→{multi_m[3]:.3f})
and weighted kappa ({okb_wk:.3f}→{multi_wk:.3f}). The cost is {ben_dangerous} new
benign-but-actually-oncogenic calls to watch, and a small drop in sensitivity
({okb_m[1]:.3f}→{multi_m[1]:.3f}).

For the research goal (find actionable variants in genes the panel never looked at),
the take-home is: **multi-evidence is useful as a noise filter, but recovering the
actionable calls needs O4/somatic-frequency data.** See `rescued_by_multievidence.tsv`.

## Honesty / limitations

This is a **v1 approximation** of SVIG-UK, not the certified WGLS pipeline:
- Implemented codes: O2(conservative, no TSG check), O3, O6, O7, O8, O10(OncoKB+ClinVar), B1, B3, B4, B6.
- Missing (data not in MAF): O4 (GENIE/COSMIC counts — the strongest discriminator), O5, O9, O11, B2, B5, B7.
- Point thresholds follow the documented framework; calibration is not validated.

Per-variant scores and triggered evidence codes are in `smart_multievidence_full.tsv`
so every call is auditable.
"""
    with open(os.path.join(HERE, "multievidence_summary.md"), "w") as g:
        g.write(md)
    print(md)
    print(f"Wrote: multievidence_summary.md, smart_multievidence_full.tsv, "
          f"rescued_by_multievidence.tsv ({len(reclassified)} rows)")


if __name__ == "__main__":
    main()
