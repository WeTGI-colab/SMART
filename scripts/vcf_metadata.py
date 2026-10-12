#!/usr/bin/env python3
"""
vcf_metadata.py — take field metadata from the VCFs themselves rather than a hand-kept list.

A VCF documents its own fields. The annotated VCF carries one ##INFO header per field that VEP
and its plugins added, with the resource file each came from; the OncoKB VCF carries one per
field it added; and the headers written by the variant caller survive annotation untouched.
SMART used to read exactly one of those headers (CSQ) and describe every other column from
Config.yaml, which is maintained by hand and therefore drifts: a field was attributed to Cancer
Hotspots when the caller had written it, and another to Manta when it came from DRAGEN.

Metadata taken from the file cannot drift, because it is written by whatever produced the data,
in the run that produced it. So:

    description   the VCF wins wherever it has one; Config.yaml fills the gaps
    source        the VCF wins
    version       the VCF wins, and for every OncoKB field it is the data version observed
                  in this run — an analysis cannot span two OncoKB releases

Columns that no VCF describes — the MAF standard columns, the values SMART computes, and the
OncoKB MafAnnotator output, which is a table and not a VCF — stay with Config.yaml.

    python3 scripts/vcf_metadata.py --annotated s_annotated.vcf --oncokb s.oncoKB.vcf \
        --input s.vcf.gz --out s.vcf_metadata.json
    python3 scripts/vcf_metadata.py --sidecar s.vcf_metadata.json --summary
"""
import argparse
import gzip
import json
import os
import re
import sys

# ``##INFO=<ID=x,Number=1,Type=String,Description="...">`` and the same shape for ##FORMAT
HEADER_RE = re.compile(
    r'^##(?P<kind>INFO|FORMAT)=<ID=(?P<id>[^,]+),'
    r'(?:Number=(?P<number>[^,]*),)?(?:Type=(?P<type>[^,]*),)?'
    r'.*?Description="(?P<description>(?:[^"\\]|\\.)*)"', re.S)
# VEP writes the subfields of CSQ as a pipe-delimited list inside its own Description
CSQ_RE = re.compile(r'^##INFO=<ID=CSQ,.*Format: (?P<fields>[^"]+)"')
# ``##VEP="v114.0" API="v114" time="..." cache="/refs/homo_sapiens_merged/114_GRCh38"``
VEP_RE = re.compile(r'^##VEP="(?P<version>[^"]+)"(?:.*?cache="(?P<cache>[^"]+)")?')
# a plugin announces itself as ``##SpliceAI_pred_DS_AG=SpliceAI (v1.3) predicted effect ...``
PLUGIN_RE = re.compile(r"^##(?P<id>[A-Za-z][\w.]*)=(?P<text>[^<].*)$")
# a custom annotation's Description names the file it came from
FILE_RE = re.compile(r"(?:from )?(?P<path>[^\s\"]*?(?P<base>[\w.\-]+\.(?:vcf|bed|tsv|txt|gz))\b)")
SIMPLE_HEADERS = ("fileformat", "reference", "source")


def _open(path):
    return gzip.open(path, "rt", errors="replace") if str(path).endswith(".gz") \
        else open(path, errors="replace")


def read_headers(path):
    """Every ##INFO/##FORMAT entry and the run-level lines, from one VCF.

    Returns (fields, run, csq_fields). Reading stops at the first data line: the headers are
    all we need and a VCF can be very large.
    """
    fields, run, csq = {}, {}, []
    if not path or not os.path.exists(path):
        return fields, run, csq
    with _open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            m = CSQ_RE.match(line)
            if m:
                csq = [f.strip() for f in m.group("fields").split("|") if f.strip()]
                continue
            m = HEADER_RE.match(line)
            if m:
                g = m.groupdict()
                fields[g["id"]] = {"kind": g["kind"], "number": g["number"], "type": g["type"],
                                   "description": g["description"].replace('\\"', '"').strip()}
                continue
            m = VEP_RE.match(line)
            if m:
                run["vep_version"] = m.group("version")
                if m.group("cache"):
                    run["vep_cache"] = m.group("cache")
                continue
            m = PLUGIN_RE.match(line)
            if m:
                key, text = m.group("id"), m.group("text").strip()
                if key.lower() in SIMPLE_HEADERS:
                    run[key.lower()] = text
                else:
                    # a plugin's self-description doubles as the description of its field
                    run.setdefault("_plugin_notes", {})[key] = text
    return fields, run, csq


def _source_from(description, default):
    """The resource a custom annotation names in its own Description, if it names one."""
    m = FILE_RE.search(description or "")
    return m.group("base") if m else default


def harvest(annotated=None, oncokb=None, input_vcf=None, oncokb_version=None):
    """Field metadata for one sample, with each field tagged by where it came from.

    origin is what decides precedence later: "input" and "vep" and "oncokb" fields are
    described by the VCF, everything else falls through to Config.yaml.
    """
    in_fields, in_run, _ = read_headers(input_vcf)
    an_fields, an_run, csq = read_headers(annotated)
    ok_fields, ok_run, _ = read_headers(oncokb)

    run = {**in_run, **an_run, **ok_run}
    run.pop("_plugin_notes", None)
    notes = {**(in_run.get("_plugin_notes") or {}), **(an_run.get("_plugin_notes") or {}),
             **(ok_run.get("_plugin_notes") or {})}
    if input_vcf:
        run["input_caller"] = in_run.get("source", "unknown")
    if oncokb_version:
        run["oncokb_data_version"] = oncokb_version
    vep = run.get("vep_version", "")

    fields = {}

    def put(name, description, source, version, origin):
        fields[name] = {"description": description, "source": source,
                        "version": version, "origin": origin}

    # 1. the caller's own fields. Only INFO entries become a column of the same name; FORMAT
    #    entries are expanded and renamed by post_analysis (AF becomes VAF), so they must not
    #    claim the column name — doing so would attribute VEP's AF to the caller.
    caller = f"Input VCF ({run.get('input_caller', 'variant caller')})"
    for name, d in in_fields.items():
        if d["kind"] == "INFO":
            put(name, d["description"], caller, "Not versioned by the caller", "input")
    format_fields = {n: {**d, "source": caller} for n, d in in_fields.items()
                     if d["kind"] == "FORMAT"}

    # 2. fields VEP and its plugins added. A custom annotation names its resource file in its
    #    own Description, which is the honest version string for that run.
    for name, d in an_fields.items():
        if name in in_fields or name == "CSQ":
            continue
        resource = _source_from(d["description"], None)
        if resource:
            # VEP names a custom annotation's fields after the annotation itself
            # (ClinVar_CLNSIG, COSMIC_CNT), so the prefix is the readable source and the
            # file it names is the version actually used in this run.
            put(name, d["description"], name.split("_")[0], resource, "vep")
        else:
            put(name, d["description"], "VEP", vep, "vep")

    # 3. CSQ subfields with no header of their own are core VEP output: no description to take,
    #    but the version is still a fact of this run.
    for name in csq:
        if name in fields:
            continue
        note = notes.get(name)
        if note:
            # "SpliceAI (v1.3) predicted effect on splicing..." names the plugin and its version
            m = re.match(r"(?P<plugin>[A-Za-z][\w.\-]*)\s*\((?P<ver>v?[\d.]+)\)", note)
            src = m.group("plugin") if m else _source_from(note, "VEP")
            ver = m.group("ver") if m else vep
            put(name, note, src, ver, "vep")
        else:
            put(name, "", "VEP", vep, "vep")

    # 4. OncoKB. Its VCF is the annotated VCF with its own fields added, so only the entries
    #    that were not already there belong to it; the rest are VEP's and keep VEP's provenance.
    #    Its data version is a property of the run, not of the tool.
    for name, d in ok_fields.items():
        if name == "CSQ" or name in in_fields or name in an_fields:
            continue
        put(name, d["description"], "OncoKB",
            run.get("oncokb_data_version", "see ONCOKB_DATA_VERSION"), "oncokb")

    return {"fields": fields, "run": run,
            "input_field_ids": sorted(in_fields),
            # only INFO entries can overwrite an annotation: vcf2table applies them after the
            # CSQ fields, while FORMAT entries are expanded separately and renamed (AF -> VAF)
            "input_info_ids": sorted(k for k, v in in_fields.items() if v["kind"] == "INFO"),
            "input_format_fields": {n: d["description"] for n, d in format_fields.items()},
            "csq_fields": csq}


def merge(sidecars):
    """Combine the per-sample sidecars of one run into the metadata for the merged tables.

    A run is supposed to be one pipeline, one set of resources. Where two samples disagree —
    a different VEP version, a different ClinVar file, a different OncoKB release — the run was
    not homogeneous and the merged header would be describing some of its own rows wrongly, so
    the disagreement is reported rather than silently resolved.
    """
    fields, run, disagreements = {}, {}, []
    info_ids, csq, fmt = set(), [], {}
    for path in sidecars:
        with open(path) as fh:
            m = json.load(fh)
        sample = os.path.basename(path).split(".")[0]
        for name, d in m.get("fields", {}).items():
            seen = fields.get(name)
            if seen and (seen["source"], seen["version"]) != (d["source"], d["version"]):
                disagreements.append((name, f"{seen['source']}/{seen['version']}",
                                      f"{d['source']}/{d['version']}", sample))
            fields.setdefault(name, d)
        for k, v in m.get("run", {}).items():
            if k in run and run[k] != v:
                disagreements.append((f"run.{k}", str(run[k]), str(v), sample))
            run.setdefault(k, v)
        info_ids |= set(m.get("input_info_ids") or [])
        fmt.update(m.get("input_format_fields") or {})
        csq = csq or m.get("csq_fields") or []
    return {"fields": fields, "run": run, "input_info_ids": sorted(info_ids),
            "input_field_ids": sorted(info_ids | set(fmt)), "input_format_fields": fmt,
            "csq_fields": csq, "disagreements": disagreements}


def check(meta, table_columns=None, liftover=False, oncokb_versions=None):
    """The things that can go wrong, in the order they matter.

    Returns a list of (level, message). "error" means a value in the output is wrong;
    "warning" means it is incomplete or its provenance is doubtful.
    """
    out = []
    fields, run = meta["fields"], meta["run"]
    csq = set(meta.get("csq_fields") or [])
    inp = set(meta.get("input_field_ids") or [])
    inp_info = set(meta.get("input_info_ids") or inp)

    # A field the caller wrote and a field an annotation added can share a name, and the INFO
    # value is applied last — so the annotation's value is silently replaced. That is a wrong
    # value, not a missing one.
    clash = sorted(inp_info & csq)
    if clash:
        out.append(("error", "name collision between the input VCF and VEP/plugin fields: "
                             + ", ".join(clash)
                             + ". The caller's value overwrites the annotation in the table."))

    if table_columns is not None:
        # only INFO entries are expected under their own name; FORMAT entries are expanded
        # and renamed by post_analysis, so their absence here is not a loss
        lost = sorted(i for i in inp_info if i not in set(table_columns))
        if lost:
            out.append(("warning", f"{len(lost)} field(s) declared in the input VCF did not "
                                   f"reach the table: {', '.join(lost[:12])}"
                                   + (" …" if len(lost) > 12 else "")))

    if liftover:
        out.append(("warning", "liftover enabled — INFO fields written by the variant caller "
                               "may not survive GATK LiftoverVcf, so field provenance taken "
                               "from the annotated VCF may be incomplete for this run."))

    seen = sorted({v for v in (oncokb_versions or []) if v})
    if len(seen) > 1:
        out.append(("error", "more than one OncoKB data version in a single run "
                             f"({', '.join(seen)}). An analysis cannot span two releases."))

    if not run.get("vep_version"):
        out.append(("warning", "no ##VEP header in the annotated VCF: VEP version unknown."))

    for name, a, b, sample in meta.get("disagreements", [])[:10]:
        out.append(("error", f"samples in this run disagree on {name}: {a} vs {b} "
                             f"(first seen differing at {sample}). The merged tables would "
                             f"describe some of their own rows wrongly."))
    return out


def _same_source(a, b):
    """'Cancer Hotspots' and 'CancerHotspots' name the same thing."""
    norm = lambda x: re.sub(r"[^a-z0-9]", "", str(x or "").lower())
    na, nb = norm(a), norm(b)
    return bool(na) and bool(nb) and (na in nb or nb in na)


def _is_generic_vep(source):
    """A core CSQ field has no header of its own, so all the file says is 'VEP delivered this'.
    That is weaker than Config.yaml naming the upstream resource (gnomAD, 1000G, CADD), so it
    does not override it and is not reported as a contradiction."""
    return str(source).strip() == "VEP"


def conflicts(meta, field_config):
    """Where Config.yaml claims a source the file contradicts. The file is right."""
    out = []
    for name, d in meta["fields"].items():
        cfg = field_config.get(name)
        if not cfg or _is_generic_vep(d["source"]):
            continue
        claimed = str(cfg.get("source", "")).strip()
        if claimed and not _same_source(claimed, d["source"]):
            out.append((name, claimed, d["source"]))
    return out


def resolve(name, meta, field_config, field_patterns=None, oncokb_version=None):
    """Metadata for one column, applying the agreed precedence.

    A field the VCFs describe is described by them. Everything else keeps Config.yaml — except
    the OncoKB MafAnnotator columns, whose version is the data version of this run.
    """
    harvested = meta["fields"].get(name)
    cfg = dict(field_config.get(name) or {})
    if not cfg and field_patterns:
        import fnmatch
        for p in field_patterns:
            if fnmatch.fnmatch(name, p.get("pattern", "")):
                cfg = {k: v for k, v in p.items() if k != "pattern"}
                break
    if harvested:
        source = harvested["source"]
        if _is_generic_vep(source) and cfg.get("source"):
            # the run still fixes the version; the config still names the upstream resource
            source = cfg["source"]
        return {"description": harvested["description"] or cfg.get("description", ""),
                "source": source, "version": harvested["version"],
                "tier": cfg.get("tier"), "origin": harvested["origin"]}
    if name.startswith("ONCOKB_") and oncokb_version:
        cfg["version"] = oncokb_version
        cfg["origin"] = "oncokb-maf"
    return cfg


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotated", help="VEP-annotated VCF")
    ap.add_argument("--oncokb", help="OncoKB-annotated VCF")
    ap.add_argument("--input", dest="input_vcf", help="the VCF as it arrived from the caller")
    ap.add_argument("--oncokb-version", help="OncoKB data version observed in this run")
    ap.add_argument("--out", help="write the harvested metadata here as JSON")
    ap.add_argument("--sidecar", help="read a previously written sidecar instead of harvesting")
    ap.add_argument("--merge", nargs="+", help="combine the per-sample sidecars of one run")
    ap.add_argument("--summary", action="store_true", help="print what was harvested")
    a = ap.parse_args(argv)

    if a.merge:
        meta = merge(a.merge)
    elif a.sidecar:
        meta = json.load(open(a.sidecar))
    elif a.annotated or a.oncokb or a.input_vcf:
        meta = harvest(a.annotated, a.oncokb, a.input_vcf, a.oncokb_version)
    else:
        ap.error("give --annotated/--oncokb/--input, --sidecar, or --merge")

    if a.out:
        with open(a.out, "w") as fh:
            json.dump(meta, fh, indent=1, sort_keys=True)
        print(f"written: {a.out}")
    if a.summary or not a.out:
        by = {}
        for d in meta["fields"].values():
            by[d["origin"]] = by.get(d["origin"], 0) + 1
        print(f"fields described by the VCFs: {len(meta['fields'])}")
        for k, v in sorted(by.items()):
            print(f"  {k:8s} {v}")
        print("run:")
        for k, v in sorted(meta["run"].items()):
            print(f"  {k:22s} {str(v)[:70]}")
        for level, msg in check(meta):
            print(f"  [{level.upper()}] {msg}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
