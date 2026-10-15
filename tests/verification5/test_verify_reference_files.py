"""
Unit tests for scripts/verify_reference_files.py

`utils/get_ref_files.sh` downloaded and built every reference file the pipeline
annotates from and verified none of them -- no checksum anywhere in it, and one
"is this HTML?" test on one download. The CIViC VCF in the lab reference directory
was structurally corrupt for six months because of it: its eight header lines and
its first data record were one physical line, so the ##INFO declarations were not
parseable as header lines and record 1:8013449 was unreachable by tabix. VEP did not
refuse the file. It annotated, matched nothing, and wrote empty CIViC columns that
look exactly like "this variant is not in CIViC".

Each test below names the defect it exists for. Several are regressions against
false positives found by running the checker against the real reference directory --
those matter as much as the true positives, because a check that fails a sound
reference set is a check that gets switched off.

Covers:
  - the CIViC defect: fused header lines, missing #CHROM line, no parseable ##INFO,
    and the data record the fusion makes unreachable
  - the bootstrapping hazard: a manifest generated from a corrupt directory must not
    bless it
  - truncation, empty files, error-page stubs, plain gzip where bgzip is required,
    stale and mismatched indexes
  - .fai consistency, including the final partial sequence line
  - probes taken from the columns the .tbi was actually built on
  - required versus optional, and the exit statuses entrypoint.sh switches on
"""
import gzip
import importlib.util
import json
import os
import shutil
import struct
import subprocess
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "..", "scripts")
_spec = importlib.util.spec_from_file_location(
    "verify_reference_files", os.path.join(SCRIPTS, "verify_reference_files.py"))
vrf = importlib.util.module_from_spec(_spec)
sys.modules["verify_reference_files"] = vrf
_spec.loader.exec_module(vrf)

HAVE_BGZIP = shutil.which("bgzip") is not None and shutil.which("tabix") is not None
needs_tabix = pytest.mark.skipif(
    not HAVE_BGZIP, reason="bgzip/tabix not on PATH; index behaviour cannot be exercised")


# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------

GOOD_HEADER = [
    "##fileformat=VCFv4.2",
    "##fileDate=2026-04-14",
    "##source=CIViC_nightly_VariantSummaries",
    "##reference=grch38",
    '##INFO=<ID=GN,Number=1,Type=String,Description="CIViC gene name (HGNC symbol)">',
    '##INFO=<ID=VT,Number=1,Type=String,Description="CIViC variant name">',
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="CIViC variant annotation">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
]

RECORDS = [
    ("1", 8013449, "C", "A", "GN=ERRFI1;VT=E384*"),
    ("1", 11106535, "T", "C", "GN=MTOR;VT=Mutation"),
    ("7", 140753336, "A", "T", "GN=BRAF;VT=V600E"),
]


def vcf_text(header=GOOD_HEADER, records=RECORDS, fuse_header=False):
    """A sites-only VCF. `fuse_header` reproduces the CIViC defect exactly.

    The fusion is not "no newlines in the file": the header lines and the first data
    record share one physical line and the remaining records are terminated normally,
    which is what the real artefact looks like (first newline at byte 1454, 1262
    lines in total, one of them matching '^#').
    """
    body = "".join(line + "\n" for line in header)
    rows = ["%s\t%d\t.\t%s\t%s\t.\t.\t%s\n" % r for r in records]
    if fuse_header:
        body = "".join(header) + (rows[0] if rows else "\n")
        rows = rows[1:]
    return body + "".join(rows)


def write_bgzip(path, text, index_args=("-p", "vcf")):
    """bgzip + tabix a text blob. Returns the .gz path."""
    plain = str(path) + ".plain"
    with open(plain, "w") as fh:
        fh.write(text)
    with open(str(path), "wb") as out:
        subprocess.run(["bgzip", "-c", plain], stdout=out, check=True)
    os.remove(plain)
    subprocess.run(["tabix", "-f"] + list(index_args) + [str(path)], check=True)
    return str(path)


def inspect(path, kind, tabix=None, scan=True):
    """Run one inspector and return (facts, defects, warnings)."""
    facts, defects, warnings = {}, [], []
    vrf.INSPECTORS[kind](str(path), tabix or shutil.which("tabix"), scan,
                         defects, facts, warnings)
    return facts, defects, warnings


def write_fasta(tmp_path, name="ref.fa", seq="ACGT" * 15 + "ACG", width=20,
                fai=True, fai_rows=None):
    """A FASTA with a hand-built .fai so the offset arithmetic is exercised.

    The default sequence is 63 bases wrapped at 20, so the last line holds 3 bases --
    the case that made the checker call hg38 truncated.
    """
    fa = tmp_path / name
    header = ">chr1\n"
    lines = [header] + [seq[i:i + width] + "\n" for i in range(0, len(seq), width)]
    fa.write_text("".join(lines))
    if fai:
        rows = fai_rows or [("chr1", len(seq), len(header), width, width + 1)]
        (tmp_path / (name + ".fai")).write_text(
            "".join("%s\t%d\t%d\t%d\t%d\n" % r for r in rows))
    return str(fa)


# ---------------------------------------------------------------------------
# The CIViC defect
# ---------------------------------------------------------------------------

@needs_tabix
def test_fused_header_is_fatal(tmp_path):
    """The defect: eight header lines and record 1 on one physical line.

    This file passes every check `get_ref_files.sh` and `entrypoint.sh` ever made --
    it exists, it is non-empty, it is bgzip-compressed and tabix indexes it without
    complaint. It must not pass this one.
    """
    path = write_bgzip(tmp_path / "civic.vcf.gz", vcf_text(fuse_header=True))
    _facts, defects, _warnings = inspect(path, "vcf")
    assert defects, "a VCF whose header is one fused line was accepted"
    joined = " | ".join(defects)
    assert "single physical line" in joined
    assert "'##' markers" in joined


@needs_tabix
def test_fused_header_loses_the_info_declarations(tmp_path):
    """Why the fusion matters at all: VEP --custom pulls named INFO fields.

    With the ##INFO lines fused into one unparseable line there is nothing for VEP to
    extract, so GN, VT and CSQ come out empty for every variant -- and an empty CIViC
    column reads as "not in CIViC", not as an error.
    """
    path = write_bgzip(tmp_path / "civic.vcf.gz", vcf_text(fuse_header=True))
    facts, defects, _ = inspect(path, "vcf")
    assert facts["info_declarations"] == 0
    assert any("##INFO" in d for d in defects)


@needs_tabix
def test_fused_header_names_the_unreachable_record(tmp_path):
    """The record the fusion swallows is reported by position, not just 'odd file'.

    On the real artefact this is 1:8013449 (ERRFI1 E384*). tabix sees a line starting
    with '#', skips it as a comment, and returns nothing -- so a report on that
    variant says "not found" instead of failing.
    """
    path = write_bgzip(tmp_path / "civic.vcf.gz", vcf_text(fuse_header=True))
    facts, defects, _ = inspect(path, "vcf")
    assert facts.get("record_hidden_in_header") == ["1", 8013449]
    assert any("8013449" in d and "cannot be retrieved" in d for d in defects)


@needs_tabix
def test_a_sound_vcf_passes(tmp_path):
    """The same content, correctly framed, must produce no findings at all.

    Paired with the test above deliberately: the two files differ only in eight
    newline characters, so this is what pins the checker to the defect rather than to
    the data.
    """
    path = write_bgzip(tmp_path / "civic.vcf.gz", vcf_text())
    facts, defects, _ = inspect(path, "vcf")
    assert defects == []
    assert facts["header_lines"] == 8
    assert facts["records"] == 3
    assert facts["info_declarations"] == 3
    assert facts["probe_hits"] >= 1


# ---------------------------------------------------------------------------
# The bootstrapping hazard
# ---------------------------------------------------------------------------

@needs_tabix
def test_generate_refuses_to_certify_a_corrupt_directory(tmp_path, monkeypatch):
    """A manifest generated from a damaged directory would record damage as correct.

    The only reference directory available to generate the first manifest from is the
    one containing the corrupt CIViC file, so this is the failure mode that would
    otherwise defeat the whole exercise. `generate` runs the same checks as `verify`
    and exits non-zero rather than quietly recording the corruption as expected.
    """
    refs = tmp_path / "refs"
    (refs / "CIVIC").mkdir(parents=True)
    write_bgzip(refs / "CIVIC" / "civic_grch38.vcf.gz", vcf_text(fuse_header=True))
    manifest = tmp_path / "m.json"

    rc = vrf.main(["generate", "--refs", str(refs), "--manifest", str(manifest),
                   "--only", "civic"])
    assert rc == vrf.EXIT_FATAL, "generate certified a corrupt reference directory"

    recorded = json.loads(manifest.read_text())
    civic = recorded["files"][0]
    assert civic["defects"], "the corruption was recorded as the expected state"


@needs_tabix
def test_verify_against_a_manifest_generated_from_corruption_still_fails(tmp_path):
    """Even if such a manifest is committed and shipped, it cannot bless the file.

    Every byte matches and every recorded fact matches -- that is the point. A
    digest-only design passes here, which is why the manifest carries the structural
    verdict as well as the hash.
    """
    refs = tmp_path / "refs"
    (refs / "CIVIC").mkdir(parents=True)
    write_bgzip(refs / "CIVIC" / "civic_grch38.vcf.gz", vcf_text(fuse_header=True))
    manifest = tmp_path / "m.json"
    vrf.main(["generate", "--refs", str(refs), "--manifest", str(manifest), "--only", "civic"])

    rc = vrf.main(["verify", "--refs", str(refs), "--manifest", str(manifest), "--only", "civic"])
    assert rc == vrf.EXIT_FATAL


@needs_tabix
def test_manifest_defect_is_not_reported_twice(tmp_path, capsys):
    """A still-broken file is reported once, not once live and once from the manifest.

    Duplicated findings make a report look like more defects than there are, and a
    clinical log that cries wolf stops being read.
    """
    refs = tmp_path / "refs"
    (refs / "CIVIC").mkdir(parents=True)
    write_bgzip(refs / "CIVIC" / "civic_grch38.vcf.gz", vcf_text(fuse_header=True))
    manifest = tmp_path / "m.json"
    vrf.main(["generate", "--refs", str(refs), "--manifest", str(manifest), "--only", "civic"])
    capsys.readouterr()

    vrf.main(["verify", "--refs", str(refs), "--manifest", str(manifest), "--only", "civic"])
    out = capsys.readouterr().out
    assert out.count("single physical line") == 1


@needs_tabix
def test_manifest_records_the_bootstrapping_warning(tmp_path):
    """The hazard is written into the manifest, not only into the source.

    A manifest outlives the knowledge of how it was made, and the first question
    asked of a failing verify is whether the manifest or the data is wrong.
    """
    refs = tmp_path / "refs"
    (refs / "CIVIC").mkdir(parents=True)
    write_bgzip(refs / "CIVIC" / "civic_grch38.vcf.gz", vcf_text())
    manifest = tmp_path / "m.json"
    vrf.main(["generate", "--refs", str(refs), "--manifest", str(manifest), "--only", "civic"])
    recorded = json.loads(manifest.read_text())
    assert "damaged" in recorded["bootstrapping_warning"]
    assert recorded["generated_from"] == os.path.abspath(str(refs))


# ---------------------------------------------------------------------------
# The other failures that actually happen to reference files
# ---------------------------------------------------------------------------

def test_empty_file_is_fatal(tmp_path):
    """A zero-byte file passes `[[ -e ]]`, which is the only test entrypoint.sh had."""
    path = tmp_path / "clinvar.vcf.gz"
    path.write_bytes(b"")
    rec = vrf.inspect_entry(str(tmp_path), "x", "clinvar.vcf.gz", "vcf", True, "",
                            None, 1 << 30, 1 << 28, True)
    assert any("0 bytes" in d for d in rec["defects"])


def test_html_error_page_saved_as_data_is_fatal(tmp_path):
    """A quota stub or error page served with HTTP 200 and written to the data path.

    `get_ref_files.sh` already guarded against this on one Google Drive download;
    the guard belongs on every file, not on the one where it was first noticed.
    """
    path = tmp_path / "revel.tsv.gz"
    path.write_bytes(b"<!doctype html><html><body>Quota exceeded</body></html>")
    defects = []
    vrf._stub_check(str(path), defects)
    assert any("error page or quota stub" in d for d in defects)


def test_truncated_gzip_is_reported_as_damage_not_as_a_crash(tmp_path):
    """A half-finished download: gzip raises mid-stream.

    Reported by name rather than allowed to surface as an unhandled traceback, which
    reads as a bug in the checker rather than damage in the reference file.
    """
    path = tmp_path / "t.tsv.gz"
    full = gzip.compress(b"#chrom\tpos\tscore\n" + b"".join(
        b"1\t%d\t0.5\n" % i for i in range(1, 5000)))
    path.write_bytes(full[:len(full) // 2])
    defects = []
    vrf.inspect_header_and_records(str(path), "tsv", True, defects)
    assert any("truncated or corrupt" in d for d in defects)


def test_plain_gzip_where_bgzip_is_required_is_fatal(tmp_path):
    """A table rebuilt with `gzip` instead of `bgzip`.

    Readable by everything that decompresses it and unqueryable by tabix, which is
    the only tool that matters: annotation then matches nothing and reports no error.
    """
    path = tmp_path / "x.vcf.gz"
    path.write_bytes(gzip.compress(vcf_text().encode()))
    assert vrf._is_bgzf(str(path)) is False
    _facts, defects, _ = inspect(path, "vcf")
    assert any("not BGZF" in d for d in defects)


@needs_tabix
def test_bgzip_is_recognised_as_bgzf(tmp_path):
    """Regression: the BGZF test must not reject real bgzip output.

    It reads the 'BC' extra subfield out of the gzip header, and a false negative here
    would fail every reference table in the set.
    """
    path = write_bgzip(tmp_path / "x.vcf.gz", vcf_text())
    assert vrf._is_bgzf(path) is True


@needs_tabix
def test_index_from_a_previous_version_of_the_file_is_detected(tmp_path):
    """Data re-downloaded, .tbi kept from the old file.

    The index's offsets then point into bytes that have moved, so queries return the
    wrong records or none, and never an error. The probe is what proves it: it asks
    for a position the file itself declares.
    """
    path = write_bgzip(tmp_path / "x.vcf.gz", vcf_text())
    idx = path + ".tbi"
    # Replace the data with a differently sized file, keeping the old index.
    other = str(tmp_path / "other.vcf.gz")
    write_bgzip(tmp_path / "other.vcf.gz", vcf_text(records=[("5", 500, "A", "G", "GN=X")]))
    shutil.copyfile(other, path)
    shutil.copystat(other, path)
    os.utime(idx, (os.path.getmtime(path) - 3600,) * 2)

    _facts, defects, warnings = inspect(path, "vcf")
    assert defects or warnings, "a stale index was accepted"
    assert any("older than the data" in w for w in warnings)


@needs_tabix
def test_missing_index_is_fatal(tmp_path):
    """Nothing in the pipeline checks that a .tbi exists beside its data."""
    path = write_bgzip(tmp_path / "x.vcf.gz", vcf_text())
    os.remove(path + ".tbi")
    _facts, defects, _ = inspect(path, "vcf")
    assert any("index" in d and "missing" in d for d in defects)


@needs_tabix
def test_index_mtime_older_than_data_is_a_warning_not_a_failure(tmp_path):
    """Regression: CADD ships its own .tbi, downloaded before its 81 GB data file.

    Both CADD indexes in the real reference directory are legitimately older than
    their data and both probe correctly. Making the mtime comparison fatal failed a
    sound reference set, so it is a warning and the probe decides.
    """
    path = write_bgzip(tmp_path / "x.vcf.gz", vcf_text())
    os.utime(path + ".tbi", (os.path.getmtime(path) - 86400,) * 2)
    _facts, defects, warnings = inspect(path, "vcf")
    assert defects == []
    assert any("older than the data" in w for w in warnings)


@needs_tabix
def test_header_fused_onto_a_tsv_column_line_is_detected(tmp_path):
    """Two tables concatenated: the second file's header lands mid-stream."""
    text = ("#chrom\tpos\tscore\n"
            "1\t100\t0.1\n"
            "#chrom\tpos\tscore\n"
            "1\t200\t0.2\n")
    path = write_bgzip(tmp_path / "t.tsv.gz", text, index_args=("-s", "1", "-b", "2", "-e", "2"))
    _facts, defects, _ = inspect(path, "tsv")
    assert any("concatenated" in d for d in defects)


@needs_tabix
def test_column_count_shift_is_detected(tmp_path):
    """A header declaring a different number of columns from its data.

    Every value annotated out of the file is then read from the wrong column, which
    produces plausible-looking wrong numbers rather than an error.
    """
    text = "#chrom\tpos\tscore\textra\n1\t100\t0.1\n1\t200\t0.2\n"
    path = write_bgzip(tmp_path / "t.tsv.gz", text, index_args=("-s", "1", "-b", "2", "-e", "2"))
    _facts, defects, _ = inspect(path, "tsv")
    assert any("column mapping" in d for d in defects)


# ---------------------------------------------------------------------------
# Probes taken from the index's own column layout
# ---------------------------------------------------------------------------

@needs_tabix
def test_probe_uses_the_columns_the_index_was_built_on(tmp_path):
    """Regression: LOEUF is indexed on columns 56/57/58, REVEL's position is column 3.

    Assuming columns 1 and 2 took LOEUF's first field -- a gene symbol -- as a contig
    and queried `OR4F5:79501`, which returns nothing. That is indistinguishable from
    the real CIViC defect, and a check that cries wolf on a sound file is worse than
    no check.
    """
    text = ("#gene\tgene_id\tchrom\tstart\tend\n"
            "OR4F5\tENSG00000186092\tchr1\t65419\t71585\n"
            "SAMD11\tENSG00000187634\tchr1\t923928\t944581\n")
    path = write_bgzip(tmp_path / "loeuf.tsv.gz", text,
                       index_args=("-s", "3", "-b", "4", "-e", "5", "-S", "1"))
    cfg = vrf._read_tbi_config(path + ".tbi")
    assert (cfg["col_seq"], cfg["col_beg"]) == (3, 4)

    facts, defects, _ = inspect(path, "tsv")
    assert defects == []
    assert facts["probe"] == ["chr1", 65419]


@needs_tabix
def test_interval_indexed_table_probe_accepts_a_covering_row(tmp_path):
    """Regression: a gene-interval table's returned row need not start at the probe.

    Comparing the returned record's second column against the queried position
    reported "the index does not agree with the data" for LOEUF, whose second column
    is a gene id. Coverage of the interval is the right test.
    """
    text = ("#gene\tgene_id\tchrom\tstart\tend\n"
            "OR4F5\tENSG00000186092\tchr1\t65419\t71585\n")
    path = write_bgzip(tmp_path / "loeuf.tsv.gz", text,
                       index_args=("-s", "3", "-b", "4", "-e", "5", "-S", "1"))
    facts, defects, _ = inspect(path, "tsv")
    assert defects == []
    assert facts["probe_covering_rows"] >= 1


@needs_tabix
def test_vcf_indexed_as_a_generic_table_is_fatal(tmp_path):
    """`tabix -s 1 -b 2` on a VCF indexes the wrong coordinate convention.

    Queries near indels then return the wrong window, silently.
    """
    path = write_bgzip(tmp_path / "x.vcf.gz", vcf_text(),
                       index_args=("-s", "1", "-b", "2", "-e", "2"))
    _facts, defects, _ = inspect(path, "vcf")
    assert any("format 2" in d for d in defects)


# ---------------------------------------------------------------------------
# FASTA / .fai
# ---------------------------------------------------------------------------

def test_fai_arithmetic_accepts_a_final_partial_line(tmp_path):
    """Regression: rounding the last contig up to a whole line over-states the file.

    On hg38 that charged a full line's newline to a 40-base final line and declared
    the genome truncated by 10 bytes -- a false positive on the one file every run
    depends on, which is the fastest way to get a check switched off.
    """
    assert vrf._fasta_contig_bytes(37240, 50, 51) == 37985
    assert vrf._fasta_contig_bytes(100, 50, 51) == 102     # exact multiple, no partial line
    fa = write_fasta(tmp_path)
    _facts, defects, _ = inspect(fa, "fasta")
    assert defects == []


def test_missing_fai_is_fatal(tmp_path):
    """VEP and scripts/check_reference.py both read the genome through its .fai."""
    fa = write_fasta(tmp_path, fai=False)
    _facts, defects, _ = inspect(fa, "fasta")
    assert any(".fai index is missing" in d for d in defects)


def test_fai_describing_more_bytes_than_the_fasta_has_is_fatal(tmp_path):
    """A truncated FASTA, or a .fai belonging to a different build."""
    fa = write_fasta(tmp_path)
    data = open(fa, "rb").read()
    open(fa, "wb").write(data[:len(data) // 2])
    _facts, defects, _ = inspect(fa, "fasta")
    assert any("truncated or the index belongs to a different file" in d for d in defects)


def test_fai_offsets_pointing_at_the_wrong_contig_are_fatal(tmp_path):
    """A consistent-looking .fai from another assembly passes the arithmetic.

    scripts/check_reference.py seeks to these offsets to read reference bases, so a
    .fai from another build makes its REF-allele check report mismatches for a sound
    VCF -- or agreement for an unsound one. The offset must land on sequence preceded
    by that contig's own '>' header.
    """
    fa = write_fasta(tmp_path, fai_rows=[("chr9", 63, 6, 20, 21)])
    _facts, defects, _ = inspect(fa, "fasta")
    assert any("does not belong to this FASTA" in d for d in defects)


def test_stale_fai_is_detected(tmp_path):
    """FASTA rebuilt, .fai left from the previous file."""
    fa = write_fasta(tmp_path)
    os.utime(fa, (os.path.getmtime(fa + ".fai") + 3600,) * 2)
    _facts, defects, _ = inspect(fa, "fasta")
    assert any(".fai is older than the FASTA" in d for d in defects)


# ---------------------------------------------------------------------------
# Other kinds
# ---------------------------------------------------------------------------

def test_chain_file_without_its_chain_header_is_fatal(tmp_path):
    """A liftover chain that lost its first line still lifts over, silently and wrongly."""
    path = tmp_path / "hg19ToHg38.over.chain"
    path.write_text("1000 chr1 248956422 + 10000 10100 chr1 249250621 + 10000 10100 1\n")
    _facts, defects, _ = inspect(path, "text")
    assert any("does not start with 'chain'" in d for d in defects)


def test_empty_resource_directory_is_fatal(tmp_path):
    """A failed extraction leaves the directory behind.

    `entrypoint.sh` tests `[[ -e "$VEP_PLUGIN_DIR" ]]`, which an empty directory
    passes, and VEP then runs with no plugins and annotates fewer columns.
    """
    d = tmp_path / "Plugins"
    d.mkdir()
    _facts, defects, _ = inspect(d, "dir")
    assert any("directory is empty" in x for x in defects)


def test_gztext_counts_rows_without_needing_an_index(tmp_path):
    """Regression: GENIE counts and the SVIG-UK tables are read whole, not queried.

    `scripts/genie_annotation.py` opens genie_o4_counts.tsv.gz directly, so it is
    plain gzip with no .tbi. Checking it as an indexed table produced three defects
    against a sound file: not BGZF, no CHROM/POS in column 1/2, and a missing index.
    """
    path = tmp_path / "genie_o4_counts.tsv.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("gene\tvariant_class\tchrom\tstart\toncotree_code\tn_patients\n")
        fh.write("BRAF\tmissense\t7\t140753336\tSKCM\t42\n")
    facts, defects, _ = inspect(path, "gztext")
    assert defects == []
    # 2, not 1: these tables carry no comment marker on their header row, so it is
    # counted among the rows. Deliberate -- the count is only ever compared against
    # itself, and being consistently off by one still detects truncation.
    assert facts["records"] == 2

    entry = [e for e in vrf.REGISTRY if e[0] == "genie_o4_counts"][0]
    assert entry[2] == "gztext", "GENIE must not be checked as a tabix-indexed table"


# ---------------------------------------------------------------------------
# Required versus optional, and the exit statuses entrypoint.sh switches on
# ---------------------------------------------------------------------------

def test_missing_required_file_is_fatal(tmp_path):
    """entrypoint.sh must not annotate with ClinVar absent."""
    rc = vrf.main(["verify", "--refs", str(tmp_path), "--only", "clinvar"])
    assert rc == vrf.EXIT_FATAL


def test_absent_optional_file_is_not_a_failure(tmp_path):
    """dbNSFP, CADD, COSMIC and GENIE are optional per entrypoint.sh.

    It enables each annotation on an `-e` test and runs unchanged without it, so
    absence is a configuration choice. Treating it as a defect would fail every
    deployment that has not licensed COSMIC.
    """
    rc = vrf.main(["verify", "--refs", str(tmp_path),
                   "--only", "cosmic", "cadd_snv", "dbnsfp", "genie_o4_counts"])
    assert rc == vrf.EXIT_OK


@needs_tabix
def test_present_but_damaged_optional_file_is_fatal(tmp_path):
    """Absence is a choice; damage is not.

    entrypoint.sh adds `--custom $COSMIC_VCF,...` on an existence test alone, so a
    damaged COSMIC VCF is used rather than skipped, and the columns it feeds come out
    empty while looking annotated.
    """
    refs = tmp_path / "refs"
    (refs / "COSMIC").mkdir(parents=True)
    write_bgzip(refs / "COSMIC" / "cosmic_cmc_grch38.vcf.gz", vcf_text(fuse_header=True))
    rc = vrf.main(["verify", "--refs", str(refs), "--only", "cosmic"])
    assert rc == vrf.EXIT_FATAL


@needs_tabix
def test_updated_but_sound_file_warns_and_fails_only_under_strict(tmp_path):
    """A changed digest with intact structure is a reference update, not damage.

    ClinVar is released weekly and CIViC nightly. Making every digest change fatal
    means a pipeline that refuses to run whenever someone updates ClinVar, and a
    check that refuses to run is a check that gets commented out. --strict is for a
    frozen validated reference set, where any change at all is the thing you want.
    """
    refs = tmp_path / "refs"
    (refs / "CIVIC").mkdir(parents=True)
    target = refs / "CIVIC" / "civic_grch38.vcf.gz"
    write_bgzip(target, vcf_text())
    manifest = tmp_path / "m.json"
    assert vrf.main(["generate", "--refs", str(refs), "--manifest", str(manifest),
                     "--only", "civic"]) == vrf.EXIT_OK

    # A legitimate nightly update: one more record, same structure.
    write_bgzip(target, vcf_text(records=RECORDS + [("12", 25245350, "C", "T", "GN=KRAS")]))

    assert vrf.main(["verify", "--refs", str(refs), "--manifest", str(manifest),
                     "--only", "civic"]) == vrf.EXIT_OK
    assert vrf.main(["verify", "--refs", str(refs), "--manifest", str(manifest),
                     "--only", "civic", "--strict"]) == vrf.EXIT_WARNINGS_STRICT


def test_missing_reference_directory_is_distinguished_from_a_bad_reference(tmp_path):
    """Exit 2 means the check could not run; exit 1 means the reference is bad.

    entrypoint.sh switches on the two separately so that a broken container is never
    reported as a broken reference, or the reverse.
    """
    rc = vrf.main(["verify", "--refs", str(tmp_path / "nope")])
    assert rc == vrf.EXIT_CANNOT_RUN


def test_named_but_missing_manifest_cannot_run(tmp_path):
    """Asking for a manifest that is not there is a usage error, not a pass."""
    rc = vrf.main(["verify", "--refs", str(tmp_path), "--manifest", str(tmp_path / "no.json")])
    assert rc == vrf.EXIT_CANNOT_RUN


def test_verify_works_with_no_manifest_at_all(tmp_path):
    """Structural checks alone still catch the CIViC class of defect.

    This is what entrypoint.sh falls back to on a deployment that predates the
    manifest, so that adding the check does not require every site to regenerate one
    before it can run -- a requirement people route around rather than meet.
    """
    rc = vrf.main(["verify", "--refs", str(tmp_path), "--only", "cosmic"])
    assert rc == vrf.EXIT_OK


# ---------------------------------------------------------------------------
# Digests
# ---------------------------------------------------------------------------

def test_large_files_get_a_partial_digest_that_is_labelled_as_one(tmp_path):
    """CADD SNVs are 81 GB; a full SHA-256 of the set is not affordable at startup.

    A check that takes an hour gets skipped, so files above the threshold get a
    head+tail+size digest. The kind is recorded so a partial digest can never be
    compared against a full one and reported as a mismatch on a sound file.
    """
    path = tmp_path / "big.bin"
    path.write_bytes(b"x" * 4096)
    full_kind, full = vrf._digest(str(path), 4096, max_hash_bytes=1 << 20)
    part_kind, part = vrf._digest(str(path), 4096, max_hash_bytes=1024)
    assert full_kind == "sha256"
    assert part_kind.startswith("sha256-head")
    assert full != part


def test_partial_digest_detects_truncation(tmp_path):
    """The head+tail+size digest is chosen for the failure that actually occurs.

    A dropped connection or a builder that died mid-write always changes the size and
    the tail, which is what the partial digest covers -- unlike a head-only sample.
    """
    path = tmp_path / "big.bin"
    path.write_bytes(b"a" * 4096 + b"b" * 4096)
    _k1, before = vrf._digest(str(path), 8192, max_hash_bytes=1024)
    path.write_bytes(b"a" * 4096)
    _k2, after = vrf._digest(str(path), 4096, max_hash_bytes=1024)
    assert before != after


def test_tbi_config_reader_rejects_a_non_tbi_file(tmp_path):
    """A .tbi that is not a .tbi must not be read as column numbers.

    Returning None makes the probe unavailable and reported as unchecked, rather than
    silently probing column 0 of nothing.
    """
    path = tmp_path / "x.tbi"
    path.write_bytes(gzip.compress(b"not an index at all" + b"\0" * 64))
    assert vrf._read_tbi_config(str(path)) is None


def test_registry_covers_every_required_resource_entrypoint_checks(tmp_path):
    """The registry must not drift from what entrypoint.sh actually annotates with.

    A reference file that entrypoint.sh requires but the registry omits is a file
    nothing verifies, which is the hole this whole module exists to close.
    """
    entrypoint = os.path.join(os.path.dirname(__file__), "..", "..", "entrypoint.sh")
    text = open(entrypoint).read()
    required = {name: rel for name, rel, _k, req, _n in vrf.REGISTRY if req}
    for rel in required.values():
        if rel in ("homo_sapiens", "Plugins"):
            continue        # referenced as $VEP_DIR / $VEP_PLUGIN_DIR, not by path
        assert rel in text, "%s is in the registry but entrypoint.sh does not use it" % rel


class TestTheContigListDescribesChromosomes:
    """The contig list is the fact that says which chromosomes a reference covers.

    It was derived from column 1 for every kind of file. That is CHROM in a VCF, but a
    tabix-indexed TSV puts the sequence name wherever its builder chose: LOEUF's is
    column 56, and column 1 is a gene symbol. So LOEUF's manifest entry recorded 18,143
    gene names as its contigs -- a 115 KiB fact that said nothing about chromosomes, and
    so could not show that gnomAD v4.1 omits chrX and chrY, which is the reason O8 is
    unevaluable for every X-linked tumour suppressor.
    """

    def tsv_with_seq_in_last_column(self, rows):
        """A LOEUF-shaped table: gene symbol first, coordinates appended last."""
        head = "#gene\ttranscript\tscore\tchrom\tstart\tend"
        body = "\n".join("%s\t%s\t0.5\t%s\t100\t900" % r for r in rows)
        return head + "\n" + body + "\n"

    @needs_tabix
    def test_contigs_come_from_the_index_not_column_one(self, tmp_path):
        rows = [("GENE%d" % i, "ENST%d" % i, c)
                for i, c in enumerate(["chr1", "chr1", "chr2", "chr17"])]
        path = write_bgzip(tmp_path / "loeuf.tsv.gz",
                           self.tsv_with_seq_in_last_column(rows),
                           index_args=("-s", "4", "-b", "5", "-e", "6", "-S", "1"))
        facts, defects = vrf.inspect_tsv(path, "tabix", True, [], {}, [])
        assert facts["index_columns"] == [4, 5, 6]
        assert facts["contigs"] == ["chr1", "chr2", "chr17"], facts["contigs"]
        # the gene symbols must not appear -- that was the defect
        assert not any(c.startswith("GENE") for c in facts["contigs"])
        assert defects == []

    @needs_tabix
    def test_a_sex_chromosome_would_now_be_visible(self, tmp_path):
        """The point of fixing it: whether a constraint table covers chrX is a fact the
        manifest must carry, because O8's behaviour for X-linked genes depends on it."""
        rows = [("A", "T1", "chr1"), ("B", "T2", "chrX")]
        path = write_bgzip(tmp_path / "loeuf.tsv.gz",
                           self.tsv_with_seq_in_last_column(rows),
                           index_args=("-s", "4", "-b", "5", "-e", "6", "-S", "1"))
        facts, _ = vrf.inspect_tsv(path, "tabix", True, [], {}, [])
        assert "chrX" in facts["contigs"]

    def test_a_vcf_still_reads_its_contigs_from_column_one(self):
        """The default must not change: CHROM is a VCF's first field."""
        defects = []
        import inspect as _inspect
        sig = _inspect.signature(vrf.inspect_header_and_records)
        assert sig.parameters["seq_col"].default == 0

    def test_an_implausible_count_is_reported_rather_than_recorded(self, tmp_path):
        """If the sequence column is ever misidentified again, the symptom is thousands
        of distinct values. Recording them hides the mistake; a defect surfaces it."""
        rows = [("G%d" % i, "T%d" % i, "chr1") for i in range(vrf.MAX_PLAUSIBLE_CONTIGS + 5)]
        plain = tmp_path / "t.tsv"
        plain.write_text(self.tsv_with_seq_in_last_column(rows))
        import gzip as _gz
        gz = tmp_path / "t.tsv.gz"
        with _gz.open(gz, "wt") as fh:
            fh.write(plain.read_text())
        defects = []
        # column 1 is the gene symbol: the mistake being guarded against
        _, _, _, n, contigs = vrf.inspect_header_and_records(
            str(gz), "tsv", True, defects, seq_col=0)
        assert contigs == []
        assert any("too many to be sequence names" in d for d in defects), defects
        assert n == vrf.MAX_PLAUSIBLE_CONTIGS + 5

    def test_the_plausible_ceiling_allows_a_real_reference(self):
        """GRCh38 with unplaced contigs and alts is a few hundred sequences."""
        assert vrf.MAX_PLAUSIBLE_CONTIGS >= 500


class TestTheCacheItRequiresIsTheCacheVepUses:
    """It required `homo_sapiens`; the pipeline runs `vep --merged`, which needs
    `homo_sapiens_merged`.

    Both directions were wrong. A reference directory holding only the merged cache --
    the correct one -- failed verification. One holding only the plain cache was
    certified and then failed at VEP, which is the worse direction, because the whole
    point of this script is to refuse a reference set before a clinical run annotates
    against it. A machine with both caches sees neither failure, which is how it lasted.
    """

    ROOT = os.path.join(os.path.dirname(__file__), "..", "..")

    def test_the_pipeline_really_runs_merged(self):
        """The premise. If `--merged` is ever dropped, this expectation inverts."""
        entrypoint = open(os.path.join(self.ROOT, "entrypoint.sh"), encoding="utf-8").read()
        assert "--merged" in entrypoint

    def test_the_required_cache_is_the_merged_one(self):
        cache = [r for r in vrf.REGISTRY if r[0] == "vep_cache"]
        assert cache, "vep_cache is no longer a declared resource"
        assert cache[0][1] == "homo_sapiens_merged", (
            f"it requires {cache[0][1]!r}, but `vep --merged` looks for "
            f"homo_sapiens_merged")

    def test_it_is_still_required_rather_than_optional(self):
        """Annotation cannot proceed without it, so a missing cache must fail rather
        than warn."""
        cache = [r for r in vrf.REGISTRY if r[0] == "vep_cache"][0]
        assert cache[3] is True

    def test_the_plain_cache_is_not_required_as_well(self):
        """The pipeline never reads it, so requiring it would reject a correctly
        provisioned directory for holding only what is needed."""
        paths = [r[1] for r in vrf.REGISTRY]
        assert "homo_sapiens" not in paths
