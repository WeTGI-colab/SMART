#!/usr/bin/env python3
"""
Verify SMART's reference files against a manifest before a run annotates from them.

Why this exists
---------------
`utils/get_ref_files.sh` downloads and builds every reference file the pipeline
annotates from and verifies none of them. There is no checksum anywhere in it
(`grep -c 'sha256\\|md5\\|shasum' utils/get_ref_files.sh` returns 0) and its only
sanity guard is an "is this HTML rather than data?" test on one Google Drive
download. Worse, every step is guarded by `[[ -f "$FILE" ]] && echo "[SKIP]"`, so a
file that exists is never looked at again, whatever is inside it.

That is not a theoretical gap. The CIViC VCF in the reference directory on the
machine this was written on is structurally corrupt and has been annotated against
for six months:

    $ gzip -dc refs/CIVIC/civic_grch38.vcf.gz | wc -l         -> 1262
    $ gzip -dc refs/CIVIC/civic_grch38.vcf.gz | grep -c '^#'  -> 1
    $ tabix refs/CIVIC/civic_grch38.vcf.gz 1:8013449-8013449  -> (empty)

Its eight header lines and its first data record are a single physical line with no
newlines between them. The consequences are exactly the kind that produce output
rather than an error: the three `##INFO` declarations are not parseable as header
lines, and the first record (ERRFI1 E384*) is unreachable because tabix sees a line
beginning with `#` and skips it as a comment. VEP does not refuse this file. It
annotates, finds nothing, and writes empty CIViC columns -- which look exactly like
"this variant is not in CIViC".

The failures that actually happen to reference files are:

    truncated download      (connection dropped; size and tail change)
    quota / error stub      (HTML or JSON where data was expected)
    partially written build (builder died mid-write)
    mis-concatenated build  (what happened to CIViC: content correct, framing wrong)
    stale index             (data re-downloaded, .tbi left from the old file)
    plain gzip not bgzip    (readable by gzip, unqueryable by tabix)

A checksum alone catches the first three and misses the last three, because a
mis-framed or stale-indexed file is a perfectly consistent sequence of bytes. So
this script records *structural* facts beside the digest -- how many header lines,
how many records, which contigs, and whether a query for a position the file itself
says it contains actually returns that record.

The bootstrapping hazard
------------------------
A manifest is only ever as good as the directory it was generated from. Generating
one from a damaged directory records the damage as the expected state, and every
later run then "passes" against corruption. This is not hypothetical here either:
the CIViC file in the only reference directory available to generate from is the
corrupt one.

So `generate` is not a recording tool. It runs the same structural checks as
`verify` and, for any file that fails them, writes the failures into the manifest
entry as `defects` and exits non-zero. `verify` then treats an entry carrying
`defects` as a failure of that file regardless of whether the bytes still match, so
a manifest generated from the corrupt directory cannot bless it even if it is
committed and shipped. The only way to get a clean manifest is to fix the file.

Usage
-----
    # record the expected state (refuses to bless a file that fails its checks)
    python3 scripts/verify_reference_files.py generate \\
        --refs /path/to/refs --manifest /path/to/refs/smart_reference_manifest.json

    # check a directory before annotating
    python3 scripts/verify_reference_files.py verify \\
        --refs /path/to/refs --manifest /path/to/refs/smart_reference_manifest.json

    # structure only -- no manifest needed, still catches the CIViC class of defect
    python3 scripts/verify_reference_files.py verify --refs /path/to/refs

Exit codes (what `entrypoint.sh` acts on)
-----------------------------------------
    0   every required file present and sound
    1   at least one FATAL finding -- do not annotate
    2   the check itself could not run (bad arguments, no reference directory,
        no tabix on PATH). Distinguished from 1 so that a broken checker is never
        mistaken for a broken reference, and vice versa.
    3   warnings only, and --strict was given

Fatal versus warning
--------------------
    required file missing                       FATAL
    required or optional file structurally bad  FATAL
    optional file absent                        INFO -- `entrypoint.sh` already
                                                treats dbNSFP, CADD, COSMIC and
                                                GENIE as optional and runs without
                                                them, so absence is a configuration
                                                choice, not a defect
    bytes changed, structure still sound        WARNING (FATAL under --strict)

That last rule is deliberate and is the one most likely to be argued with. ClinVar
is released weekly and CIViC nightly, so a reference directory's digests legitimately
change. Making every such change fatal would mean a pipeline that refuses to run
every time someone updates ClinVar, and a check that refuses to run is a check that
gets commented out. A changed digest with intact structure means "this is not the
file that was validated" -- which belongs in the run log, loudly, next to the
provenance block -- while a changed digest with broken structure means "this file is
damaged", which stops the run. Use --strict for a frozen validated reference set,
where any change at all is the thing you want to hear about.

A note on size
--------------
The reference set is ~190 GB: CADD SNVs alone are 81 GB and SpliceAI indels 65 GB.
A full SHA-256 of all of it is not something a clinical run can afford at startup,
and a check that takes an hour will be skipped. Files above --max-hash-bytes get a
`head+tail+size` digest instead: the first and last 8 MiB plus the exact byte count.
That is chosen to match the failure modes above -- a truncated or partially written
file always changes the size and the tail -- and the manifest records which kind of
digest each file carries so a partial digest is never mistaken for a full one.
Likewise, record counting requires decompressing the whole file, so it is done only
below --max-scan-bytes; larger files rely on the index's own contig list and on the
query probe, both of which are O(1).
"""
import argparse
import gzip
import hashlib
import io
import json
import math
import os
import shutil
import re
import struct
import subprocess
import sys
import time
import zlib

MANIFEST_VERSION = 2

# Default digest policy. 1 GiB covers every reference file except the five large
# score tables (CADD x2, SpliceAI x2, dbNSFP); those get the head+tail+size digest.
DEFAULT_MAX_HASH_BYTES = 1024 * 1024 * 1024
DEFAULT_MAX_SCAN_BYTES = 256 * 1024 * 1024     # compressed size above which we do not count records
HASH_EDGE_BYTES = 8 * 1024 * 1024              # bytes taken from each end for a partial digest
HASH_BLOCK = 4 * 1024 * 1024

# A VCF header is small. Reading past this means we never found a data record, which
# is itself the finding -- so bound the read rather than streaming a 65 GB file in
# search of a line that does not start with '#'.
MAX_HEADER_BYTES = 64 * 1024 * 1024

# A reference genome has a few hundred sequences at most, counting unplaced contigs and
# alts. Anything beyond this is not a sequence column, and recording it as one hides the
# fact the list exists to carry: which chromosomes the file actually covers.
MAX_PLAUSIBLE_CONTIGS = 1000

SEVERITY_FATAL = "FATAL"
SEVERITY_WARN = "WARNING"
SEVERITY_INFO = "INFO"

EXIT_OK = 0
EXIT_FATAL = 1
EXIT_CANNOT_RUN = 2
EXIT_WARNINGS_STRICT = 3

# Markers of an error page or quota stub served with a 200 and saved as data. This is
# the one guard `get_ref_files.sh` already had (on a single Google Drive download);
# it is kept and applied to every file rather than to one.
STUB_MARKERS = (
    b"<!doctype html", b"<html", b"<?xml", b"<HTML",
    b"Quota exceeded", b"quota exceeded",
    b"Google Drive - Virus scan warning",
    b"AccessDenied", b"<Error>", b"404 Not Found", b"Not Found",
)


# ---------------------------------------------------------------------------
# The registry: which files the pipeline annotates from, and what each one is.
#
# Paths are relative to the reference root and are taken from `entrypoint.sh`'s
# own definitions, not from `pipeline_config.sh` -- the two disagree about the
# genome FASTA (`liftover/hg38.fa` versus `reference/hg38.fa`) and `entrypoint.sh`
# is the one that annotates, so it is the one that decides what "present" means.
# ---------------------------------------------------------------------------
REGISTRY = [
    # name                 relative path                                               kind          required  note
    ("reference_fasta",    "liftover/hg38.fa",                                          "fasta",      True,  "GRCh38 genome; VEP and scripts/check_reference.py read it through its .fai"),
    ("liftover_chain",     "liftover/hg19ToHg38.over.chain",                            "text",       True,  "hg19->hg38 chain used by the CIViC and COSMIC builders"),
    ("clinvar",            "ClinVar/clinvar.vcf.gz",                                    "vcf",        True,  "VEP --custom; 44 output columns"),
    ("civic",              "CIVIC/civic_grch38.vcf.gz",                                 "vcf",        True,  "VEP --custom GN,VT,CSQ; built by utils/civic_formating.py"),
    ("cancer_hotspots",    "CancerHotSpots/hg38.hotspots_changv2_gao_nc.vcf.gz",        "vcf",        True,  "VEP --custom; 10 output columns"),
    ("spliceai_snv",       "SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz",              "vcf",        True,  "SpliceAI plugin"),
    ("spliceai_indel",     "SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz",            "vcf",        True,  "SpliceAI plugin"),
    ("revel",              "REVEL/new_tabbed_revel_grch38.tsv.gz",                      "tsv",        True,  "REVEL plugin"),
    ("loeuf",              "gnomAD_constraints/loeuf_dataset_grch38.tsv.gz",            "tsv",        True,  "LOEUF plugin"),
    ("vep_plugins",        "Plugins",                                                   "dir",        True,  "VEP plugin source directory"),
    # `homo_sapiens_merged`, not `homo_sapiens`: entrypoint.sh runs `vep --merged`, so
    # that is the cache VEP looks for. Requiring the plain one got both directions wrong
    # -- a reference directory holding only the merged cache, which is the correct one,
    # failed verification, while one holding only the plain cache was certified and then
    # failed at VEP. Invisible on a machine that has both, which is how it survived.
    ("vep_cache",          "homo_sapiens_merged",                                       "dir",        True,  "VEP offline cache (merged, as `vep --merged` requires)"),

    # Optional per entrypoint.sh: each is licence-gated or very large, and the
    # pipeline runs without it. Absent is not a defect. Present but damaged is,
    # because entrypoint.sh enables the annotation on an `-e` test alone and a
    # damaged file then silently empties the columns it feeds.
    ("cosmic",             "COSMIC/cosmic_cmc_grch38.vcf.gz",                           "vcf",        False, "optional, COSMIC licence"),
    ("cadd_snv",           "CADD/whole_genome_SNVs.tsv.gz",                             "tsv",        False, "optional, ~81 GB"),
    ("cadd_indel",         "CADD/gnomad.genomes.r4.0.indel.tsv.gz",                     "tsv",        False, "optional, ~1.2 GB"),
    ("dbnsfp",             "dbNSFP/dbNSFP5.0a_grch38.smart.gz",                         "tsv",        False, "optional, academic licence; built by utils/build_dbnsfp.py"),
    # Read whole by scripts/genie_annotation.py rather than queried by position, so
    # it is plain gzip with no .tbi and must not be checked as an indexed table --
    # doing so reported three defects against a sound file.
    ("genie_o4_counts",    "GENIE/genie_o4_counts.tsv.gz",                              "gztext",     False, "optional, GENIE data-use agreement; built by utils/build_genie_o4.py"),

    # SVIG-UK classifier resources. entrypoint.sh notes that a code whose resource
    # is missing is simply not applied, so these are optional too -- but a damaged
    # one changes a classification rather than withholding it, which is worse than
    # absence and so is still fatal when present.
    ("svig_cgc_roles",     "SVIG_UK/cgc_gene_roles.tsv",                                "text",       False, "optional, SVIG-UK gene roles"),
    ("svig_cosmic_o4",     "SVIG_UK/cosmic_o4_counts.tsv.gz",                           "gztext",     False, "optional, SVIG-UK O4 counts"),
    ("svig_chip_ab",       "SVIG_UK/gnomad_v41_chip_allele_balance.tsv.gz",             "gztext",     False, "optional, SVIG-UK CHIP allele balance"),
    ("svig_o5_variants",   "SVIG_UK/o5_reference_variants.tsv",                         "text",       False, "optional, SVIG-UK O5 reference variants"),
    ("oncokb_gene_roles",  "oncokb_gene_roles.tsv",                                     "text",       False, "optional, OncoKB gene roles"),
    ("svig_canonical",     "svig_uk_canonical_variants.tsv",                            "text",       False, "optional, SVIG-UK canonical variants"),
]


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _is_bgzf(path):
    """True if `path` is BGZF (bgzip), not merely gzip.

    tabix can only index BGZF. A file compressed with plain `gzip` is readable by
    every tool that decompresses it and unqueryable by the one that matters, so a
    reference table rebuilt with `gzip` instead of `bgzip` annotates nothing and
    reports no error. BGZF marks itself with an 'BC' extra subfield in the gzip
    header, which is what is looked for here.
    """
    try:
        with open(path, "rb") as fh:
            head = fh.read(18)
    except OSError:
        return False
    if len(head) < 18 or head[0:2] != b"\x1f\x8b":
        return False
    flg = head[3]
    if not flg & 0x04:          # FEXTRA not set -> cannot be BGZF
        return False
    xlen = struct.unpack("<H", head[10:12])[0]
    try:
        with open(path, "rb") as fh:
            fh.seek(12)
            extra = fh.read(xlen)
    except OSError:
        return False
    i = 0
    while i + 4 <= len(extra):
        si1, si2, slen = extra[i], extra[i + 1], struct.unpack("<H", extra[i + 2:i + 4])[0]
        if si1 == ord("B") and si2 == ord("C"):
            return True
        i += 4 + slen
    return False


def _read_tbi_config(idx):
    """Read the column layout out of a .tbi header.

    The probe query has to be built from the columns the index was actually built
    with, not from an assumption. They are not 1 and 2 across the reference set:

        clinvar / civic / spliceai / hotspots / cosmic   fmt=2 (VCF)  seq=1 beg=2
        CADD / dbNSFP                                    fmt=0        seq=1 beg=2
        REVEL                                            fmt=0        seq=1 beg=3
        LOEUF                                            fmt=0        seq=56 beg=57

    Guessing columns 1 and 2 took LOEUF's first field -- a gene symbol -- as a contig
    and queried `OR4F5:79501`, which returns nothing. That is indistinguishable from
    the real CIViC defect this script exists for, and a check that cries wolf on a
    sound file is worse than no check. So the layout is read from the index.

    Returns a dict, or None if the header cannot be read.
    """
    try:
        with gzip.open(idx, "rb") as fh:
            head = fh.read(36)
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error):
        return None
    if len(head) < 36 or head[:4] != b"TBI\x01":
        return None
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<8i", head[4:36])
    return {
        "n_ref": n_ref, "format": fmt,
        "col_seq": col_seq, "col_beg": col_beg, "col_end": col_end,
        "meta": chr(meta) if 0 < meta < 128 else "#",
        "skip": skip,
    }


def _open_text(path):
    return gzip.open(path, "rt", errors="replace") if path.endswith(".gz") else open(path, "r", errors="replace")


def _digest(path, size, max_hash_bytes):
    """SHA-256 of the whole file, or of head+tail+size when it is too big to read.

    Returns (kind, hexdigest). `kind` is recorded in the manifest so that a partial
    digest can never be read as a full one -- a verify that compared a full digest
    against a partial one would report a mismatch on a sound file, and a verify that
    silently accepted either would be weaker than it claims to be.
    """
    h = hashlib.sha256()
    if size <= max_hash_bytes:
        with open(path, "rb") as fh:
            while True:
                block = fh.read(HASH_BLOCK)
                if not block:
                    break
                h.update(block)
        return "sha256", h.hexdigest()

    # head + tail + exact size. Truncation, a quota stub and a partially written
    # build all change at least one of the three.
    edge = min(HASH_EDGE_BYTES, size // 2)
    with open(path, "rb") as fh:
        h.update(fh.read(edge))
        fh.seek(size - edge)
        h.update(fh.read(edge))
    h.update(str(size).encode())
    return "sha256-head%d-tail%d-size" % (edge, edge), h.hexdigest()


def _stub_check(path, defects):
    """Report a file whose first bytes are an error page rather than data."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(2048)
    except OSError as exc:
        defects.append("unreadable: %s" % exc)
        return
    if not head:
        defects.append("file is empty (0 bytes)")
        return
    for marker in STUB_MARKERS:
        if marker in head:
            defects.append(
                "first 2 KiB contain %r -- this looks like an error page or quota "
                "stub saved as data, not the reference file" % marker.decode("latin-1")
            )
            return


def _run_tabix(tabix, args, timeout=120):
    """Run tabix and return (returncode, stdout, stderr). Never raises on failure."""
    try:
        proc = subprocess.run([tabix] + args, capture_output=True, text=True, timeout=timeout)
        return proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired:
        return 124, "", "tabix timed out after %ss" % timeout
    except OSError as exc:
        return 127, "", str(exc)


# ---------------------------------------------------------------------------
# Structural inspection, one function per kind of file
#
# Each returns (facts, defects): `facts` are the cheap properties recorded in the
# manifest and compared on verify; `defects` are the things that are wrong now,
# independent of any manifest.
# ---------------------------------------------------------------------------

def inspect_header_and_records(path, kind, scan, defects, seq_col=0):
    """Read the leading comment block and, when `scan`, the whole file.

    Returns (header_lines, header_last, first_record, n_records, contigs).

    `seq_col` is the zero-based column holding the sequence name. It is 0 for a VCF,
    where CHROM is the first field, but a tabix-indexed TSV puts it wherever its index
    says: LOEUF's is column 56, and its *first* column is a gene symbol. Taking column
    0 regardless recorded 18,143 gene names as this file's contig list -- a fact that
    was useless (it said nothing about which chromosomes the file covers, so a release
    that dropped the sex chromosomes looked no different from one that did not) and
    expensive (a 115 KiB manifest entry). The caller reads the column from the index.

    `seq_col=None` means the file has no sequence column to read -- the SVIG-UK tables
    are keyed by gene and never queried by position -- so no contig list is collected
    and none is expected. Without this, the implausible-count check below would fire on
    every one of them and call a sound file damaged.
    """
    header_lines = []
    first_record = None
    n_records = None
    contigs = []

    try:
        with _open_text(path) as fh:
            read_bytes = 0
            for line in fh:
                read_bytes += len(line)
                if line.startswith("#"):
                    header_lines.append(line.rstrip("\n"))
                    if read_bytes > MAX_HEADER_BYTES:
                        defects.append(
                            "no data record found in the first %d MiB -- the file is "
                            "all header, or has no newline separating header from data"
                            % (MAX_HEADER_BYTES // (1024 * 1024))
                        )
                        return header_lines, None, None, None, []
                    continue
                first_record = line.rstrip("\n")
                break

            if first_record is None:
                defects.append("no data records: the file contains a header and nothing else")
                return header_lines, None, None, 0, []

            if scan:
                n_records = 1
                seen = {}

                def _seq(line):
                    f = line.split("\t")
                    return f[seq_col] if seq_col < len(f) else ""

                order = []
                if seq_col is not None:
                    chrom = _seq(first_record)
                    seen[chrom] = None
                    order = [chrom]
                for line in fh:
                    if line.startswith("#"):
                        # A comment after data starts is how a mis-concatenation of two
                        # files shows up: the second file's header lands mid-record-stream.
                        defects.append(
                            "comment line found after data began -- two files look "
                            "concatenated: %r" % line[:80].rstrip("\n")
                        )
                        continue
                    if not line.strip():
                        continue
                    n_records += 1
                    if seq_col is None:
                        continue
                    chrom = _seq(line.rstrip("\n"))
                    if chrom not in seen:
                        seen[chrom] = None
                        order.append(chrom)
                contigs = order
                # A reference file has tens of sequences, not thousands. Thousands
                # means this is not a sequence column at all -- the mistake above --
                # so say so rather than storing the list and calling it contigs.
                if seq_col is not None and len(contigs) > MAX_PLAUSIBLE_CONTIGS:
                    defects.append(
                        "column %d yielded %d distinct values, which is too many to be "
                        "sequence names -- the index's sequence column looks "
                        "misidentified, so no contig list was recorded"
                        % (seq_col + 1, len(contigs))
                    )
                    contigs = []
    except (EOFError, OSError, gzip.BadGzipFile, zlib.error) as exc:
        # A truncated .gz raises here rather than returning short data, which is the
        # one failure mode gzip does report -- so say so precisely instead of
        # letting it surface as an unhandled traceback that looks like a bug in this
        # script rather than damage in the reference file.
        defects.append("decompression failed (file is truncated or corrupt): %s: %s"
                       % (type(exc).__name__, exc))
        return header_lines, None, first_record, None, []

    return header_lines, (header_lines[-1] if header_lines else None), first_record, n_records, contigs


def check_fused_header(header_lines, defects):
    """Catch the defect that motivated this script.

    The CIViC VCF's eight header lines and its first data record are one physical
    line, so `header_lines` is 1, that line carries eight '##' markers, and the
    '#CHROM' column line is in the middle of it instead of starting a line of its
    own. All three are checked, because each is independently diagnostic and
    together they name the fault precisely in the run log.
    """
    if len(header_lines) == 1:
        defects.append(
            "header is a single physical line -- VCF meta-information lines are not "
            "parseable and the first data record is unreachable (this is the CIViC "
            "defect: see the module docstring)"
        )
    for i, line in enumerate(header_lines):
        markers = line.count("##")
        if markers > 1:
            defects.append(
                "header line %d carries %d '##' markers -- meta-information lines are "
                "fused together with no newline between them" % (i + 1, markers)
            )
        if "#CHROM" in line and not line.startswith("#CHROM"):
            defects.append(
                "header line %d contains '#CHROM' but does not start with it -- the "
                "column-header line is fused to whatever precedes it" % (i + 1)
            )


def _probe_from_index_columns(path, cfg, defects, facts):
    """Take a query position out of the file's first indexed data row.

    The position is read from the file itself rather than hard-coded, so the probe is
    always a position the file asserts it contains. A query that then returns nothing
    proves the file is unqueryable -- which is exactly the CIViC failure, where the
    first record is swallowed by a header line that never ended.

    `cfg` is the .tbi header, so `skip` leading lines and `meta`-prefixed comments are
    honoured and the contig/position are taken from the columns the index was built on.
    """
    meta = cfg["meta"]
    skip = max(0, cfg["skip"])
    col_seq, col_beg = cfg["col_seq"], cfg["col_beg"]
    if col_seq < 1 or col_beg < 1:
        defects.append("the .tbi header declares sequence column %d and begin column %d; "
                       "neither can be read, so the index cannot be probed"
                       % (col_seq, col_beg))
        return None
    try:
        with _open_text(path) as fh:
            skipped = 0
            for line in fh:
                if skipped < skip:
                    skipped += 1
                    continue
                if line.startswith(meta):
                    continue
                if not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                facts["probe_columns"] = [col_seq, col_beg]
                if len(cols) < max(col_seq, col_beg):
                    defects.append("first indexed data row has %d columns but the .tbi "
                                   "was built on columns %d/%d -- the index and the data "
                                   "do not describe the same table"
                                   % (len(cols), col_seq, col_beg))
                    return None
                chrom, pos = cols[col_seq - 1], cols[col_beg - 1]
                if not pos.isdigit():
                    defects.append("first indexed data row's position column %d is %r, "
                                   "not an integer -- the row is damaged or the index "
                                   "describes a different table" % (col_beg, pos[:40]))
                    return None
                return (chrom, int(pos))
    except (EOFError, OSError, gzip.BadGzipFile, zlib.error) as exc:
        defects.append("decompression failed while reading a probe row: %s: %s"
                       % (type(exc).__name__, exc))
        return None
    return None


# A genomic record's opening two fields, used to recover the data record that a fused
# header line swallowed. Deliberately loose: it only has to find a contig-looking token
# followed by a tab and an integer somewhere in the tail of the fused line.
_EMBEDDED_RECORD = re.compile(r"(?:^|\t|[A-Za-z])((?:chr)?(?:\d{1,2}|[XYxy]|MT?|mt?))\t(\d+)\t")


def _probe_record_hidden_in_header(header_line, path, tabix, defects, facts):
    """Show that the record a fused header line swallowed cannot be queried.

    The fusion checks already fail the file, so this is not what decides the verdict.
    It is here because the fusion's *consequence* is the part that matters clinically
    and the part nobody believes until it is printed: the CIViC VCF's first record
    (ERRFI1 E384*, 1:8013449) is in the file, is indexed, and comes back empty, so a
    report on that variant says "not in CIViC" rather than failing. Naming the lost
    record in the run log is the difference between "this file is odd" and "this file
    silently drops variant X".
    """
    at = header_line.find("#CHROM")
    if at < 0:
        return
    m = _EMBEDDED_RECORD.search(header_line, at)
    if not m:
        return
    chrom, pos = m.group(1), int(m.group(2))
    facts["record_hidden_in_header"] = [chrom, pos]
    if tabix is None:
        return
    rc, out, _err = _run_tabix(tabix, [path, "%s:%d-%d" % (chrom, pos, pos)])
    hits = [l for l in out.splitlines() if l.strip() and not l.startswith("#")] if rc == 0 else []
    if not hits:
        defects.append(
            "the data record fused into the header (%s:%d) cannot be retrieved: tabix "
            "returns nothing for it because the line it lives on starts with '#' and is "
            "skipped as a comment. Every lookup of that variant reports 'not found'"
            % (chrom, pos))


def inspect_vcf(path, tabix, scan, defects, facts, warnings):
    if not _is_bgzf(path):
        defects.append("not BGZF-compressed -- gzip can read it but tabix cannot index it")

    header_lines, _, first_record, n_records, contigs = inspect_header_and_records(
        path, "vcf", scan, defects)

    facts["header_lines"] = len(header_lines)
    facts["records"] = n_records
    facts["contigs"] = contigs or None

    check_fused_header(header_lines, defects)
    for line in header_lines:
        if "#CHROM" in line and not line.startswith("#CHROM"):
            _probe_record_hidden_in_header(line, path, tabix, defects, facts)

    if len(header_lines) < 2:
        # fileformat + #CHROM is the minimum a VCF can legally have.
        if len(header_lines) != 1:   # the 1 case is already reported above, precisely
            defects.append("VCF header has %d lines; a VCF needs at least "
                           "##fileformat and #CHROM" % len(header_lines))
    else:
        if not header_lines[0].startswith("##fileformat=VCF"):
            defects.append("first line is not ##fileformat=VCF: %r" % header_lines[0][:80])
        if not header_lines[-1].startswith("#CHROM"):
            defects.append("last header line is not the #CHROM column line: %r"
                           % header_lines[-1][:80])

    n_info = sum(1 for l in header_lines if l.startswith("##INFO="))
    facts["info_declarations"] = n_info
    if n_info == 0:
        # VEP --custom pulls named INFO fields. With no parseable ##INFO lines it
        # extracts nothing and reports no error -- empty columns indistinguishable
        # from "variant not in this database".
        defects.append("no parseable ##INFO declarations -- VEP --custom has nothing to extract")

    if first_record:
        cols = first_record.split("\t")
        facts["first_record_columns"] = len(cols)
        if len(cols) < 8:
            defects.append("first data record has %d tab-separated columns; a VCF needs "
                           "at least 8: %r" % (len(cols), first_record[:120]))
        elif not cols[1].isdigit():
            defects.append("first data record's POS is not an integer (%r) -- the record "
                           "is damaged or fused to the header" % cols[1][:40])

    idx = path + ".tbi"
    cfg = _read_tbi_config(idx)
    if cfg is not None:
        facts["index_format"] = cfg["format"]
        if cfg["format"] != 2:
            # A VCF indexed as a generic table (fmt 0) is indexed on the wrong
            # coordinate convention, so queries near indels return the wrong window.
            defects.append("the .tbi says format %d; a VCF index must be format 2 "
                           "(built with `tabix -p vcf`)" % cfg["format"])
    probe = _probe_from_index_columns(path, cfg, defects, facts) if cfg else None

    _check_index(path, idx, tabix, probe, defects, facts, "VCF", warnings, cfg)
    return facts, defects


def inspect_tsv(path, tabix, scan, defects, facts, warnings):
    """A tabix-indexed TSV score table (REVEL, LOEUF, CADD, dbNSFP, GENIE)."""
    if not _is_bgzf(path):
        defects.append("not BGZF-compressed -- gzip can read it but tabix cannot index it")

    # Read the index first: it is what says which column holds the sequence name, and a
    # tabix TSV puts it wherever the builder chose. LOEUF's is column 56 while its first
    # column is a gene symbol, so scanning column 1 recorded 18,143 gene names as the
    # contig list -- which said nothing about the chromosomes the file covers, and so
    # could not show that gnomAD v4.1 omits chrX and chrY.
    idx = path + ".tbi"
    cfg = _read_tbi_config(idx)
    seq_col = (cfg["col_seq"] - 1) if cfg and cfg.get("col_seq") else 0

    header_lines, _, first_record, n_records, contigs = inspect_header_and_records(
        path, "tsv", scan, defects, seq_col=seq_col)

    facts["header_lines"] = len(header_lines)
    facts["records"] = n_records
    facts["contigs"] = contigs or None

    # These tables carry one '#'-prefixed column line. More than one is the
    # concatenation signature; the '##' fusion test catches a mis-joined header the
    # same way it does for VCFs.
    for i, line in enumerate(header_lines):
        if "\t#" in line:
            defects.append("header line %d has a field starting with '#' -- a second "
                           "header line looks fused onto this one: %r" % (i + 1, line[:100]))

    if first_record:
        cols = first_record.split("\t")
        facts["columns"] = len(cols)
        if header_lines:
            hcols = len(header_lines[-1].lstrip("#").split("\t"))
            facts["header_columns"] = hcols
            if hcols != len(cols):
                defects.append("header declares %d columns, first data row has %d -- the "
                               "column mapping of every annotated value from this file "
                               "is shifted" % (hcols, len(cols)))

    if cfg is not None:
        facts["index_format"] = cfg["format"]
        facts["index_columns"] = [cfg["col_seq"], cfg["col_beg"], cfg["col_end"]]
    probe = _probe_from_index_columns(path, cfg, defects, facts) if cfg else None

    _check_index(path, idx, tabix, probe, defects, facts, "TSV", warnings, cfg)
    return facts, defects


def inspect_gztext(path, tabix, scan, defects, facts, warnings):
    """A gzipped table that is read whole rather than queried (SVIG-UK resources).

    No index, so the only structural facts available are that it decompresses to the
    end and that its columns are consistent -- which is what a truncated or
    half-written build breaks.

    These tables have no comment marker on their header row (GENIE's begins `gene`),
    so `records` counts the header among the rows. That is left as it is rather than
    guessed at: the number is being compared against the same number computed the same
    way at generate time, and a wrong-by-one count that is wrong consistently still
    detects a changed or truncated file, whereas a heuristic "first row looks like a
    header" would detect it only sometimes.
    """
    # No index and no coordinates: these tables are keyed by gene and read whole, so
    # there is no sequence column to collect contigs from.
    header_lines, _, first_record, n_records, _ = inspect_header_and_records(
        path, "gztext", True, defects, seq_col=None)
    facts["header_lines"] = len(header_lines)
    facts["records"] = n_records
    if first_record:
        facts["columns"] = len(first_record.split("\t"))
    return facts, defects


def inspect_text(path, tabix, scan, defects, facts, warnings):
    """A plain text table or chain file: line count, column count, stub check."""
    lines = 0
    first = None
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                if first is None and line.strip():
                    first = line.rstrip("\n")
                lines += 1
    except OSError as exc:
        defects.append("unreadable: %s" % exc)
        return facts, defects
    facts["lines"] = lines
    if lines == 0:
        defects.append("file is empty")
    if first is not None:
        facts["columns"] = len(first.split("\t"))
    if path.endswith(".chain"):
        # A chain file that lost its first line still lifts over, silently and wrongly.
        if first is None or not first.startswith("chain"):
            defects.append("first non-blank line of a liftover chain file does not start "
                           "with 'chain': %r" % ((first or "")[:80]))
    return facts, defects


def inspect_fasta(path, tabix, scan, defects, facts, warnings):
    """A FASTA plus its .fai, checked for mutual consistency.

    scripts/check_reference.py reads the genome through the .fai by seeking to a
    recorded offset. Nothing confirms that the offsets belong to this FASTA, so a
    .fai left behind from a different build reads the wrong bases at every position
    and the REF-allele check it powers then reports mismatches for a sound VCF -- or,
    worse, agreement for an unsound one.
    """
    size = os.path.getsize(path)
    fai = path + ".fai"
    facts["fai_present"] = os.path.exists(fai)
    if not facts["fai_present"]:
        defects.append(".fai index is missing -- VEP and scripts/check_reference.py "
                       "cannot read this FASTA")
        return facts, defects

    if os.path.getmtime(fai) < os.path.getmtime(path):
        defects.append(".fai is older than the FASTA -- the index is stale and its "
                       "offsets point into the previous file")

    entries = []
    try:
        with open(fai) as fh:
            for lineno, line in enumerate(fh, 1):
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    defects.append(".fai line %d has %d fields, expected 5" % (lineno, len(parts)))
                    continue
                name, length, offset, linebases, linewidth = parts[:5]
                entries.append((name, int(length), int(offset), int(linebases), int(linewidth)))
    except (OSError, ValueError) as exc:
        defects.append(".fai is unreadable or malformed: %s" % exc)
        return facts, defects

    facts["fai_contigs"] = len(entries)
    facts["fai_total_bases"] = sum(e[1] for e in entries)
    facts["size"] = size

    if not entries:
        defects.append(".fai is empty")
        return facts, defects

    # The last contig's data must end at (or just before) the end of the file. If the
    # FASTA was truncated or the .fai belongs to another build, this is where it shows.
    #
    # The arithmetic has to count the final, partial sequence line explicitly. Rounding
    # the whole contig up to `ceil(length / linebases) * linewidth` charges a full line's
    # newline padding to a final line that is shorter than `linebases`, which on hg38
    # over-states the file by 10 bytes and reported the genome as truncated when it is
    # not -- a false positive on the one file every run depends on, which is the fastest
    # way to get a check switched off.
    name, length, offset, linebases, linewidth = entries[-1]
    expected_end = offset + _fasta_contig_bytes(length, linebases, linewidth)
    facts["fai_expected_end"] = expected_end
    # A FASTA whose final line has no terminating newline is legal, so allow the
    # end to fall short by the line terminator's width.
    slack = linewidth - linebases
    if expected_end > size:
        defects.append(".fai describes %d bytes but the FASTA is %d bytes -- the FASTA is "
                       "truncated or the index belongs to a different file"
                       % (expected_end, size))
    elif expected_end < size - slack:
        defects.append(".fai describes %d bytes but the FASTA is %d bytes -- there are "
                       "%d unindexed bytes after the last contig, which is what a second "
                       "FASTA concatenated onto the first looks like"
                       % (expected_end, size, size - expected_end))

    # Spot-check that the first contig's recorded offset really lands on sequence
    # preceded by its own '>' header. A consistent-looking .fai from another assembly
    # passes the arithmetic above and fails this.
    try:
        with open(path, "rb") as fh:
            first = entries[0]
            back = min(first[2], 4096)
            fh.seek(first[2] - back)
            before = fh.read(back)
            at = fh.read(min(first[3], 200))
        if b">" + first[0].encode() not in before:
            defects.append(".fai says contig %r begins at byte %d but no '>%s' header "
                           "precedes that offset -- the index does not belong to this FASTA"
                           % (first[0], first[2], first[0]))
        if at and not all(chr(c).upper() in "ACGTNRYKMSWBDHV" for c in at.strip()):
            defects.append(".fai offset for contig %r does not land on sequence: %r"
                           % (first[0], at[:40]))
    except OSError as exc:
        defects.append("could not seek into the FASTA to check the .fai: %s" % exc)

    return facts, defects


def _fasta_contig_bytes(length, linebases, linewidth):
    """On-disk byte span of one contig's sequence, final partial line included."""
    full, rem = divmod(length, linebases)
    return full * linewidth + (rem + (linewidth - linebases) if rem else 0)


def inspect_dir(path, tabix, scan, defects, facts, warnings):
    """A directory resource (VEP cache, VEP plugins): present and non-empty."""
    try:
        n = len(os.listdir(path))
    except OSError as exc:
        defects.append("unreadable directory: %s" % exc)
        return facts, defects
    facts["entries"] = n
    if n == 0:
        defects.append("directory is empty -- a failed extraction leaves the directory "
                       "behind, and entrypoint.sh only tests that it exists")
    return facts, defects


def _check_index(path, idx, tabix, probe, defects, facts, label, warnings, cfg=None):
    """Index presence, freshness, contig list, and a query that must return a record."""
    facts["index_present"] = os.path.exists(idx)
    if not facts["index_present"]:
        defects.append("%s index %s is missing" % (label, os.path.basename(idx)))
        return

    data_mtime = os.path.getmtime(path)
    idx_mtime = os.path.getmtime(idx)
    facts["index_newer_than_data"] = idx_mtime >= data_mtime
    if idx_mtime < data_mtime:
        # The data was re-downloaded or rebuilt and the old index kept. Queries then
        # read offsets into bytes that have moved: the file answers with the wrong
        # records, or with none, and never with an error.
        #
        # A warning rather than a defect, because the honest signal here is weak. CADD
        # ships its own .tbi, which `get_ref_files.sh` downloads *before* the 81 GB data
        # file it belongs to, so on this machine both CADD indexes are legitimately
        # older than their data and the probe query against them succeeds. Making this
        # fatal would fail a sound reference set. The probe below is the authoritative
        # test of whether the index agrees with the data; this only says look closer.
        warnings.append("index is older than the data (%s vs %s). Usually harmless -- a "
                        "shipped index downloaded before its data file does this -- but "
                        "it is also what a re-downloaded file with a stale index looks "
                        "like, so confirm the probe query below passed"
                        % (time.strftime("%Y-%m-%d %H:%M", time.localtime(idx_mtime)),
                           time.strftime("%Y-%m-%d %H:%M", time.localtime(data_mtime))))

    if tabix is None:
        defects.append("tabix is not available, so the index could not be queried")
        return

    rc, out, err = _run_tabix(tabix, ["-l", path])
    if rc != 0:
        defects.append("tabix -l failed (rc=%d): %s" % (rc, err.strip()[:200]))
    else:
        facts["index_contigs"] = [c for c in out.split() if c]

    if probe is None:
        defects.append("no probe position could be taken from the file, so no query was "
                       "made against the index")
        return

    chrom, pos = probe
    facts["probe"] = [chrom, pos]
    rc, out, err = _run_tabix(tabix, [path, "%s:%d-%d" % (chrom, pos, pos)])
    if rc != 0:
        defects.append("tabix query %s:%d failed (rc=%d): %s" % (chrom, pos, rc, err.strip()[:200]))
        return
    hits = [l for l in out.splitlines() if l.strip() and not l.startswith("#")]
    facts["probe_hits"] = len(hits)
    if not hits:
        # This is the check that catches the CIViC file. The probe position is read
        # out of the file's own first data record, so an empty result means the file
        # cannot return a record it demonstrably contains -- it is unqueryable, and
        # every lookup against it comes back "not found".
        defects.append(
            "tabix query for %s:%d returns nothing, but that position comes from this "
            "file's own first data record -- records in this file are unreachable and "
            "every annotation from it will be silently empty" % (chrom, pos))
        return

    # At least one returned row must actually cover the probe. Checked through the
    # index's own columns and as an interval, not as an equality on column 2: LOEUF is
    # indexed on columns 56/57/58 as gene intervals, so a query at the first gene's
    # start legitimately returns a row whose second column (a gene id) is not the
    # position at all. Comparing column 2 reported "the index does not agree with the
    # data" against a sound file.
    if cfg is None:
        return
    c_seq, c_beg, c_end = cfg["col_seq"], cfg["col_beg"], cfg["col_end"]
    covering = 0
    for line in hits:
        cols = line.split("\t")
        if len(cols) < max(c_seq, c_beg):
            continue
        if cols[c_seq - 1] != chrom:
            continue
        try:
            beg = int(cols[c_beg - 1])
        except ValueError:
            continue
        end = beg
        if c_end > 0 and len(cols) >= c_end:
            try:
                end = int(cols[c_end - 1])
            except ValueError:
                end = beg
        # VCF (format 2) and BED-like generic formats differ by one in where a record
        # starts, so allow a single base of slack rather than encoding each convention.
        if beg - 1 <= pos <= max(end, beg):
            covering += 1
    facts["probe_covering_rows"] = covering
    if covering == 0:
        defects.append(
            "tabix query for %s:%d returned %d row(s), none of which covers that "
            "position on that contig -- the index does not describe this data. First "
            "row returned: %r" % (chrom, pos, len(hits), hits[0][:160]))


INSPECTORS = {
    "vcf": inspect_vcf,
    "tsv": inspect_tsv,
    "gztext": inspect_gztext,
    "text": inspect_text,
    "fasta": inspect_fasta,
    "dir": inspect_dir,
}


# ---------------------------------------------------------------------------
# Per-entry driver
# ---------------------------------------------------------------------------

def inspect_entry(refs, name, rel, kind, required, note, tabix,
                  max_hash_bytes, max_scan_bytes, quick):
    """Return a manifest-shaped record for one reference file, defects included."""
    path = os.path.join(refs, rel)
    rec = {
        "name": name,
        "path": rel,
        "kind": kind,
        "required": required,
        "note": note,
        "present": os.path.exists(path),
        "facts": {},
        "defects": [],
        # Kept apart from `defects` so that a finding whose signal is genuinely weak
        # -- a shipped index older than its own data file -- does not stop a clinical
        # run. Everything in `defects` is fatal; see "Fatal versus warning" above.
        "warnings": [],
    }
    if not rec["present"]:
        return rec

    defects = rec["defects"]
    facts = rec["facts"]
    warnings = rec["warnings"]

    if kind == "dir":
        INSPECTORS[kind](path, tabix, False, defects, facts, warnings)
        return rec

    size = os.path.getsize(path)
    facts["size"] = size
    if size == 0:
        defects.append("file is 0 bytes")
        return rec

    _stub_check(path, defects)

    scan = size <= max_scan_bytes
    facts["scanned"] = scan
    if not scan:
        # Recorded rather than left implicit: a reader of the manifest must be able
        # to tell "0 records" apart from "we did not count".
        facts["scan_skipped_reason"] = (
            "%d bytes exceeds --max-scan-bytes %d; record count and contig list come "
            "from the index instead" % (size, max_scan_bytes))

    INSPECTORS[kind](path, tabix, scan, defects, facts, warnings)

    if not quick:
        kind_label, digest = _digest(path, size, max_hash_bytes)
        facts["digest_kind"] = kind_label
        facts["digest"] = digest
    return rec


# ---------------------------------------------------------------------------
# generate
# ---------------------------------------------------------------------------

def cmd_generate(args):
    tabix = resolve_tabix(args)
    entries = select_entries(args)

    records = []
    for name, rel, kind, required, note in entries:
        print("[scan] %-20s %s" % (name, rel), flush=True)
        records.append(inspect_entry(args.refs, name, rel, kind, required, note, tabix,
                                     args.max_hash_bytes, args.max_scan_bytes, args.quick))

    manifest = {
        "manifest_version": MANIFEST_VERSION,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "generated_from": os.path.abspath(args.refs),
        "generator": "scripts/verify_reference_files.py",
        "max_hash_bytes": args.max_hash_bytes,
        "max_scan_bytes": args.max_scan_bytes,
        "quick": bool(args.quick),
        # Spelled out in the file itself, because a manifest outlives the knowledge of
        # how it was made and the first question anyone asks of a failing verify is
        # whether the manifest or the data is wrong.
        "bootstrapping_warning": (
            "This manifest records the state of the directory above at the time shown. "
            "A manifest generated from a damaged directory would record the damage as "
            "correct, so any entry that failed its structural checks carries a "
            "'defects' list and verify fails that file regardless of whether its bytes "
            "still match. Do not clear a 'defects' list to make verify pass: rebuild "
            "the file."
        ),
        "files": records,
    }

    bad = [r for r in records if r["defects"]]
    warned = [r for r in records if r["warnings"] and not r["defects"]]
    missing_required = [r for r in records if r["required"] and not r["present"]]

    with open(args.manifest, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=False)
        fh.write("\n")
    print("\n[OK] manifest written: %s" % args.manifest)

    for r in warned:
        for w in r["warnings"]:
            print("[WARN] %-20s %s" % (r["name"], w))

    if bad or missing_required:
        print("\n" + "=" * 78)
        print("REFUSING TO CERTIFY THIS REFERENCE DIRECTORY")
        print("=" * 78)
        for r in missing_required:
            print("  MISSING  %-20s %s" % (r["name"], r["path"]))
        for r in bad:
            print("  DEFECTIVE %-19s %s" % (r["name"], r["path"]))
            for d in r["defects"]:
                print("      - %s" % d)
        print("\nThe manifest above records these as defects, not as the expected state,")
        print("so `verify` will fail these files even against this manifest. Fix the")
        print("files and regenerate. See the module docstring, 'The bootstrapping hazard'.")
        return EXIT_FATAL
    return EXIT_OK


# ---------------------------------------------------------------------------
# verify
# ---------------------------------------------------------------------------

# Facts compared against the manifest. Deliberately not every recorded fact: mtimes
# and probe hit counts are properties of the moment, not of the file's identity.
COMPARED_FACTS = (
    "size", "header_lines", "records", "columns", "header_columns",
    "info_declarations", "first_record_columns", "contigs",
    "fai_contigs", "fai_total_bases", "lines", "entries",
)


def cmd_verify(args):
    tabix = resolve_tabix(args)
    entries = select_entries(args)

    manifest = None
    expected = {}
    if args.manifest and os.path.exists(args.manifest):
        try:
            with open(args.manifest) as fh:
                manifest = json.load(fh)
        except (OSError, ValueError) as exc:
            print("ERROR: manifest %s is unreadable: %s" % (args.manifest, exc))
            return EXIT_CANNOT_RUN
        expected = {r["name"]: r for r in manifest.get("files", [])}
    elif args.manifest:
        print("ERROR: manifest not found: %s" % args.manifest)
        print("       Run `generate` first, or omit --manifest to run the structural")
        print("       checks alone (which still catch a mis-framed or unqueryable file).")
        return EXIT_CANNOT_RUN

    findings = []     # (severity, name, message)
    checked = 0

    for name, rel, kind, required, note in entries:
        exp = expected.get(name)
        got = inspect_entry(args.refs, name, rel, kind, required, note, tabix,
                            (manifest or {}).get("max_hash_bytes", args.max_hash_bytes),
                            (manifest or {}).get("max_scan_bytes", args.max_scan_bytes),
                            args.quick or (manifest or {}).get("quick", False))
        checked += 1

        if not got["present"]:
            if required:
                findings.append((SEVERITY_FATAL, name, "required file missing: %s" % rel))
            elif exp is not None and exp.get("present"):
                # It was there when the manifest was made and is not there now. The
                # pipeline will quietly run without the annotation, so say it out loud.
                findings.append((SEVERITY_WARN, name,
                                 "optional file was present when the manifest was "
                                 "generated and is now missing: %s" % rel))
            else:
                findings.append((SEVERITY_INFO, name, "optional file not provided (%s)" % note))
            continue

        # Structural defects are fatal for required files and for optional files that
        # are present, because entrypoint.sh enables an optional annotation on an
        # existence test alone: a damaged optional file is used, not skipped.
        for d in got["defects"]:
            findings.append((SEVERITY_FATAL, name, d))
        for w in got["warnings"]:
            findings.append((SEVERITY_WARN, name, w))

        if exp is None:
            if manifest is not None:
                findings.append((SEVERITY_WARN, name,
                                 "not listed in the manifest -- it was added to the "
                                 "registry after the manifest was generated"))
            continue

        # A manifest entry that recorded defects stays failed. This is what stops a
        # manifest generated from the corrupt directory from blessing it: the bytes
        # match, every fact matches, and the file still fails.
        #
        # Only defects the live check did not already find are repeated, so that a
        # still-broken file is reported once rather than twice. A manifest defect with
        # no live counterpart is the interesting case and gets its own wording: the
        # file now looks sound but was not sound when the manifest was certified, so
        # the manifest cannot vouch for what is there now.
        live = set(got["defects"])
        for d in exp.get("defects", []):
            if d in live:
                continue
            findings.append((SEVERITY_FATAL, name,
                             "the manifest was generated from a defective copy of this "
                             "file and so cannot vouch for it; the recorded defect no "
                             "longer reproduces, which means the file changed: %s. "
                             "Rebuild the file and regenerate the manifest." % d))

        ef, gf = exp.get("facts", {}), got["facts"]

        digest_changed = False
        if "digest" in ef and "digest" in gf:
            if ef.get("digest_kind") != gf.get("digest_kind"):
                findings.append((SEVERITY_WARN, name,
                                 "digest kind differs (manifest %s, now %s) -- the two "
                                 "cannot be compared; re-run with matching "
                                 "--max-hash-bytes" % (ef.get("digest_kind"), gf.get("digest_kind"))))
            elif ef["digest"] != gf["digest"]:
                digest_changed = True

        fact_diffs = []
        for key in COMPARED_FACTS:
            if key in ef and key in gf and ef[key] != gf[key]:
                fact_diffs.append("%s: manifest %s, now %s" % (key, _short(ef[key]), _short(gf[key])))

        if digest_changed or fact_diffs:
            detail = "; ".join(fact_diffs) if fact_diffs else "contents differ"
            if digest_changed:
                detail = "SHA-256 differs" + ("; " + detail if fact_diffs else "")
            # Structure is sound (no defects were appended above), so this is a
            # different file rather than a damaged one -- normally a reference update.
            # Recorded as a warning whatever --strict says, and escalated by report()
            # into its own exit status (3). Keeping the severity honest here means the
            # report still calls it what it is, and entrypoint.sh can tell "changed but
            # sound under a frozen reference set" apart from "damaged".
            findings.append((SEVERITY_WARN, name,
                             "this is not the file the manifest was generated from (%s). "
                             "Its structure is sound, so it is probably an updated "
                             "release -- regenerate the manifest to record it." % detail))

    return report(findings, checked, manifest, args)


def _short(value):
    if isinstance(value, list):
        return "%d items" % len(value) if len(value) > 6 else repr(value)
    return repr(value)


def report(findings, checked, manifest, args):
    fatal = [f for f in findings if f[0] == SEVERITY_FATAL]
    warn = [f for f in findings if f[0] == SEVERITY_WARN]
    info = [f for f in findings if f[0] == SEVERITY_INFO]

    print("=" * 78)
    print("SMART reference verification")
    print("=" * 78)
    print("  reference directory : %s" % os.path.abspath(args.refs))
    print("  manifest            : %s" % (args.manifest if manifest else
                                          "(none -- structural checks only)"))
    if manifest:
        print("  manifest generated  : %s from %s"
              % (manifest.get("generated_at"), manifest.get("generated_from")))
    print("  files checked       : %d" % checked)
    print("")

    for label, group in (("FATAL", fatal), ("WARNING", warn), ("INFO", info)):
        if not group:
            continue
        if label == "INFO" and not args.verbose:
            print("INFO    : %d optional source(s) not provided (use --verbose to list)" % len(group))
            continue
        print("%s:" % label)
        for _sev, name, msg in group:
            print("  [%s] %s" % (name, msg))
        print("")

    if fatal:
        print("=" * 78)
        print("RESULT: FAIL -- %d fatal finding(s). Do not annotate against this "
              "reference directory." % len(fatal))
        print("A damaged reference does not make the pipeline fail; it makes the")
        print("pipeline write empty or wrong columns that look like real results.")
        print("=" * 78)
        return EXIT_FATAL

    if warn and args.strict:
        print("RESULT: FAIL (--strict) -- %d warning(s) treated as fatal." % len(warn))
        return EXIT_WARNINGS_STRICT

    print("RESULT: PASS%s" % (" with %d warning(s)" % len(warn) if warn else ""))
    return EXIT_OK


# ---------------------------------------------------------------------------
# Argument handling
# ---------------------------------------------------------------------------

def select_entries(args):
    entries = [(n, r, k, req, note) for (n, r, k, req, note) in REGISTRY]
    if args.only:
        wanted = set(args.only)
        unknown = wanted - {e[0] for e in entries}
        if unknown:
            raise SystemExit("ERROR: unknown --only name(s): %s\nKnown: %s"
                             % (", ".join(sorted(unknown)),
                                ", ".join(e[0] for e in entries)))
        entries = [e for e in entries if e[0] in wanted]
    return entries


def resolve_tabix(args):
    tabix = args.tabix or shutil.which("tabix")
    if tabix is None:
        # Not fatal by itself: the header, framing, digest and index-freshness checks
        # all still run and the missing query is reported per file as a defect. Being
        # unable to query is itself worth knowing about, so it is not silently skipped.
        print("[WARN] tabix not found on PATH -- index queries will be reported as "
              "unchecked. Pass --tabix /path/to/tabix.", file=sys.stderr)
    return tabix


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Verify SMART reference files against a manifest.",
        epilog="Exit codes: 0 pass, 1 fatal finding, 2 could not run, 3 warnings under --strict.")
    sub = parser.add_subparsers(dest="mode", required=True)

    for mode in ("generate", "verify"):
        p = sub.add_parser(mode)
        p.add_argument("--refs", required=True, help="reference directory root")
        p.add_argument("--manifest",
                       required=(mode == "generate"),
                       help="manifest path (verify: omit to run structural checks only)")
        p.add_argument("--tabix", default=None, help="path to tabix (default: from PATH)")
        p.add_argument("--only", nargs="+", default=None,
                       help="registry names to check instead of all")
        p.add_argument("--max-hash-bytes", type=int, default=DEFAULT_MAX_HASH_BYTES,
                       help="files larger than this get a head+tail+size digest "
                            "(default %d)" % DEFAULT_MAX_HASH_BYTES)
        p.add_argument("--max-scan-bytes", type=int, default=DEFAULT_MAX_SCAN_BYTES,
                       help="files larger than this are not decompressed to count "
                            "records (default %d)" % DEFAULT_MAX_SCAN_BYTES)
        p.add_argument("--quick", action="store_true",
                       help="skip digests; structural checks only")
        p.add_argument("--verbose", action="store_true", help="list INFO findings")
        if mode == "verify":
            p.add_argument("--strict", action="store_true",
                           help="treat a changed-but-sound file as fatal (frozen "
                                "validated reference set)")

    args = parser.parse_args(argv)
    if not hasattr(args, "strict"):
        args.strict = False

    if not os.path.isdir(args.refs):
        print("ERROR: reference directory not found: %s" % args.refs)
        return EXIT_CANNOT_RUN

    return cmd_generate(args) if args.mode == "generate" else cmd_verify(args)


if __name__ == "__main__":
    sys.exit(main())
