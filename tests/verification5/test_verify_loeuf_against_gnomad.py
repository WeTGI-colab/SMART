"""Tests for scripts/verify_loeuf_against_gnomad.py.

The script exists because the LOEUF table's only integrity check was that the download
was not an HTML error page, which would pass for a different gnomAD release, for a
truncation that still decompressed, and for any edited value. These tests assert it
catches each of those, and that it does not report a problem when there is none -- a
verifier that always passes is what left ONCOKB_VUS empty for years.

No network and no reference files: every fixture is built here.
"""

import csv
import gzip
import importlib.util
import io
import os
import sys
import urllib.error

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPT = os.path.join(_HERE, "..", "..", "scripts", "verify_loeuf_against_gnomad.py")
_spec = importlib.util.spec_from_file_location("verify_loeuf_against_gnomad", _SCRIPT)
vl = importlib.util.module_from_spec(_spec)
sys.modules["verify_loeuf_against_gnomad"] = vl
_spec.loader.exec_module(vl)


# The columns the real release carries, trimmed to what the script actually touches
# plus enough neighbours that a column-offset bug would show.
RELEASE_COLUMNS = ["gene", "gene_id", "transcript", "canonical", "mane_select",
                   "lof.oe_ci.upper", "mis.z_score", "chromosome"]


def release_rows(n=5):
    return [["GENE%d" % i, "ENSG%05d" % i, "ENST%05d" % i,
             "true" if i == 0 else "false", "true" if i == 0 else "false",
             "%.4e" % (0.1 * (i + 1)), "%.4e" % (i - 1), "chr%d" % (i + 1)]
            for i in range(n)]


def write_release(path, rows, columns=RELEASE_COLUMNS):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(columns)
        w.writerows(rows)
    return str(path)


def write_table(path, rows, columns=RELEASE_COLUMNS, contigs=None):
    """The indexed table: release columns, '#' on the first, coordinates appended."""
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["#" + columns[0]] + columns[1:] + ["chrom", "start", "end"])
        for i, rec in enumerate(rows):
            chrom = (contigs[i] if contigs else "chr%d" % (i + 1))
            w.writerow(list(rec) + [chrom, "100", "900"])
    return str(path)


class TestAgreement:
    def test_an_unaltered_table_passes(self, tmp_path, capsys):
        rows = release_rows()
        rel = write_release(tmp_path / "rel.tsv", rows)
        tab = write_table(tmp_path / "t.tsv.gz", rows)
        rc = vl.main([tab, "--release", "4.1", "--release-file", rel])
        out = capsys.readouterr().out
        assert rc == 0, out
        assert "[OK]" in out
        assert "equals gnomAD v4.1" in out

    def test_the_appended_index_columns_are_not_compared(self, tmp_path, capsys):
        """chrom/start/end have no counterpart in the release, so comparing them would
        make every row differ -- a verifier that fails on everything is as useless as
        one that passes on everything."""
        rows = release_rows()
        rel = write_release(tmp_path / "rel.tsv", rows)
        tab = write_table(tmp_path / "t.tsv.gz", rows)
        assert vl.main([tab, "--release-file", rel]) == 0
        out = capsys.readouterr().out
        assert "skipped the index columns chrom, end, start" in out
        # and it did compare the rest, rather than comparing nothing and declaring success
        assert "compared 8 columns" in out


class TestItCatchesWhatItExistsFor:
    def test_one_altered_score_fails(self, tmp_path, capsys):
        rows = release_rows()
        rel = write_release(tmp_path / "rel.tsv", rows)
        bad = [list(r) for r in rows]
        bad[2][6] = "9.9999e+00"            # mis.z_score, the value that drives O8
        tab = write_table(tmp_path / "t.tsv.gz", bad)
        rc = vl.main([tab, "--release-file", rel])
        out = capsys.readouterr().out
        assert rc == 1, out
        assert "different            : 1" in out
        assert "mis.z_score" in out
        assert "ENST00002" in out          # it names the row, so it can be investigated

    def test_an_altered_loeuf_fails(self, tmp_path, capsys):
        rows = release_rows()
        rel = write_release(tmp_path / "rel.tsv", rows)
        bad = [list(r) for r in rows]
        bad[0][5] = "0.0010e+00"
        tab = write_table(tmp_path / "t.tsv.gz", bad)
        assert vl.main([tab, "--release-file", rel]) == 1
        assert "lof.oe_ci.upper" in capsys.readouterr().out

    def test_a_truncated_table_is_not_silently_accepted(self, tmp_path, capsys):
        """Half a file still decompresses. The rows that survive all match, so a
        comparison that only looked at shared rows would pass it."""
        rows = release_rows(5)
        rel = write_release(tmp_path / "rel.tsv", rows)
        tab = write_table(tmp_path / "t.tsv.gz", rows[:2])
        rc = vl.main([tab, "--release-file", rel])
        out = capsys.readouterr().out
        # It passes the value comparison -- correctly, the values are right -- but the
        # count of release rows it does not hold is reported, which is the symptom.
        assert "only in release      : 3" in out
        assert rc == 0, out

    def test_a_different_schema_is_reported_as_such(self, tmp_path, capsys):
        """gnomAD v4.0 really does have a different column set from v4.1. Comparing
        against the wrong release must say so rather than report every row as differing."""
        rows = release_rows()
        rel = write_release(tmp_path / "rel.tsv", [r[:3] + r[5:7] for r in rows],
                            columns=["gene", "gene_id", "transcript",
                                     "lof.oe_ci.upper", "mis.z_score"])
        tab = write_table(tmp_path / "t.tsv.gz", rows)
        rc = vl.main([tab, "--release-file", rel])
        out = capsys.readouterr().out
        assert rc == 1
        assert "has no column" in out
        assert "canonical" in out

    def test_no_shared_transcript_is_a_failure_not_a_pass(self, tmp_path, capsys):
        rows = release_rows()
        other = [["X", "ENSG9", "ENSTZZZ%d" % i, "false", "false", "1", "1", "chr1"]
                 for i in range(3)]
        rel = write_release(tmp_path / "rel.tsv", other)
        tab = write_table(tmp_path / "t.tsv.gz", rows)
        rc = vl.main([tab, "--release-file", rel])
        out = capsys.readouterr().out
        assert rc == 1, out
        assert "share no transcript" in out

    def test_an_unreadable_table_exits_2(self, tmp_path, capsys):
        rc = vl.main([str(tmp_path / "nope.tsv.gz"), "--release-file", "x"])
        assert rc == 2
        assert "cannot read" in capsys.readouterr().out

    def test_a_table_without_a_transcript_column_exits_2(self, tmp_path, capsys):
        p = tmp_path / "t.tsv.gz"
        with gzip.open(p, "wt", newline="") as fh:
            fh.write("#gene\tvalue\nA\t1\n")
        assert vl.main([str(p), "--release-file", "x"]) == 2
        assert "no 'transcript' column" in capsys.readouterr().out


class TestUnreachableIsDistinctFromWrong:
    def test_network_failure_exits_3_not_1(self, tmp_path, capsys, monkeypatch):
        """An unverifiable reference and a wrong one need different answers: the first
        is a judgement for the caller, the second must stop the run."""
        rows = release_rows()
        tab = write_table(tmp_path / "t.tsv.gz", rows)

        def boom(*a, **k):
            raise urllib.error.URLError("no route to host")
        monkeypatch.setattr(vl.urllib.request, "urlopen", boom)
        rc = vl.main([tab, "--cache-dir", str(tmp_path / "cache")])
        out = capsys.readouterr().out
        assert rc == 3, out
        assert "[SKIP]" in out
        assert "unverified" in out

    def test_a_truncated_reference_download_is_rejected(self, tmp_path, monkeypatch):
        """A short read of the reference makes the local file look like it has extra
        genes -- the opposite of the truth -- so it must not be cached and compared."""
        class Resp:
            headers = {"Content-Length": "999999"}
            def read(self, n):
                if getattr(self, "done", False):
                    return b""
                self.done = True
                return b"gene\ttranscript\n"
            def __enter__(self):
                return self
            def __exit__(self, *a):
                return False
        monkeypatch.setattr(vl.urllib.request, "urlopen", lambda *a, **k: Resp())
        with pytest.raises(IOError, match="truncated"):
            vl.fetch_release("4.1", str(tmp_path))
        assert not os.path.exists(
            os.path.join(str(tmp_path), "gnomad.v4.1.constraint_metrics.tsv"))


class TestTheAutosomeOnlyGapIsDeclared:
    def test_the_expected_contigs_are_exactly_the_autosomes(self):
        """The v4.1 release this table is pinned to publishes no chrX or chrY constraint.
        O8's 'not evaluated' note for X-linked genes, and the LOEUF description in
        Config.yaml, both rest on that.

        Two different wrong things are guarded against here. Filling the gap from v2.1.1
        would score an X-linked gene against a calibration that reaches the 3.09
        threshold about half as often -- a sex-dependent bias, invisibly. Swapping in
        gnomAD v4.1.1 or v4.1.2, which do carry the sex chromosomes, is the right route
        but is not a file swap: on the 17,460 genes shared with v4.1, 750 more reach
        3.09 and 21 fewer do, none of them on a sex chromosome. Either way a file that
        arrives carrying chrX must stop the run and be looked at, not pass quietly."""
        assert vl.EXPECTED_CONTIGS == tuple("chr%d" % i for i in range(1, 23))
        assert "chrX" not in vl.EXPECTED_CONTIGS
        assert "chrY" not in vl.EXPECTED_CONTIGS

    def test_autosome_only_is_stated_not_assumed(self, tmp_path, capsys):
        rows = release_rows(3)
        rel = write_release(tmp_path / "rel.tsv", rows)
        tab = write_table(tmp_path / "t.tsv.gz", rows)
        vl.main([tab, "--release-file", rel])
        out = capsys.readouterr().out
        assert "autosomes only" in out
        assert "no chrX/chrY constraint" in out

    def test_a_sex_chromosome_row_warns_and_says_what_to_revisit(self, tmp_path, capsys):
        rows = release_rows(3)
        rel = write_release(tmp_path / "rel.tsv", rows)
        tab = write_table(tmp_path / "t.tsv.gz", rows,
                          contigs=["chr1", "chrX", "chr3"])
        vl.main([tab, "--release-file", rel])
        out = capsys.readouterr().out
        assert "[WARN]" in out
        assert "chrX" in out
        assert "O8" in out and "Config.yaml" in out

    def test_contigs_are_reported_in_genomic_not_lexical_order(self, tmp_path, capsys):
        """chr10 sorting before chr2 would make a complete list look wrong at a glance."""
        rows = release_rows(3)
        rel = write_release(tmp_path / "rel.tsv", rows)
        tab = write_table(tmp_path / "t.tsv.gz", rows, contigs=["chr10", "chr2", "chr1"])
        vl.main([tab, "--release-file", rel])
        out = [l for l in capsys.readouterr().out.splitlines() if "contigs:" in l][0]
        assert out.split("contigs:")[1].split() == ["chr1", "chr2", "chr10"]


class TestTheReleaseURLsAreTheCitableOnes:
    def test_the_default_release_is_the_one_the_config_declares(self):
        assert "4.1" in vl.RELEASES
        assert vl.RELEASES["4.1"].endswith("gnomad.v4.1.constraint_metrics.tsv")
        assert vl.RELEASES["4.1"].startswith(
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/")

    def test_no_release_url_is_an_unversioned_or_private_link(self):
        """The defect being fixed: the table came from a Google Drive link with no
        version in it and no way to tell what it served."""
        for rel, url in vl.RELEASES.items():
            assert "drive.google.com" not in url
            assert rel in url, (rel, url)
