"""
Unit tests for scripts/verification_report.py.

A full run leaves its evidence in half a dozen files. The report composes them into one
document: every variant with the annotation it received, every check that passed and against
which authority, every check that did not and why, and the versions of everything involved.

The property that matters most is negative: a section whose input is missing must say so
rather than be omitted. A report that quietly drops what it could not find is worse than no
report, because it reads as complete.
"""
import importlib.util
import io
import json
import os
import sys

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "..", "scripts")
_spec = importlib.util.spec_from_file_location(
    "verification_report", os.path.join(SCRIPTS, "verification_report.py"))
vr = importlib.util.module_from_spec(_spec)
sys.modules["verification_report"] = vr
_spec.loader.exec_module(vr)

REAL_V1 = os.path.join(os.path.dirname(__file__), "..", "verification1")
REAL_MAF = os.path.join(REAL_V1, "output", "output", "Final_result_tier1.maf")


class Args:
    def __init__(self, **kw):
        self.verification1 = kw.get("verification1")
        self.vcf = kw.get("vcf")
        self.verification3 = kw.get("verification3")
        self.allowlist = kw.get("allowlist")
        self.out = None


def make_run(tmp_path, maf_rows=None, api=None, local=None, provenance=None, meta=None):
    """A minimal verification1 directory, with only the parts a test needs."""
    v1 = tmp_path / "verification1"
    (v1 / "output" / "output").mkdir(parents=True)
    (v1 / "output" / "logs").mkdir(parents=True)
    hdr = ["Hugo_Symbol", "ID", "HGVSc", "HGVSp_Short", "Variant_Classification",
           "Consequence", "IMPACT", "NM_Transcript", "TRANSCRIPT_MATCH"]
    rows = maf_rows if maf_rows is not None else [
        ["BRAF", "BRAF_V600E", "c.1799T>A", "p.V600E", "Missense_Mutation",
         "missense_variant", "MODERATE", "NM_004333", "match"]]
    with open(v1 / "output" / "output" / "Final_result_tier1.maf", "w") as fh:
        fh.write("#SMART_VERSION 2.1.0\n")
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    if api is not None:
        with open(v1 / "results.tsv", "w") as fh:
            fh.write("Module\tSample\tVariant_ID\tGene\tMAF_Field\tMAF_Value\tAPI_Value\tResult\n")
            for r in api:
                fh.write("\t".join(r) + "\n")
    if local is not None:
        with open(v1 / "local_source_results.tsv", "w") as fh:
            fh.write("Module\tSample\tVariant_ID\tGene\tMAF_Field\tMAF_Value\tSource_Value\tResult\tNote\n")
            for r in local:
                fh.write("\t".join(r) + "\n")
    if provenance:
        (v1 / "output" / "logs" / "run_provenance.txt").write_text(provenance)
    if meta:
        (v1 / "output" / "logs" / "s_metadata.json").write_text(json.dumps(meta))
    return str(v1)


def render(**kw):
    """The Markdown rendering, which is what most assertions below read."""
    blocks, _ = vr.collect(Args(**kw))
    return vr.render_markdown(blocks)


class TestMissingInputsAreStatedNotOmitted:
    def test_absent_provenance_is_reported(self, tmp_path):
        text = render(verification1=make_run(tmp_path))
        assert "run_provenance.txt" in text and "absent" in text, (
            "a missing provenance file was omitted rather than reported, so the report "
            "reads as though the versions were simply not relevant")

    def test_absent_api_results_are_reported(self, tmp_path):
        text = render(verification1=make_run(tmp_path))
        assert "Not available" in text

    def test_absent_allowlist_is_reported(self, tmp_path):
        text = render(verification1=make_run(tmp_path))
        assert "Columns that carry no value" in text
        assert "Not available" in text

    def test_an_absent_maf_is_a_failure_not_an_empty_report(self, tmp_path):
        v1 = tmp_path / "empty"
        (v1 / "output" / "output").mkdir(parents=True)
        blocks, ok = vr.collect(Args(verification1=str(v1)))
        assert ok is False, "a run that produced no MAF should not look like a clean report"
        assert "did not complete" in vr.render_markdown(blocks)


class TestTheCohortSection:
    def test_every_variant_appears(self, tmp_path):
        rows = [["BRAF", "BRAF_V600E", "c.1799T>A", "p.V600E", "Missense_Mutation",
                 "missense_variant", "MODERATE", "NM_004333", "match"],
                ["TERT", "TERT_C228T", "", "", "5'Flank",
                 "upstream_gene_variant", "MODIFIER", "NM_198253", "match"]]
        text = render(verification1=make_run(tmp_path, maf_rows=rows))
        assert "BRAF_V600E" in text and "TERT_C228T" in text
        assert "2 variants" in text, "the annotated-variants count is missing"

    def test_sentinels_are_not_printed_as_annotations(self, tmp_path):
        """A MAF read as text carries 'nan' and '.'; printing either would present a file
        format's artefact as if it were an annotation."""
        rows = [["CDK4", "MantaDUP:CDK4_AMP", "nan", ".", "5'UTR",
                 "coding_sequence_variant", "MODIFIER", "NM_000075", "match"]]
        text = render(verification1=make_run(tmp_path, maf_rows=rows))
        assert "| nan |" not in text and "| . |" not in text
        assert "—" in text

    def test_the_transcript_guarantee_is_counted(self, tmp_path):
        rows = [["BRAF", "v1", "c.1T>A", "p.X1Y", "Missense_Mutation", "missense_variant",
                 "MODERATE", "NM_004333", "match"],
                ["XYZ", "v2", "c.2T>A", "p.X2Y", "Missense_Mutation", "missense_variant",
                 "MODERATE", "", "fallback:NM_999999.1"]]
        text = render(verification1=make_run(tmp_path, maf_rows=rows))
        assert "1 of 2 rows were annotated on a transcript that was supplied" in text
        assert "fallback:NM_999999.1" in text, (
            "a row annotated on a transcript nobody asked for must be named")


class TestTheVerificationSections:
    API = [["OncoKB", "s", "v1", "BRAF", "ONCOGENIC", "Oncogenic", "Oncogenic", "PASS"],
           ["VEP", "s", "v1", "BRAF", "Consequence", "a", "b", "MISMATCH"],
           ["VEP", "s", "v2", "TERT", "HGVSc", "", "", "SKIP: no HGVSc/Feature"],
           ["CIViC", "s", "v1", "BRAF", "PRESENCE", "", "id=1", "COVERAGE_GAP"]]

    def test_the_tally_is_per_module_with_a_total(self, tmp_path):
        text = render(verification1=make_run(tmp_path, api=self.API))
        assert "| OncoKB | 1" in text
        assert "**total**" in text

    def test_a_mismatch_is_named_with_both_values(self, tmp_path):
        text = render(verification1=make_run(tmp_path, api=self.API))
        assert "1 mismatch(es)" in text
        assert "Consequence" in text

    def test_no_mismatches_is_stated_positively(self, tmp_path):
        api = [["OncoKB", "s", "v1", "BRAF", "ONCOGENIC", "x", "x", "PASS"]]
        text = render(verification1=make_run(tmp_path, api=api))
        assert "No mismatches" in text

    def test_skip_reasons_are_given_not_just_counted(self, tmp_path):
        """A skip is a comparison that cannot be made; without the reason it reads as a gap."""
        text = render(verification1=make_run(tmp_path, api=self.API))
        assert "no HGVSc/Feature" in text

    def test_errors_are_separated_from_disagreements(self, tmp_path):
        api = [["VEP", "s", "v1", "BRAF", "HGVSc", "", "", "ERROR: 500 Server Error"]]
        text = render(verification1=make_run(tmp_path, api=api))
        assert "did not answer" in text
        assert "mismatch" not in text.split("## Checked against each")[1].split("##")[0].lower() \
            or "No mismatches" in text

    def test_both_authorities_get_their_own_section(self, tmp_path):
        local = [["dbNSFP", "s", "v1", "BRAF", "REVEL", "0.9", "0.9", "PASS", ""]]
        text = render(verification1=make_run(tmp_path, api=self.API, local=local))
        assert "Checked against each source's own API" in text
        assert "Checked against the local source files" in text
        assert "| dbNSFP | 1" in text


class TestProvenance:
    def test_the_provenance_block_is_quoted_verbatim(self, tmp_path):
        prov = "Tool versions:\n  VEP:  114.2\nReference versions:\n  ClinVar: clinvar.vcf.gz"
        text = render(verification1=make_run(tmp_path, provenance=prov))
        assert "VEP:  114.2" in text and "clinvar.vcf.gz" in text

    def test_field_sources_are_grouped_with_their_versions(self, tmp_path):
        meta = {"fields": {"REVEL": {"source": "VEP", "version": "v114.0"},
                           "ONCOGENIC": {"source": "OncoKB", "version": "v7.6"},
                           "CADD_PHRED": {"source": "VEP", "version": "v114.0"}}}
        text = render(verification1=make_run(tmp_path, meta=meta))
        assert "3 annotated fields" in text
        assert "| VEP | v114.0 | 2 |" in text

    def test_the_smart_version_is_taken_from_the_maf(self, tmp_path):
        text = render(verification1=make_run(tmp_path))
        assert "SMART version 2.1.0" in text, (
            "the version must come from the output itself, not from the VERSION file, so a "
            "report cannot claim a version the run did not produce")


class TestAgainstTheRealRun:
    @pytest.mark.skipif(not os.path.exists(REAL_MAF),
                        reason="verification1 has not been run")
    def test_a_report_can_be_built_from_the_real_artefacts(self):
        blocks, ok = vr.collect(Args(
            verification1=REAL_V1,
            verification3=os.path.join(os.path.dirname(__file__), "..",
                                       "verification3", "results.tsv"),
            allowlist=os.path.join(REAL_V1, "empty_columns_allowed.tsv")))
        text = vr.render_markdown(blocks)
        assert ok is True
        for section in ("## What was run, and what was verified",
                        "## The variants put into the test, as submitted",
                        "## What each variant was annotated as",
                        "## Checked against each source's own API",
                        "## Checked against the local source files",
                        "## Transcript prioritisation",
                        "## Columns that carry no value",
                        "## Where each annotation came from",
                        "## Checks that could not be made, and why",
                        "## Appendix A", "## Appendix B"):
            assert section in text, f"{section} missing from the report"
        assert "variants." in text


class TestThePdfRendering:
    """The PDF is the publishable artefact, so it is rendered from the same blocks.

    Converting the Markdown to PDF instead would mean re-parsing prose to recover structure
    that had just been discarded, and the two outputs could then disagree about what the run
    found. These tests assert the single-source property as much as the file itself.
    """

    def test_a_pdf_is_produced_and_is_a_pdf(self, tmp_path):
        blocks, _ = vr.collect(Args(verification1=make_run(tmp_path)))
        out = tmp_path / "r.pdf"
        assert vr.render_pdf(blocks, str(out)) is True
        data = out.read_bytes()
        assert data.startswith(b"%PDF"), "not a PDF"
        assert b"/Type /Page" in data
        assert len(data) > 1000, "a PDF this small cannot contain the report"

    def test_the_cohort_table_reaches_the_pdf(self, tmp_path):
        """A table silently dropped by the renderer is the failure that matters here."""
        rows = [["BRAF", "BRAF_V600E", "c.1799T>A", "p.V600E", "Missense_Mutation",
                 "missense_variant", "MODERATE", "NM_004333", "match"]]
        blocks, _ = vr.collect(Args(verification1=make_run(tmp_path, maf_rows=rows)))
        n_tables = sum(1 for b in blocks if b[0] == "table")
        assert n_tables >= 1
        out = tmp_path / "r.pdf"
        assert vr.render_pdf(blocks, str(out)) is True
        # reportlab compresses content streams, so assert on the structure rather than by
        # searching for the text: a table becomes many small drawing operations.
        assert out.stat().st_size > 1500

    def test_both_renderings_come_from_one_block_list(self, tmp_path):
        """If a section exists in one output it exists in the other, by construction."""
        blocks, _ = vr.collect(Args(verification1=make_run(tmp_path)))
        md = vr.render_markdown(blocks)
        headings = [b[1] for b in blocks if b[0] == "h2"]
        assert headings, "no sections were produced"
        for h in headings:
            assert f"## {h}" in md
        out = tmp_path / "r.pdf"
        assert vr.render_pdf(blocks, str(out)) is True

    def test_a_missing_reportlab_is_reported_not_fatal(self, tmp_path, monkeypatch):
        """The Markdown is already written by then, and an optional dependency should not
        fail a verification run that otherwise succeeded."""
        import builtins
        real = builtins.__import__

        def no_reportlab(name, *a, **kw):
            if name.startswith("reportlab"):
                raise ImportError("simulated")
            return real(name, *a, **kw)

        monkeypatch.setattr(builtins, "__import__", no_reportlab)
        blocks, _ = vr.collect(Args(verification1=make_run(tmp_path)))
        assert vr.render_pdf(blocks, str(tmp_path / "r.pdf")) is False

    def test_a_pipe_in_a_value_cannot_break_the_markdown_table(self):
        """'|' inside a cell would end the column and shift every later one."""
        md = vr.render_markdown([("table", ["a", "b"], [["x|y", "z"]])])
        assert md.count("|") == 3 + 3 + 3 + 1, md   # header, rule, row, plus the escaped one
        assert "x\\|y" in md

    @pytest.mark.skipif(not os.path.exists(REAL_MAF),
                        reason="verification1 has not been run")
    def test_the_real_report_renders_to_several_pages(self, tmp_path):
        blocks, ok = vr.collect(Args(
            verification1=REAL_V1,
            verification3=os.path.join(os.path.dirname(__file__), "..",
                                       "verification3", "results.tsv"),
            allowlist=os.path.join(REAL_V1, "empty_columns_allowed.tsv")))
        assert ok
        out = tmp_path / "real.pdf"
        assert vr.render_pdf(blocks, str(out)) is True
        data = out.read_bytes()
        pages = data.count(b"/Type /Page") - data.count(b"/Type /Pages")
        assert pages >= 3, f"the real report rendered to only {pages} page(s)"


class TestTheSubmittedVariantsPage:
    """The input VCF, shown as submitted.

    What was *put into* the test is a different question from what came out, and a reader
    checking "was my variant tested, at the coordinate I use" wants the five columns that
    identify a variant rather than our column conventions.
    """

    @staticmethod
    def _vcf(tmp_path, records):
        p = tmp_path / "in.vcf"
        p.write_text("##fileformat=VCFv4.2\n"
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                     "".join("\t".join(r) + "\t.\tPASS\t.\n" for r in records))
        return str(p)

    def test_the_records_are_shown_verbatim(self, tmp_path):
        vcf = self._vcf(tmp_path, [["chr7", "140753336", "BRAF_V600E", "A", "T"]])
        text = render(verification1=make_run(tmp_path), vcf=vcf)
        assert "#CHROM" in text
        assert "chr7\t140753336\tBRAF_V600E\tA\tT" in text

    def test_a_symbolic_alt_survives(self, tmp_path):
        """A CNA has no sequence to state; the record must still be shown."""
        vcf = self._vcf(tmp_path, [["chr7", "116672196", "MantaDUP:MET_AMP", "N", "<DUP>"]])
        text = render(verification1=make_run(tmp_path), vcf=vcf)
        assert "<DUP>" in text

    def test_records_are_ordered_by_position(self, tmp_path):
        vcf = self._vcf(tmp_path, [["chr12", "1", "b", "A", "T"],
                                   ["chr2", "9", "a", "C", "G"]])
        rows = vr.read_vcf(vcf)
        assert [r[0] for r in rows] == ["chr2", "chr12"], "contigs sorted as strings"

    def test_an_absent_vcf_is_stated(self, tmp_path):
        text = render(verification1=make_run(tmp_path))
        assert "was not given to this report" in text

    def test_a_bgzipped_vcf_is_read(self, tmp_path):
        import gzip as _gz
        p = tmp_path / "in.vcf.gz"
        with _gz.open(p, "wt") as fh:
            fh.write("#CHROM\tPOS\tID\tREF\tALT\n"
                     "chr1\t100\tv\tA\tG\t.\tPASS\t.\n")
        assert vr.read_vcf(str(p)) == [["chr1", "100", "v", "A", "G"]]


class TestTheErrorsSection:
    def test_errors_are_grouped_by_kind_not_by_url(self, tmp_path):
        """Each row carries a different query string; grouping on it would turn one outage
        into as many distinct problems as there were checks."""
        api = [["VEP", "s", f"v{i}", "G", "HGVSc", "", "",
                f"ERROR: 500 Server Error for url: https://rest.ensembl.org/vep/x{i}"]
               for i in range(5)]
        text = render(verification1=make_run(tmp_path, api=api))
        assert "5 check(s) did not complete" in text
        assert "| 5 |" in text, "five instances of one failure were not collapsed"

    def test_an_error_is_not_presented_as_a_disagreement(self, tmp_path):
        api = [["VEP", "s", "v1", "G", "HGVSc", "", "", "ERROR: timed out"]]
        text = render(verification1=make_run(tmp_path, api=api))
        assert "None of these is a disagreement" in text

    def test_the_affected_modules_are_named(self, tmp_path):
        api = [["VEP", "s", "v1", "G", "HGVSc", "", "", "ERROR: timed out"]]
        local = [["gnomAD", "s", "v1", "G", "AF", "", "", "ERROR: timed out", ""]]
        text = render(verification1=make_run(tmp_path, api=api, local=local))
        assert "**VEP**" in text and "**gnomAD**" in text

    def test_no_errors_is_stated_rather_than_omitted(self, tmp_path):
        api = [["OncoKB", "s", "v1", "G", "ONCOGENIC", "x", "x", "PASS"]]
        text = render(verification1=make_run(tmp_path, api=api))
        assert "None: every check completed." in text


class TestTheAppendix:
    API = [["OncoKB", "s", "NRAS_Q61R", "NRAS", "ONCOKB_VUS", "False", "False", "PASS"],
           ["OncoKB", "s", "NRAS_Q61R", "NRAS", "ONCOGENIC", "Oncogenic", "Oncogenic", "PASS"],
           ["VEP", "s", "NRAS_Q61R", "NRAS", "Consequence", "a", "b", "MISMATCH"]]

    def test_every_comparison_appears(self, tmp_path):
        text = render(verification1=make_run(tmp_path, api=self.API))
        assert "Appendix A" in text
        for field in ("ONCOKB_VUS", "ONCOGENIC", "Consequence"):
            assert field in text.split("Appendix A")[1]

    def test_a_pass_and_a_failure_are_marked_differently(self, tmp_path):
        tail = render(verification1=make_run(tmp_path, api=self.API)).split("Appendix A")[1]
        assert "\u2713 ONCOKB_VUS" in tail
        assert "\u2717 Consequence" in tail

    def test_both_values_are_shown_so_a_reader_can_judge(self, tmp_path):
        tail = render(verification1=make_run(tmp_path, api=self.API)).split("Appendix A")[1]
        assert "MAF='False'" in tail and "API='False'" in tail

    def test_comparisons_are_grouped_under_their_variant_and_module(self, tmp_path):
        tail = render(verification1=make_run(tmp_path, api=self.API)).split("Appendix A")[1]
        assert "[OncoKB]  NRAS NRAS_Q61R" in tail
        assert "[VEP]  NRAS NRAS_Q61R" in tail

    def test_the_local_source_appendix_names_its_authority(self, tmp_path):
        local = [["dbNSFP", "s", "v1", "G", "REVEL", "0.9", "0.9", "PASS", ""]]
        text = render(verification1=make_run(tmp_path, local=local))
        assert "Appendix B" in text
        assert "source='0.9'" in text.split("Appendix B")[1]
