"""
Which sample column is the tumour, and what the matched normal contributes.

SMART reports one sample, and every value it derives from a sample column — GT, AD, DP, VAF —
describes the tumour. It used to take whichever column came first, without checking that there
was only one. For a Strelka or SomaticSniper VCF, where the normal precedes the tumour, that
reported the normal: a variant at 55% in the tumour came out at **VAF 0.02 with GT 0/0**,
which reads as a subclonal call rather than as the wrong sample. Nothing in the output said
two samples had been seen.

The three paired callers in `tests/verification2` disagree about the order — Mutect2 3.8 writes
TUMOR first, Strelka and SomaticSniper write NORMAL first — and none of them declares
`##tumor_sample`. So the order cannot be assumed and the name is what is available.

These tests cover the four outcomes and the numbers that come out of each.
"""
import importlib.util
import os
import re
import subprocess
import sys

import pandas as pd
import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
SCRIPTS = os.path.join(ROOT, "scripts")

_v = importlib.util.spec_from_file_location("vcf2table", os.path.join(SCRIPTS, "vcf2table.py"))
v2t = importlib.util.module_from_spec(_v)
sys.modules["vcf2table"] = v2t
_v.loader.exec_module(v2t)

_p = importlib.util.spec_from_file_location("post", os.path.join(SCRIPTS, "post_analysis.py"))
post = importlib.util.module_from_spec(_p)
sys.modules["post"] = post
_p.loader.exec_module(post)

CSQ_HDR = ('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. '
           'Format: Allele|Consequence|SYMBOL|Feature|HGVSc|HGVSp">')
CSQ = "CSQ=T|missense_variant|BRAF|NM_004333.6|c.1799T>A|p.Val600Glu"


def write_vcf(tmp_path, samples, data, meta=()):
    p = tmp_path / "in.vcf"
    lines = ["##fileformat=VCFv4.2", *meta, CSQ_HDR,
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
             "chr7\t140753336\t.\tA\tT\t.\tPASS\t" + CSQ + "\tGT:AD\t" + "\t".join(data)]
    p.write_text("\n".join(lines) + "\n")
    tx = tmp_path / "tx.txt"
    tx.write_text("NM_004333.6\n")
    return str(p), str(tx)


def run(tmp_path, samples, data, meta=(), extra=()):
    """vcf2table then the post_analysis steps that read a sample column."""
    vcf, tx = write_vcf(tmp_path, samples, data, meta)
    out = tmp_path / "out.tsv"
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "vcf2table.py"),
                        "--vcf", vcf, "--transcripts", tx, "--output", str(out), *extra],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr)
    df = pd.read_csv(out, sep="\t", dtype=str)
    df = post.standardize_format_column(df)
    df = post.normalize_format_fields(df)
    df = post.add_vaf_column(df)
    return df.iloc[0], r.stdout


# A tumour at 55% and a clean normal, so a mix-up is unmistakable in the numbers.
TUMOUR = "0/1:45,55"
NORMAL = "0/0:98,2"


class TestWhichColumnIsTheTumour:
    def test_a_single_sample_column_is_the_tumour(self, tmp_path):
        """The laboratory's routine DRAGEN input, which must be unaffected by all of this."""
        row, log = run(tmp_path, ["D00_Thyroid"], [TUMOUR])
        assert row["Tumor_Sample_Barcode"] == "D00_Thyroid"
        assert "the only sample column" in log
        assert float(row["VAF"]) == pytest.approx(0.55)

    def test_the_tumour_is_found_when_the_normal_comes_first(self, tmp_path):
        """Strelka and SomaticSniper order. This is the defect these tests exist for."""
        row, log = run(tmp_path, ["NORMAL", "TUMOR"], [NORMAL, TUMOUR])
        assert row["Tumor_Sample_Barcode"] == "TUMOR"
        assert row["GT"] == "0/1", "the normal's genotype was reported"
        assert float(row["VAF"]) == pytest.approx(0.55), "the normal's VAF was reported"

    def test_the_tumour_is_found_when_it_comes_first(self, tmp_path):
        """Mutect2 3.8 order. Both orders must give the same answer."""
        row, _ = run(tmp_path, ["TUMOR", "NORMAL"], [TUMOUR, NORMAL])
        assert row["Tumor_Sample_Barcode"] == "TUMOR"
        assert float(row["VAF"]) == pytest.approx(0.55)

    def test_the_order_does_not_change_any_value(self, tmp_path):
        """Stated as its own test: column order must be irrelevant to every derived number."""
        a, _ = run(tmp_path, ["NORMAL", "TUMOR"], [NORMAL, TUMOUR])
        b, _ = run(tmp_path, ["TUMOR", "NORMAL"], [TUMOUR, NORMAL])
        for c in ("GT", "AD", "VAF", "t_ref_count", "t_alt_count", "t_depth",
                  "n_alt_count", "n_depth", "VAF_NORMAL"):
            assert str(a[c]) == str(b[c]), f"{c} depends on column order"

    def test_a_declared_tumour_sample_wins_over_the_name(self, tmp_path):
        """Mutect2 4.x writes ##tumor_sample. An explicit statement beats any guess."""
        row, log = run(tmp_path, ["SampA", "SampB"], [NORMAL, TUMOUR],
                       meta=["##tumor_sample=SampB", "##normal_sample=SampA"])
        assert row["Tumor_Sample_Barcode"] == "SampB"
        assert row["Matched_Norm_Sample_Barcode"] == "SampA"
        assert "##tumor_sample" in log
        assert float(row["VAF"]) == pytest.approx(0.55)

    def test_unrecognisable_names_are_refused_not_guessed(self, tmp_path):
        """A wrong sample is a plausible low allele fraction, so refusing is the safer failure."""
        with pytest.raises(RuntimeError) as e:
            run(tmp_path, ["SampA", "SampB"], [TUMOUR, NORMAL])
        msg = str(e.value)
        assert "cannot tell which sample column is the tumour" in msg
        assert "'SampA'" in msg and "'SampB'" in msg, "the message must name the columns"
        assert "--tumour-column" in msg, "the message must say how to resolve it"

    def test_the_operator_can_name_the_column(self, tmp_path):
        row, log = run(tmp_path, ["SampA", "SampB"], [TUMOUR, NORMAL],
                       extra=["--tumour-column", "SampA", "--normal-column", "SampB"])
        assert row["Tumor_Sample_Barcode"] == "SampA"
        assert "--tumour-column" in log
        assert float(row["VAF"]) == pytest.approx(0.55)

    def test_naming_a_column_that_does_not_exist_is_refused(self, tmp_path):
        with pytest.raises(RuntimeError) as e:
            run(tmp_path, ["SampA"], [TUMOUR], extra=["--tumour-column", "Nope"])
        assert "is not a column of this VCF" in str(e.value)

    def test_the_identification_ignores_derived_columns(self, tmp_path):
        """It must read the VCF's own columns, not the frame's.

        An earlier version passed the DataFrame's columns, which by that point also hold
        NM_Transcript, every CSQ field and the derived TRANSCRIPT_* and END columns — so
        anything outside the nine fixed VCF columns looked like a sample and even a
        single-sample VCF became unidentifiable.
        """
        row, log = run(tmp_path, ["D00_Thyroid"], [TUMOUR])
        assert row["Tumor_Sample_Barcode"] == "D00_Thyroid"
        for derived in ("NM_Transcript", "Consequence", "SYMBOL", "END"):
            assert derived not in log


class TestTheMatchedNormalColumns:
    def test_the_standard_maf_counts_are_emitted(self, tmp_path):
        """t_* and n_* are the MAF format's own answer to this, so they are not invented here.
        None of them existed before, not even in the tumour-only case."""
        row, _ = run(tmp_path, ["NORMAL", "TUMOR"], [NORMAL, TUMOUR])
        assert int(row["t_ref_count"]) == 45
        assert int(row["t_alt_count"]) == 55
        assert int(row["t_depth"]) == 100
        assert int(row["n_ref_count"]) == 98
        assert int(row["n_alt_count"]) == 2
        assert int(row["n_depth"]) == 100

    def test_vaf_normal_distinguishes_somatic_from_germline(self, tmp_path):
        """The reason to carry the normal at all, rather than only refusing a paired VCF.

        A fraction near 0 in the normal says somatic; near 0.5 or 1.0 says germline. SMART
        could not show this before, because the normal's data was either discarded or
        reported as the tumour's.
        """
        somatic, _ = run(tmp_path, ["NORMAL", "TUMOR"], ["0/0:98,2", "0/1:45,55"])
        assert float(somatic["VAF_NORMAL"]) == pytest.approx(0.02)

        germline, _ = run(tmp_path, ["NORMAL", "TUMOR"], ["0/1:50,50", "0/1:45,55"])
        assert float(germline["VAF_NORMAL"]) == pytest.approx(0.50)

    def test_the_normal_fields_are_prefixed_so_they_cannot_be_mistaken(self, tmp_path):
        """N_AD is the normal's allelic depth and nothing else can claim that name.

        The tumour's FORMAT keys already contend with INFO and with VEP's CSQ for column
        names; expanding a second sample into the same names would make which sample a cell
        describes depend on which source claimed the name first.
        """
        row, _ = run(tmp_path, ["NORMAL", "TUMOR"], [NORMAL, TUMOUR])
        assert row["N_AD"] == "98,2"
        assert row["N_GT"] == "0/0"
        assert row["AD"] == "45,55", "the unprefixed column must be the tumour's"

    def test_a_tumour_only_vcf_leaves_the_normal_columns_empty(self, tmp_path):
        """Empty, not absent and not zero: there was no normal, which is different from a
        normal with no reads."""
        row, _ = run(tmp_path, ["D00_Thyroid"], [TUMOUR])
        assert row["Matched_Norm_Sample_Barcode"] in ("", None) or pd.isna(
            row["Matched_Norm_Sample_Barcode"])
        for c in ("n_ref_count", "n_alt_count", "n_depth", "VAF_NORMAL"):
            assert pd.isna(row[c]) or str(row[c]) == "", f"{c} should be empty, got {row[c]!r}"
        assert int(row["t_depth"]) == 100, "the tumour counts must still be there"


class TestTheFlagIsReachableFromDocker:
    """`vcf2table` refuses an ambiguous tumour/normal VCF and names `--tumour-column` in
    the refusal. `entrypoint.sh` did not parse that flag and invoked vcf2table with only
    `--vcf --transcripts --output --debug`, so a container user was pointed at a remedy
    they had no way to apply -- and the CHANGELOG offered it as the escape hatch.
    """

    ENTRYPOINT = os.path.join(os.path.dirname(__file__), "..", "..", "entrypoint.sh")

    def text(self):
        return open(self.ENTRYPOINT, encoding="utf-8").read()

    def test_the_flags_are_parsed(self):
        s = self.text()
        assert "--tumour-column|--tumor-column)" in s, "not in the argument loop"
        assert "--normal-column)" in s

    def test_they_reach_vcf2table(self):
        s = self.text()
        i = s.index('vcf2table.py" \\')
        call = s[i:i + 500]
        assert "--tumour-column" in call, "parsed but never passed on"
        assert "--normal-column" in call

    def test_they_are_omitted_when_unset_rather_than_passed_empty(self):
        """`--tumour-column ''` is not the same request as no flag: vcf2table would take
        the empty string as a column name and refuse for the wrong reason."""
        s = self.text()
        i = s.index('vcf2table.py" \\')
        call = s[i:i + 500]
        assert "${TUMOUR_COLUMN:+" in call, "must use the :+ form so unset means absent"
        assert "${NORMAL_COLUMN:+" in call

    def test_they_are_documented_in_the_usage_text(self):
        s = self.text()
        usage = s[s.index("Options:"):s.index("Options:") + 2500]
        assert "--tumour-column" in usage, (
            "a flag a user is told to use must appear in --help")

    def test_the_refusal_still_names_the_flag(self):
        """The two halves have to agree: the message names the flag, the entrypoint
        accepts it."""
        v2t = open(os.path.join(os.path.dirname(__file__), "..", "..",
                                "scripts", "vcf2table.py"), encoding="utf-8").read()
        assert "--tumour-column" in v2t


class TestTheImageCarriesEveryClassifierReferenceList:
    """Two of the four SVIG-UK tables were absent from the Docker image.

    The classifier resolves a default path for each, under `reference_lists/`. The image
    copied only `svig_uk_canonical_variants.tsv` and
    `svig_uk_vcep_population_thresholds.tsv`, so in a container:

      - O10/B6 could not load `svig_uk_functional_assays.tsv` and fell back to the OncoKB
        mutation-effect proxy, capped at supporting -- 1 point where the shipped assay
        table gives 2;
      - the TP53 exon 4/6/9 synonymous provision in `svig_uk_gene_specific_splice.tsv`
        never applied.

    Neither emitted a note, so the same variant could be classified differently inside
    and outside the container with nothing in the output to explain it. These are the
    guideline's own content, already in the repository, with no licence restriction --
    there was no reason for them to be left out.
    """

    ROOT = os.path.join(os.path.dirname(__file__), "..", "..")

    def dockerfile(self):
        return open(os.path.join(self.ROOT, "Dockerfile"), encoding="utf-8").read()

    def defaulted_lists(self):
        """Reference lists the classifier names a default path for."""
        src = open(os.path.join(self.ROOT, "scripts", "svig_uk_classifier.py"),
                   encoding="utf-8").read()
        return sorted(set(re.findall(r"reference_lists/(\S+?\.tsv)", src)))

    def test_the_classifier_names_defaults_at_all(self):
        """The premise: if it stops defaulting to reference_lists/, this test is moot."""
        assert self.defaulted_lists(), "no reference_lists/ defaults found"

    def test_every_default_is_copied_into_the_image(self):
        copied = self.dockerfile()
        missing = [f for f in self.defaulted_lists()
                   if f"reference_lists/{f}" not in copied]
        assert not missing, (
            f"the classifier defaults to {missing} but the Dockerfile does not copy "
            f"them, so a container run silently takes a different route")

    def test_the_scores_file_travels_with_its_parameters(self):
        """`svig_uk_functional_assays.tsv` is useless without the scores file -- the
        loader takes both -- so copying one and not the other would reintroduce the
        silent fallback."""
        copied = self.dockerfile()
        assert ("svig_uk_functional_assays.tsv" in copied) == \
               ("svig_uk_functional_assay_scores.tsv" in copied)

    def test_every_copied_list_exists_in_the_repository(self):
        copied = re.findall(r"reference_lists/(\S+?\.tsv)", self.dockerfile())
        assert copied
        for f in copied:
            assert os.path.exists(os.path.join(self.ROOT, "reference_lists", f)), \
                f"Dockerfile copies reference_lists/{f}, which does not exist"
