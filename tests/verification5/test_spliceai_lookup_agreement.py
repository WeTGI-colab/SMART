"""
Check SMART's splice annotation against an independent implementation of the same model.

The Broad's SpliceAI-lookup runs SpliceAI on its own infrastructure, so it is the one thing
that can say whether our numbers are right rather than merely self-consistent. Its answers for
a small set of **public** variants are frozen in
`data/spliceai_lookup_frozen.json`; nothing patient-derived is listed, and the fixture records
when it was queried and against which service.

Frozen rather than queried on each run, for three reasons. The service is rate limited and asks
that batch work be done locally. It returned HTTP 500 on four of 93 rows during this work, so a
live test would fail for reasons that are not ours. And a frozen answer is a dated record, which
is what makes a later difference meaningful: if these values stop matching, either our annotation
changed or their Gencode release did, and both are worth knowing. Refresh with
`utils/freeze_spliceai_lookup.py` deliberately, not automatically.

Note that the service runs a **modified** SpliceAI (github.com/bw2/SpliceAI) with MNV support and
raw REF/ALT scores that the released package does not expose, so agreement is expected for
single-nucleotide substitutions and not for multi-nucleotide ones.
"""
import json
import os
import shutil
import subprocess

import pytest

HERE = os.path.dirname(__file__)
FROZEN = os.path.join(HERE, "data", "spliceai_lookup_frozen.json")
REFS = "/Users/monkiky/Documents/external/refs"
PRECOMPUTED = os.path.join(REFS, "SpliceAI", "spliceai_scores.raw.snv.hg38.vcf.gz")

DS = ("DS_AG", "DS_AL", "DS_DG", "DS_DL")


def frozen():
    with open(FROZEN) as fh:
        return json.load(fh)


def best(scores, refseq=None):
    """The record for `refseq`, else the one with the highest delta."""
    pool = [s for s in scores if refseq and refseq in (s.get("t_refseq_ids") or [])]
    pool = pool or scores
    def peak(s):
        return max([float(s[k]) for k in DS if s.get(k) is not None] or [-1.0])
    return max(pool, key=peak)


class TestTheFixtureItself:
    """What the frozen answers assert about the published tool, independent of our code."""

    def test_provenance_is_recorded(self):
        p = frozen()["_provenance"]
        for k in ("source", "endpoint", "queried", "note"):
            assert p.get(k), f"the fixture must say {k}"

    def test_every_variant_says_why_it_is_there(self):
        for v, d in frozen()["variants"].items():
            assert d.get("why"), f"{v} has no stated purpose"
            assert d.get("scores"), f"{v} has no scores"

    def test_the_service_was_queried_at_500_bp(self):
        # the window is the point of comparison; at 50 bp there would be nothing to compare
        for v, d in frozen()["variants"].items():
            assert d["distance"] == 500, f"{v} queried at {d['distance']} bp"

    def test_abca3_signal_lies_beyond_the_precomputed_window(self):
        # the finding this fixture exists to lock in: the strongest signal is at 106 bp,
        # outside the +/-50 bp the precomputed tables were generated with
        s = best(frozen()["variants"]["16-2317763-T-A"]["scores"], "NM_001089.3")
        assert abs(int(s["DP_AL"])) > 50
        assert float(s["DS_AL"]) >= 0.2      # and it reaches the threshold we act on

    def test_met_returns_our_reported_transcript_labelled(self):
        # MET is one of the 16 panel transcripts a gene-level score describes wrongly
        scores = frozen()["variants"]["7-116771900-T-C"]["scores"]
        refseqs = {r for s in scores for r in (s.get("t_refseq_ids") or [])}
        assert "NM_001127500.3" in refseqs      # ours
        assert "NM_000245.4" in refseqs         # the MANE Select

    def test_cdkn2a_returns_both_reported_transcripts(self):
        scores = frozen()["variants"]["9-21974826-A-G"]["scores"]
        refseqs = {r for s in scores for r in (s.get("t_refseq_ids") or [])}
        assert {"NM_000077.5", "NM_058195.4"} <= refseqs   # p16INK4a and p14ARF

    def test_raw_ref_and_alt_scores_are_captured(self):
        # the data the issue-21 criterion needs: whether a site was already there, which a
        # delta cannot say. Exposed by the fork, not by the released package.
        s = best(frozen()["variants"]["16-2317763-T-A"]["scores"], "NM_001089.3")
        assert s.get("DS_AL_REF") is not None and s.get("DS_AL_ALT") is not None
        assert float(s["DS_AL_REF"]) > float(s["DS_AL_ALT"]), \
            "an acceptor loss should weaken a site that was already there"


@pytest.mark.skipif(not os.path.exists(PRECOMPUTED),
                    reason="precomputed SpliceAI table not present")
class TestPrecomputedDivergence:
    """Lock in how the precomputed tables differ, so a change of table is noticed."""

    @staticmethod
    def precomputed(chrom, pos, ref, alt):
        r = subprocess.run(["tabix", PRECOMPUTED, f"{chrom}:{pos}-{pos}"],
                           capture_output=True, text=True)
        out = []
        for line in r.stdout.splitlines():
            f = line.split("\t")
            if len(f) > 7 and f[3] == ref and alt in f[4].split(","):
                for rec in f[7].split("SpliceAI=")[-1].split(","):
                    p = rec.split("|")
                    if len(p) == 10:
                        out.append(dict(zip(("ALLELE", "SYMBOL") + DS +
                                            ("DP_AG", "DP_AL", "DP_DG", "DP_DL"), p)))
        return out

    def test_abca3_is_below_the_threshold_in_the_tables_and_above_it_in_the_service(self):
        recs = self.precomputed("16", 2317763, "T", "A")
        assert recs, "the variant should be in the table"
        pre = max(float(r[k]) for r in recs for k in DS)
        srv = max(float(best(frozen()["variants"]["16-2317763-T-A"]["scores"],
                             "NM_001089.3")[k]) for k in DS)
        assert pre < 0.2 <= srv, (
            f"the documented divergence has changed: table {pre}, service {srv}. If the "
            f"tables were replaced this test is the thing that noticed.")

    def test_every_delta_position_in_the_tables_is_within_the_window(self):
        for r in self.precomputed("16", 2317763, "T", "A"):
            for k in ("DP_AG", "DP_AL", "DP_DG", "DP_DL"):
                assert abs(int(r[k])) <= 50


SPLICEAI = shutil.which("spliceai") or os.path.join(
    HERE, "..", "..", "external_tools", "venv_spliceai", "bin", "spliceai")
FASTA = os.path.join(REFS, "SpliceAI-lookup", "references_files", "hg38.fa")
ANNOTATION = os.path.join(REFS, "SpliceAI-lookup", "tso500_grch38.spliceai_annotation.txt")


@pytest.mark.skipif(not (os.path.exists(SPLICEAI) and os.path.exists(FASTA)),
                    reason="no local SpliceAI model or reference FASTA")
class TestLocalComputationAgrees:
    """Our own run of the released model against the frozen service answers.

    One variant, because this loads five Keras models and a 3 GB FASTA. It is the validation
    the local-compute branch rests on, so it is asserted rather than described.
    """

    def test_abca3_matches_the_service_to_two_decimals(self, tmp_path):
        import re
        length = None
        with open(FASTA + ".fai") as fh:
            for line in fh:
                f = line.split("\t")
                if f[0] == "chr16":
                    length = f[1]
        vcf = tmp_path / "in.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n"
                       f"##contig=<ID=chr16,length={length}>\n"
                       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                       "chr16\t2317763\t.\tT\tA\t.\t.\t.\n")
        # ABCA3 is not a cancer gene, so it is not in the panel annotation. A one-row
        # annotation for it is committed beside this test, built by
        # scripts/build_spliceai_annotation.py from the same GENCODE v43 table the published
        # service uses -- which is what makes the comparison a comparison.
        ann = os.path.join(HERE, "data", "spliceai_annotation_abca3.txt")
        if not os.path.exists(ann):
            pytest.skip("no ABCA3 annotation fixture")
        out = tmp_path / "out.vcf"
        r = subprocess.run([SPLICEAI, "-I", str(vcf), "-O", str(out), "-R", FASTA,
                            "-A", ann, "-D", "500", "-M", "0"],
                           capture_output=True, text=True)
        if r.returncode != 0 or not out.exists():
            pytest.skip(f"local spliceai did not run: {r.stderr[-200:]}")
        m = re.search(r"SpliceAI=([^;\s]+)", out.read_text())
        assert m, "no SpliceAI annotation produced"
        p = m.group(1).split("|")
        local = {k: float(p[2 + i]) for i, k in enumerate(DS)}
        srv = best(frozen()["variants"]["16-2317763-T-A"]["scores"], "NM_001089.3")
        for k in DS:
            assert abs(local[k] - float(srv[k])) <= 0.01, (
                f"{k}: local {local[k]} against the service's {srv[k]}")
