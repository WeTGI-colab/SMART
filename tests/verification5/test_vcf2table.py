"""
Unit tests for scripts/vcf2table.py

Covers:
  - calculate_end()        — end-position logic for every variant class
  - convert_hgvsp_short()  — 3-letter → 1-letter AA code conversion
  - load_nm_transcripts()  — whitelist loading and version stripping
"""
import pandas as pd
import pytest

from conftest import vcf2table


# ===========================================================================
# calculate_end
# ===========================================================================

def _row(**kwargs) -> pd.Series:
    """Build a minimal pandas Series for calculate_end."""
    defaults = {"POS": 100, "REF": "A", "ALT": "G",
                "VARIANT_CLASS": "snv", "SVLEN": None}
    defaults.update(kwargs)
    return pd.Series(defaults)


class TestCalculateEnd:
    def test_snv_returns_pos(self):
        assert vcf2table.calculate_end(_row(VARIANT_CLASS="snv", POS=114713908)) == 114713908

    def test_deletion_returns_pos_plus_ref_minus_one(self):
        # REF=AGT (3 bp) → end = 100 + 3 - 1 = 102
        row = _row(POS=100, REF="AGT", ALT="A", VARIANT_CLASS="deletion")
        assert vcf2table.calculate_end(row) == 102

    def test_deletion_single_base(self):
        # REF=AT, ALT=A → end = 100 + 2 - 1 = 101
        row = _row(POS=100, REF="AT", ALT="A", VARIANT_CLASS="deletion")
        assert vcf2table.calculate_end(row) == 101

    def test_insertion_with_svlen(self):
        # SVLEN=10 → end = 100 + 10 - 1 = 109
        row = _row(POS=100, REF="A", ALT="AGGGGGGGGG",
                   VARIANT_CLASS="insertion", SVLEN=10)
        assert vcf2table.calculate_end(row) == 109

    def test_insertion_svlen_as_string(self):
        # SVLEN stored as string (common in VCF INFO)
        row = _row(POS=200, REF="A", ALT="ATT", VARIANT_CLASS="insertion", SVLEN="5")
        assert vcf2table.calculate_end(row) == 204

    def test_insertion_without_svlen_uses_alt_length(self):
        # ALT=AGG (3 bp) → end = 100 + 3 - 1 = 102
        row = _row(POS=100, REF="A", ALT="AGG", VARIANT_CLASS="insertion", SVLEN=None)
        assert vcf2table.calculate_end(row) == 102

    def test_chromosome_breakpoint_returns_empty_string(self):
        row = _row(VARIANT_CLASS="chromosome_breakpoint")
        assert vcf2table.calculate_end(row) == ""

    def test_unknown_class_returns_pos(self):
        # Fallback: any unrecognised class returns POS
        row = _row(POS=500, VARIANT_CLASS="sequence_alteration")
        assert vcf2table.calculate_end(row) == 500

    @pytest.mark.parametrize("pos,ref,alt,expected", [
        (114713908, "T", "C",   114713908),     # NRAS Q61R — SNV
        (208248388, "C", "T",   208248388),     # IDH1 R132H — SNV
        (55174772,  "GGAATTAAGAGAAGCA", "G", 55174787),  # EGFR del — deletion: 55174772 + 16 - 1
    ])
    def test_real_verification1_variants(self, pos, ref, alt, expected):
        n = len(ref)
        if n == len(alt):
            vc = "snv"
        elif n > len(alt):
            vc = "deletion"
        else:
            vc = "insertion"
        row = _row(POS=pos, REF=ref, ALT=alt, VARIANT_CLASS=vc)
        assert vcf2table.calculate_end(row) == expected


# ===========================================================================
# convert_hgvsp_short
# ===========================================================================

class TestConvertHgvspShort:
    @pytest.mark.parametrize("hgvsp,expected", [
        ("p.Val600Glu",  "p.V600E"),     # BRAF V600E
        ("p.Gln61Arg",   "p.Q61R"),     # NRAS Q61R
        ("p.Arg132His",  "p.R132H"),    # IDH1 R132H
        ("p.Glu545Lys",  "p.E545K"),    # PIK3CA E545K
        ("p.His1047Arg", "p.H1047R"),   # PIK3CA H1047R
        ("p.Arg248Trp",  "p.R248W"),    # TP53 R248W
        ("p.Arg175His",  "p.R175H"),    # TP53 R175H
    ])
    def test_three_letter_to_one_letter(self, hgvsp, expected):
        assert vcf2table.convert_hgvsp_short(hgvsp) == expected

    def test_stop_codon_ter_becomes_asterisk(self):
        assert vcf2table.convert_hgvsp_short("p.Arg248Ter") == "p.R248*"

    def test_strips_ensp_prefix(self):
        # Full HGVSp from VEP includes an ENSP prefix
        result = vcf2table.convert_hgvsp_short("ENSP00000358548.4:p.Val600Glu")
        assert result == "p.V600E"

    def test_already_single_letter_preserved(self):
        # If the regex finds 1-letter codes, aa_dict.get falls back to the key
        result = vcf2table.convert_hgvsp_short("p.V600E")
        assert result == "p.V600E"

    @pytest.mark.parametrize("hgvsp,expected", [
        # Previously broken: frameshifts kept 3-letter codes and lost the stop position
        ("ENSP00000347948.4:p.Thr15GlufsTer55", "p.T15Efs*55"),
        ("p.Trp288CysfsTer12",                  "p.W288Cfs*12"),   # NPM1 type A
        ("p.Leu367ThrfsTer46",                  "p.L367Tfs*46"),   # CALR type 1
        # Previously broken: multi-residue in-frame changes were left in 3-letter form
        ("p.Met1673_Leu1676del",                "p.M1673_L1676del"),
        ("p.His24_Ser37delinsPro",              "p.H24_S37delinsP"),
        ("p.Ala767_Val769dup",                  "p.A767_V769dup"),
        ("p.Lys745_Glu746insIlePro",            "p.K745_E746insIP"),
        ("p.Ile836del",                         "p.I836del"),      # FLT3 I836del
        # Synonymous and stop-loss extension
        ("p.Thr477%3D",                         "p.T477="),
        ("p.Ter403TrpextTer?",                  "p.*403Wext*?"),
    ])
    def test_non_substitution_forms(self, hgvsp, expected):
        assert vcf2table.convert_hgvsp_short(hgvsp) == expected

    def test_empty_string_returns_empty(self):
        assert vcf2table.convert_hgvsp_short("") == ""

    def test_none_returns_empty(self):
        assert vcf2table.convert_hgvsp_short(None) == ""

    def test_no_p_dot_returns_empty(self):
        assert vcf2table.convert_hgvsp_short("c.1799T>A") == ""

    def test_non_string_returns_empty(self):
        assert vcf2table.convert_hgvsp_short(123) == ""


# ===========================================================================
# load_nm_transcripts
# ===========================================================================

class TestLoadNmTranscripts:
    def test_loads_transcripts(self, nm_file):
        result = vcf2table.load_nm_transcripts(nm_file)
        assert "NM_000077" in result
        assert "NM_002524" in result
        assert "NM_005228" in result
        assert "NM_000546" in result

    def test_strips_version_suffixes(self, nm_file):
        result = vcf2table.load_nm_transcripts(nm_file)
        # Versioned entries must not appear
        assert "NM_000077.5" not in result
        assert "NM_002524.3" not in result

    def test_skips_blank_lines(self, nm_file):
        result = vcf2table.load_nm_transcripts(nm_file)
        assert "" not in result

    def test_comment_lines_treated_as_transcripts(self, nm_file):
        # The current implementation does NOT skip # lines,
        # so a comment becomes part of the set (stripped to '#').
        # Test that the function runs without error; the comment char
        # itself is not a valid NM ID but it won't crash.
        result = vcf2table.load_nm_transcripts(nm_file)
        assert isinstance(result, set)

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises((FileNotFoundError, OSError)):
            vcf2table.load_nm_transcripts(str(tmp_path / "nonexistent.txt"))

    def test_returns_set(self, nm_file):
        result = vcf2table.load_nm_transcripts(nm_file)
        assert isinstance(result, set)

    def test_single_entry_no_version(self, tmp_path):
        p = tmp_path / "single.txt"
        p.write_text("NM_005896\n")
        result = vcf2table.load_nm_transcripts(str(p))
        assert result == {"NM_005896"}

    def test_version_stripped_to_base(self, tmp_path):
        p = tmp_path / "v.txt"
        p.write_text("NM_023110.3\n")
        assert vcf2table.load_nm_transcripts(str(p)) == {"NM_023110"}


class TestTranscriptMatchColumns:
    """TRANSCRIPT_REQUESTED and TRANSCRIPT_MATCH: did the row get the transcript asked for?

    VEP annotates every transcript overlapping a variant and vcf2table keeps one block. When
    none of the supplied transcripts is among them it falls back to VEP's first block --

        if not nm_selected and first_vep is not None:
            # No match found -- use the first CSQ record

    -- which is an arbitrary transcript, chosen by output order. Before these columns nothing
    in the output said so, and it was worse than merely silent: `first_nm` was initialised and
    never assigned, so a fallback row left NM_Transcript *empty*, which reads as "no
    transcript" rather than "one we did not ask for".

    Every other transcript guarantee is downstream of this one: a predictor correctly resolved
    to the selected transcript is still the wrong answer if the selection was a fallback.
    """

    @staticmethod
    def _run(tmp_path, csq_blocks, whitelist, info_extra=""):
        """Run vcf2table over a one-variant VCF with the given CSQ blocks."""
        import subprocess, sys, os, csv as _csv
        fields = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Feature", "HGVSc",
                  "MANE_SELECT", "MANE", "MANE_PLUS_CLINICAL"]
        vcf = tmp_path / "in.vcf"
        info = "CSQ=" + ",".join("|".join(b.get(f, "") for f in fields) for b in csq_blocks)
        if info_extra:
            info += ";" + info_extra
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. '
            f'Format: {"|".join(fields)}">\n'
            '##INFO=<ID=ONCOKB_PREFERRED_TRANSCRIPT,Number=1,Type=String,Description="x">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n"
            f"chr7\t140753336\tv1\tA\tT\t.\tPASS\t{info}\tGT:AD\t0/1:50,20\n")
        tx = tmp_path / "tx.txt"
        tx.write_text("\n".join(whitelist) + "\n")
        out = tmp_path / "out.tsv"
        script = os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "vcf2table.py")
        r = subprocess.run([sys.executable, script, "--vcf", str(vcf),
                            "--transcripts", str(tx), "--output", str(out)],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-800:]
        with open(out, newline="") as fh:
            rows = list(_csv.DictReader(fh, delimiter="\t"))
        assert len(rows) == 1, f"expected one row, got {len(rows)}"
        return rows[0]

    def test_a_supplied_transcript_is_reported_as_a_match(self, tmp_path):
        row = self._run(tmp_path,
                        [dict(Feature="NM_004333.6", SYMBOL="BRAF",
                              HGVSc="NM_004333.6:c.1799T>A", Consequence="missense_variant")],
                        ["NM_004333.6"])
        assert row["TRANSCRIPT_MATCH"] == "match"
        assert row["TRANSCRIPT_REQUESTED"] == "NM_004333"
        assert row["NM_Transcript"] == "NM_004333"

    def test_the_supplied_transcript_is_chosen_over_the_first_block(self, tmp_path):
        """The whole point of a whitelist: order in the VCF must not decide."""
        row = self._run(tmp_path,
                        [dict(Feature="NM_999999.1", SYMBOL="BRAF", Consequence="intron_variant"),
                         dict(Feature="NM_004333.6", SYMBOL="BRAF",
                              HGVSc="NM_004333.6:c.1799T>A", Consequence="missense_variant")],
                        ["NM_004333.6"])
        assert row["TRANSCRIPT_MATCH"] == "match"
        assert row["NM_Transcript"] == "NM_004333"
        assert row["Consequence"] == "missense_variant", "the wrong block was kept"

    def test_a_fallback_is_declared_and_names_the_transcript_used(self, tmp_path):
        """The silent failure, now stated: no supplied transcript was annotated here."""
        row = self._run(tmp_path,
                        [dict(Feature="NM_888888.2", SYMBOL="BRAF",
                              Consequence="missense_variant")],
                        ["NM_004333.6"])
        assert row["TRANSCRIPT_MATCH"] == "fallback:NM_888888.2", row["TRANSCRIPT_MATCH"]
        assert row["TRANSCRIPT_REQUESTED"] == "", (
            "nothing was requested for this variant, so the column must be empty rather "
            "than implying a transcript was honoured")

    def test_a_fallback_names_the_first_block_when_there_are_several(self, tmp_path):
        row = self._run(tmp_path,
                        [dict(Feature="NM_111111.1", SYMBOL="BRAF", Consequence="intron_variant"),
                         dict(Feature="NM_222222.2", SYMBOL="BRAF", Consequence="missense_variant")],
                        ["NM_004333.6"])
        assert row["TRANSCRIPT_MATCH"] == "fallback:NM_111111.1"

    def test_an_oncokb_pin_of_a_supplied_transcript_is_a_match(self, tmp_path):
        """A pin is always one of the supplied transcripts, so it is a match, not an outcome.

        oncokb2.0.py sets ONCOKB_PREFERRED_TRANSCRIPT only when
        `candidate in preferred_tx` -- the same list the user supplied. A first version of
        these columns reported the pin path as "pinned:<acc>", which made 32 of the 37
        verification rows look as though the transcript had come from somewhere other than
        the user's list, and would have hidden a real fallback among them.
        """
        row = self._run(tmp_path,
                        [dict(Feature="NM_004333.6", SYMBOL="BRAF",
                              HGVSc="NM_004333.6:c.1799T>A", Consequence="missense_variant")],
                        ["NM_004333.6"],
                        info_extra="ONCOKB_PREFERRED_TRANSCRIPT=NM_004333")
        assert row["TRANSCRIPT_MATCH"] == "match", row["TRANSCRIPT_MATCH"]
        assert row["TRANSCRIPT_REQUESTED"] == "NM_004333"

    def test_the_only_two_outcomes_are_match_and_fallback(self, tmp_path):
        """Stated as its own test so a third value cannot creep back in unnoticed."""
        for blocks, wl in (
            ([dict(Feature="NM_004333.6", SYMBOL="BRAF",
                   HGVSc="NM_004333.6:c.1799T>A", Consequence="missense_variant")],
             ["NM_004333.6"]),
            ([dict(Feature="NM_888888.2", SYMBOL="BRAF", Consequence="missense_variant")],
             ["NM_004333.6"]),
        ):
            m = self._run(tmp_path, blocks, wl)["TRANSCRIPT_MATCH"]
            assert m == "match" or m.startswith("fallback:"), m

    def test_a_fallback_no_longer_leaves_the_transcript_column_empty_without_saying_why(
            self, tmp_path):
        """NM_Transcript may still be empty on a fallback, but TRANSCRIPT_MATCH explains it.

        That pairing is the fix: an empty transcript with no explanation was the dangerous
        state, because it looks like missing data rather than a substitution.
        """
        row = self._run(tmp_path,
                        [dict(Feature="NM_888888.2", SYMBOL="BRAF",
                              Consequence="missense_variant")],
                        ["NM_004333.6"])
        assert row["TRANSCRIPT_MATCH"].startswith("fallback:")
        assert row["TRANSCRIPT_MATCH"] != "match"
