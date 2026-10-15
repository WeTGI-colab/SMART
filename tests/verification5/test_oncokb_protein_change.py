"""
The protein change sent to OncoKB must be the variant on the row, not a prefix of it.

`oncokb2.0.py` built it with
`re.match(r"([A-Za-z*]{3})(\\d+)([A-Za-z*]{3})")` and returned `f"{ref1}{pos}{alt1}"`.
`re.match` anchors at the start but does not require reaching the end, so
`Trp288CysfsTer12` matched `Trp` + `288` + `Cys` and **the `fs*12` was discarded**. OncoKB
was then asked about `W288C`, a missense substitution at that codon, which is a different
variant.

Measured against the live API while diagnosing this:

    W288Cfs*12   Oncogenic, LEVEL_Dx1
    W288Cfs      Oncogenic, LEVEL_Dx1
    W288fs       Oncogenic, LEVEL_Dx1
    W288C        Unknown, "has not specifically been reviewed"

So NPM1's AML-defining diagnostic level was lost, and BRCA2 `S1982Rfs*22` -- a truncating
variant in a tumour suppressor -- came back as a VUS with "no available functional data about
the BRCA2 S1982R mutation". Every frameshift in every run was affected. It was invisible until
a frameshift was added to the test cohort, which had none.

A second defect in the same expression: anything that did not match that pattern was returned
unconverted, so `Met1673_Leu1676del` was sent to OncoKB with three-letter codes.

This file tests the conversion directly, so it does not need a token, a network or cyvcf2.
"""
import os
import re

import pytest

SRC = os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "oncokb2.0.py")


@pytest.fixture(scope="module")
def convert():
    """The function, lifted out of the module.

    `oncokb2.0.py` imports cyvcf2, which is in the image and not necessarily on a
    developer's machine, so importing the module would make this suite depend on an
    optional C extension. The function is pure, so its source is executed in isolation.
    """
    src = open(SRC).read()
    start = src.index("def _extract_one_letter_protein(")
    end = src.index("\ndef ", start + 10)
    ns = {"re": re}
    exec(compile(src[start:end], SRC, "exec"), ns)
    return lambda hgvsp: ns["_extract_one_letter_protein"](hgvsp, 0)


class TestFrameshifts:
    """The defect: a frameshift became a missense substitution at the same codon."""

    def test_the_npm1_type_a_frameshift_keeps_its_suffix(self, convert):
        assert convert("NP_1:p.Trp288CysfsTer12") == "W288Cfs*12"

    def test_a_truncating_brca2_frameshift_keeps_its_suffix(self, convert):
        assert convert("NP_1:p.Ser1982ArgfsTer22") == "S1982Rfs*22"

    def test_the_asxl1_frameshift_keeps_its_suffix(self, convert):
        assert convert("NP_1:p.Gly646TrpfsTer12") == "G646Wfs*12"

    @pytest.mark.parametrize("hgvsp,want", [
        ("NP_1:p.Trp288CysfsTer12", "W288Cfs*12"),
        ("NP_1:p.Ser1982ArgfsTer22", "S1982Rfs*22"),
    ])
    def test_the_substitution_prefix_alone_is_never_returned(self, convert, hgvsp, want):
        """What OncoKB received before: the first codon change, with the frameshift gone."""
        got = convert(hgvsp)
        assert got == want
        assert "fs" in got, f"the frameshift was dropped: {got!r}"
        assert got != want.split("fs")[0], "only the substitution prefix was returned"


class TestEveryOtherProteinForm:
    """The second defect: anything that was not a plain substitution kept three-letter codes."""

    @pytest.mark.parametrize("hgvsp,want", [
        ("NP_1:p.Val600Glu", "V600E"),
        ("NP_1:p.Arg248Ter", "R248*"),
        ("NP_1:p.Glu746_Ala750del", "E746_A750del"),
        ("NP_1:p.Tyr772_Ala775dup", "Y772_A775dup"),
        ("NP_1:p.Met1673_Leu1676del", "M1673_L1676del"),
        ("NP_1:p.His24_Ser37delinsPro", "H24_S37delinsP"),
        ("NP_1:p.Thr125=", "T125="),
    ])
    def test_the_form_is_converted_and_preserved(self, convert, hgvsp, want):
        assert convert(hgvsp) == want

    def test_no_three_letter_code_survives(self, convert):
        """A three-letter code reaching OncoKB is a query it will not recognise."""
        codes = ("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
                 "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
        for hgvsp in ("NP_1:p.Trp288CysfsTer12", "NP_1:p.His24_Ser37delinsPro",
                      "NP_1:p.Met1673_Leu1676del", "NP_1:p.Glu746_Ala750del"):
            got = convert(hgvsp)
            left = [c for c in codes if c in got]
            assert not left, f"{hgvsp} -> {got!r} still carries {left}"

    def test_the_lower_case_hgvs_keywords_are_untouched(self, convert):
        """del, ins, dup, fs and ext must survive: a three-letter code starts with a capital,
        which is the whole reason a case-sensitive whole-string replacement is safe here.
        An earlier attempt added lookarounds and they refused the `Cys` in `CysfsTer12`
        because `f` follows it, and the `Pro` in `delinsPro` because `s` precedes it."""
        assert convert("NP_1:p.Glu746_Ala750del").endswith("del")
        assert convert("NP_1:p.Tyr772_Ala775dup").endswith("dup")
        assert "delins" in convert("NP_1:p.His24_Ser37delinsPro")
        assert "fs" in convert("NP_1:p.Trp288CysfsTer12")


class TestTheEdges:
    def test_a_value_without_a_protein_change_gives_nothing(self, convert):
        """Falsy, not specifically None: the caller gates on `if gene and protein_change`,
        so '' and None are equivalent and asserting which would test an implementation
        detail rather than the behaviour that matters."""
        assert not convert("")
        assert not convert("NP_1:")

    def test_the_accession_prefix_is_stripped(self, convert):
        assert convert("ENSP00000493543.1:p.Val600Glu") == "V600E"

    def test_an_unprefixed_value_works_too(self, convert):
        assert convert("p.Val600Glu") == "V600E"
