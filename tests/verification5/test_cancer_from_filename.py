"""Three documented tumour types could never be matched from a filename.

`get_cancer_from_filename` split the basename on `[_.-]` and required a single part to
equal a CANCER_MAP key. Three of the 17 keys contain hyphens, so after splitting on `-`
no part could ever equal them:

    head-and-neck-squamous-cell-carcinoma-hnsc -> HNSC
    head-and-neck-snuc                         -> SNUC
    thyroid-carcinoma                          -> THCA

All three fell through to UNKNOWN. README names HNSC and SNUC as supported, spelled as
the token a user would put in a filename -- so the documentation told a user to do
something that could not work, and the failure was silent: an UNKNOWN tumour type makes
OncoKB answer pan-cancer instead of refusing.

Spans are now tried **longest first**, and that ordering is the part worth testing:
`thyroid` is THY and `thyroid-carcinoma` is THCA, two different OncoTree codes. Shortest
first would resolve `sample_thyroid-carcinoma.vcf` to THY and send OncoKB a tumour type
the sample does not have, changing every tumour-specific level it returns.
"""

import pytest

from conftest import oncokb


class TestTheHyphenatedKeysAreReachable:
    @pytest.mark.parametrize("filename,expected", [
        ("S1_head-and-neck-snuc.vcf", "SNUC"),
        ("S1_head-and-neck-squamous-cell-carcinoma-hnsc.vcf", "HNSC"),
        ("S1_thyroid-carcinoma.vcf", "THCA"),
    ])
    def test_they_resolve(self, filename, expected):
        assert oncokb.get_cancer_from_filename(filename) == expected

    def test_every_key_in_the_map_is_reachable(self):
        """The general form of the bug: a key nobody can match is a key that lies about
        what the pipeline supports."""
        unreachable = []
        for key, code in oncokb.CANCER_MAP.items():
            if oncokb.get_cancer_from_filename(f"sample_{key}.vcf") != code:
                unreachable.append(key)
        assert not unreachable, f"these keys cannot be matched from a filename: {unreachable}"


class TestTheLongestMatchWins:
    """`thyroid` -> THY and `thyroid-carcinoma` -> THCA are different tumour types."""

    def test_the_longer_key_is_not_shadowed_by_the_shorter(self):
        assert oncokb.get_cancer_from_filename("S1_thyroid-carcinoma.vcf") == "THCA"

    def test_the_shorter_key_still_works_alone(self):
        assert oncokb.get_cancer_from_filename("S1_thyroid.vcf") == "THY"

    def test_the_two_codes_really_do_differ(self):
        """If these ever become the same code the ordering stops mattering, and this
        test should be deleted rather than silently kept."""
        assert oncokb.CANCER_MAP["thyroid"] != oncokb.CANCER_MAP["thyroid-carcinoma"]


class TestTheOrdinaryCasesAreUnaffected:
    @pytest.mark.parametrize("filename,expected", [
        ("S1_lung.vcf", "LUNG"),
        ("S1_brain.vcf.gz", "CNS"),
        ("patient.breast.vcf", "BRCA"),
        ("S1_other.vcf", "OTHER"),
        ("S1_melanoma_rerun.vcf", "MEL"),
    ])
    def test_single_token_types_still_match(self, filename, expected):
        assert oncokb.get_cancer_from_filename(filename) == expected

    @pytest.mark.parametrize("filename", ["S1.vcf", "sample_001.vcf.gz", "", "no_type_here.vcf"])
    def test_an_unnamed_type_is_unknown(self, filename):
        assert oncokb.get_cancer_from_filename(filename) == "UNKNOWN"

    def test_a_directory_in_the_path_is_not_searched(self):
        """Only the basename should be read, or a sample sitting in a directory called
        `lung/` would be annotated as lung whatever it is."""
        assert oncokb.get_cancer_from_filename("/data/lung/S1.vcf") == "UNKNOWN"
