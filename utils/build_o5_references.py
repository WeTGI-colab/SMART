#!/usr/bin/env python3

# =============================================================================
# build_o5_references.py
# =============================================================================
#
# DESCRIPTION:
#   Builds an optional table of (likely) oncogenic missense reference variants for
#   the SVIG-UK O5 missense route in svig_uk_classifier.py. The SVIG-UK canonical
#   variants list is always used by the classifier; this table adds the Cancer
#   Genome Interpreter (CGI) catalogue of validated oncogenic mutations, read from
#   the gene:protein dictionary distributed with OncoVI
#   (resources/cgi_dictionary.txt, https://github.com/MGCarta/oncovi).
#
#   CGI entries carry no transcript, so the classifier only uses them when the
#   residue number and reference amino acid agree with the transcript SMART
#   annotated. CGI data are subject to CGI's terms of use; keep the output with
#   your reference data.
#
# OUTPUT (tab-separated): gene transcript protein assessment source
# =============================================================================

import argparse
import json
import re

MISSENSE = re.compile(r"^p\.[A-Z]\d+[A-Z]$")


def main(argv=None):
    ap = argparse.ArgumentParser(description="O5 missense reference variants")
    ap.add_argument("--cgi-dictionary", required=True, help="OncoVI resources/cgi_dictionary.txt")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    with open(args.cgi_dictionary) as fh:
        entries = json.load(fh)
    n = 0
    with open(args.out, "w") as out:
        out.write("gene\ttranscript\tprotein\tassessment\tsource\n")
        for key in sorted(entries):
            gene, change = key.split(":", 1)
            if MISSENSE.match(change):
                out.write(f"{gene}\t\t{change}\tOncogenic\tCGI validated oncogenic mutations\n")
                n += 1
    print(f"wrote {n} missense references -> {args.out}")


if __name__ == "__main__":
    main()
