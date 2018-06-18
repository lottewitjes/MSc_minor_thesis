#! /usr/bin/evn python

"""A Python script to assign Pfam domains to GO terms and count how many Pfam domains were assigned to that GO term.

python pfam2go_converter.py <pfam2go_mapping> <pfam_counts> <go_counts>

Keywords:
pfam2go_mapping --> a GO mapping file containing Pfam-GO term matches
pfam_counts --> a TSV file containing average (over samples) read counts for each Pfam domain
go_counts --> desired name of output file

Returns:
go_counts --> A TSV containing the number of Pfam domains that were mapped to each GO term
"""

import sys
import subprocess
import os.path

if __name__ == "__main__":
    pfam2go = sys.argv[1]
    pfam_counts = sys.argv[2]
    go_counts = sys.argv[3]


