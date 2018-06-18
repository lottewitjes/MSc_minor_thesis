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

def parse_pfam2go(pfam2go_mapping):
    dic = {}
    with open(pfam2go_mapping, "r") as thefile:
        for line in thefile:
            if not line.startswith("!"):
                pfam, go_term = line.strip().split(">")
                pfam = pfam.split("Pfam:")[1].split(" ")[0]
                go_term = go_term.split(";")[1].strip()
                if pfam in dic:
                    dic[pfam].append(go_term)
                else:
                    dic[pfam] = [go_term]
    return dic

def go_counter(pfam2go_dic, pfam_counts):
    dic = {}
    with open(pfam_counts, "r") as thefile:
        for line in thefile:
            pfam = line.strip().split("\t")[0]
            if pfam in pfam2go_dic:
                for go_term in pfam2go_dic[pfam]:
                    if go_term in dic:
                        dic[go_term] += 1
                    else:
                        dic[go_term] = 1
    return dic

def count_writer(go_count_dic, go_counts):
    alist = []
    for key in go_count_dic:
        go_count = [key, go_count_dic[key]]
        alist.append(go_count)
    alist.sort(key=lambda x: x[1], reverse=True)
    with open(go_counts, "w") as thefile:
        for element in alist:
            line = "\t".join([element[0], element[1]])
            thefile.write(line + "\n")

if __name__ == "__main__":
    pfam2go = sys.argv[1]
    pfam_counts = sys.argv[2]
    go_counts = sys.argv[3]

    pfam2go_dic = parse_pfam2go(pfam2go)

    go_count_dic = go_counter(pfam2go_dic, pfam_counts)

    count_writer(go_count_dic, go_counts)
