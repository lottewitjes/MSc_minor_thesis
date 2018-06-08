#! /usr/bin/evn python

"""A Python script to parse BLASTN results in tabular (6) format returning phage counts

python phage_counter.py <blast_output> <output_file>

Keyword arguments:
blast_output --> BLASTN results in tabular (6) format
output_file --> desired name of output_file

Returns:
output_file --> a TSV with the following columns: RefSeq accession (NC), phage name, count
"""

import sys
import subprocess
import os.path
from Bio import Entrez

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "8th of June 2018"

def blast_output_parser(blast_output):
    dic = {}
    with open(blast_output, "r") as thefile:
        for line in thefile:
            elements = line.strip().split("\t")
            accession = elements[1].split("_")
            accession = "".join([accession[0], accession[1]])
            accession = accession.split(".")[0]
            bit_score = float(elements[11])
            e_value = float(elements[10])
            if e_value <= 0.0000001 and bit_score >= 110:
                if accession in dic:
                    dic[accession][0] += 1
                else:
                    dic[accession] = [1]
    return dic

def give_phage_names(blast_dic):
    Entrez.email = "lottewitjes@outlook.com"
    for phage in blast_dic:
        handle = Entrez.esummary(db="nucleotide", id=phage)
        record = Entrez.read(handle)
        blast_dic[phage].append(record[0]["Title"])
    return blast_dic

def write_phage_count_table(count_dic, output_file)
    with open(output_file, "w") as thefile:
        for phage in count_dic:
            line_elements = [phage, count_dic[phage][1], count_dic[phage][0]]
            line = "\t".join(line_elements)
            thefile.write(line + "\n")

if __name__ == "__main__":
    #Load arguments from command line
    blast_output = sys.argv[1]
    output_file = sys.argv[2]

    #Count phages
    blast_dic = blast_output_parser(blast_output)

    #Parse phages names into dic
    blast_dic = give_phage_names(blast_dic)

    #Write phage counts to TSV
    write_phage_count_table(blast_dic, output_file)

