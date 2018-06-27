#! /usr/bin/evn python

"""A Python script to parse BLASTN/X results in tabular (6) format to create a subset sequence database.

python make_blast_subset_db.py <blast_result_dir> <gene_protein_db>

Keyword arguments:

Returns:

"""

import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "27th of June 2018"

def blast_output_parser(blast_output):
    blast_set = set()
    with open(blast_output, "r") as thefile:
        for line in thefile:
            elements = line.strip().split("\t")
            gene_protein = elements[1]
            bit_score = float(elements[11])
            e_value = float(elements[10])
            if e_value <= 0.0000001 and bit_score >= 110: #set bitscore and e-value cutoff
                if not gene_protein in blast_set
                    blast_set.add(gene_protein)
    return blast_set

if __name__ == "__main__":
    #Get variables from command line
    blast_result_dir = sys.argv[1]
    gene_protein_db = sys.argv[2]

    blast_total = set()
    file_list = os.listdir(blast_result_dir)
    for file in file_list:
        if file.startswith("NG-5"):
            file_name = "".join([blast_result_dir, file])
            blast_subset = blast_output_parser(file_name)
            print len(blast_subset)
            blast_total = blast_total.union(blast_subset)
    print len(blast_total)
