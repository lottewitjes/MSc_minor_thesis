"""A Python script to parse BLASTN/X results in tabular (6) format returning EC counts per sample.

python blast_ec_counter.py <blast_output> <gene_protein_db>

"""

import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "26th of June 2018"

def blast_output_parser(blast_output):
    dic = {}
    with open(blast_output, "r") as thefile:
        for line in thefile:
            elements = line.strip().split("\t")
            gene_protein = elements[1]
            bit_score = float(elements[11])
            e_value = float(elements[10])
            if e_value <= 0.0000001 and bit_score >= 110: #set bitscore and e-value cutoff
                if gene_protein in dic:
                    dic[gene_protein] += 1
                else:
                    dic[gene_protein] = 1
    return dic


if __name__ == "__main__":

