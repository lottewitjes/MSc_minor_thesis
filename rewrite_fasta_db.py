#!/usr/bin/evn/python

"""A Python script that rewrites a gene or protein database resulting from a SAPP SPARQL query.

python rewrite_fasta_db.py <input> <output>

Keyword arguments:
- input
- output
Returns:
- written output
"""

import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "2nd of May 2018"
__version__ = "1.0"

def parse_fasta(input):
    dic = {}
    with open(input, "r") as thefile:
        for line in thefile:
            if line.startswith('fasta'):
                continue
            if line.startswith('">'):
	        label = line.strip().replace('"', '')
		if label in dic:
		    next(thefile)
                else:
                    dic[label] = ""
            else:
	        dic[label] += line.strip().strip('"')

    return dic

def write_unique_fasta(fasta_dic, output):
    with open(output, "w") as thefile:
        for key in fasta_dic:
            thefile.write(str(key) + "\n")
	    thefile.write(str(fasta_dic[key]) + "\n")

if __name__ == "__main__":
    input = sys.argv[1]
    output = sys.argv[2]

    fasta_dic = parse_fasta(input)
    write_unique_fasta(fasta_dic, output)
