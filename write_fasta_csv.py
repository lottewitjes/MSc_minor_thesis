#!/usr/bin/evn/python

"""A Python script that writes a fasta file from a gene or protein csv database resulting from a SAPP SPARQL query.

python write_fasta_csv.py <input> <output>

Keyword arguments:
- input --> csv where the first column is the header and the second column the sequence
- output --> fasta file name

Returns:
- written output --> csv converted to fasta
"""

import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "2nd of May 2018"
__version__ = "1.0"

def parse_csv_write_fasta(input, output):
    with open(input, "r") as file1:
        with open(output, "w") as file2:
            next(file1) #skip the column headers
            for line in file1:
                header, sequence = line.split(",")
                header = header.strip('"')
                sequence = sequence.strip('"')
                file2.write(">" + str(header) + "\n" + str(sequence))

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    parse_csv_write_fasta(input_file, output_file)
