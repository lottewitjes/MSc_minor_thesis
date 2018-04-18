#!/usr/bin/evn/python

"""A Python script that runs DIAMOND.

python run_diamond.py <blast_type> <input_directory> <database ><output_directory>

Keyword arguments:
- input_directory
- database
- output_directory
Returns:
- filled output_directory containing DIAMOND results
"""

import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "17th of April 2018"
__version__ = "1.0"

def run_diamond(blast_type, fasta_file, database, output_directory, output_file):
    """A function that runs DIAMOND on the input files.
    """
    cmd = "/metagenomics/lottewitjes/programs/diamond {} --query {} --db {} --threads 20 --max-target-seqs 1 --outfmt tab --out {}".format(blast_type, fasta_file, database, output_file)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    else:
        try:
            results = subprocess.check_call(cmd, shell=True)
            return results
        except subprocess.CalledProcessError as err:
            print err.output
            sys.exit()

if __name__ == "__main__":
    blast_type = argv[1]
    input_directory = argv[2]
    database = argv[3]
    output_directory = argv[4]

    filelist = os.listdir(input_directory)
    for file in filelist:
        output_file = file.strip(".fastq")
        run_diamond(blast_type, file, database, output_directory, output_file)

