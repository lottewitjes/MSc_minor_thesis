#!/usr/bin/evn/python

"""A Python script that runs DIAMOND.

python run_diamond.py <input_directory> <database ><output_directory>

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
__version__ = "1.0"

def run_diamond(fasta_file, database, output_directory, output_file):
    """A function that runs DIAMOND on the input files.
    """
    cmd = "/metagenomics/lottewitjes/diamond blastx --query {} --db {} --threads 20 --max-target-seqs 1 --outfmt tab --out {}".format(fasta_file, database, output_file)
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
    input_directory = argv[1]
    database = argv[2]
    output_directory = argv[3]

    filelist = os.listdir(input_directory)
    for file in filelist:
        output_file = file.strip(".fastq")
        run_diamond(file, database, output_directory, output_file)

