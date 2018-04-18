#!/usr/bin/evn/python

"""A Python script to run fastq-dump on a list of accession contained in a txt file

python fastq_dump_batch.py <accessions_file>

"""

import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "17th of April 2018"
__version__ = "1.0"

def run_fastq_dump(accession):
    cmd = "/metagenomics/lottewitjes/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump {}".format(accession)
    try:
        results = subprocess.check_call(cmd, shell=True)
        return results
    except subprocess.CalledProcessError as err:
        print err.output
        sys.exit()

if __name__ == "__main__":
    accessions_file = sys.argv[1]

    with open(accessions_file, "r") as thefile:
        for line in thefile:
            run_fastq_dump(line.strip())
