#!/usr/bin/evn/python

"""A Python script that runs FastQC.

python run_fastqc_batch.py <input directory> <output_directory>

Keyword arguments:
- input_directory
- output_directory
Returns:
- filled output_directory containing filtered read files
"""

import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "18th of April 2018"
__version__ = "1.0"

def run_fastqc(file_list, output_directory):
    file_string = " ".join(file_list)
    cmd = "/metagenomics/lottewitjes/programs/FastQC/fastqc {} -o {}".format(file_string, output_directory)
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
    output_directory = argv[2]

    file_list = os.listdir(input_directory)
    run_fastqc(file_list, output_directory)
