#!/usr/bin/evn/python

"""A Python script that runs Qtrim.
python run_qtrim_batch.py <input_directory> <output_directory>
Keyword arguments:
- input_director
- output_directory
Returns:
- filled output_directory containing Qtrim results
"""

import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "18th of April 2018"
__version__ = "1.0"

def run_qtrim(file, output_directory):
    """A function that runs Qtrim on the input files.
    """
    file_output = str(output_directory) + str(file.split(".")[0]) + "_trimmed.fastq"
    cmd = "/metagenomics/lottewitjes/programs/QTrim_v1_1/QTrim_v1_1 -fastq {} -o {}".format(file, file_output)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    try:
        results = subprocess.check_call(cmd, shell=True)
        return results
    except subprocess.CalledProcessError as err:
        print err.output
        sys.exit()

if __name__ == "__main__":
    input_directory = sys.argv[1]
    output_directory = sys.argv[2]

    filelist = os.listdir(input_directory)
    for file in filelist:
        if not file.startswith("fastqc"):
            run_qtrim(file, output_directory)
