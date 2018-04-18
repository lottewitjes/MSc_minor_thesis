#!/usr/bin/evn/python

"""A Python script that runs Qtrim.
python run_diamond.py <input_directory> <output_directory>
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
    """A function that runs DIAMOND on the input files.
    """
    file_output = str(output_directory) + str(file.split(".")[0]) + "_trimmed.fastq"
    cmd = "/metagenomics/lottewitjes/programs/QTrim_v1_1/QTrim_v1_1 -fastq {} -o {}".format(file, file_output)
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

    filelist = os.listdir(input_directory)
    for file in filelist:
        run_qtrim(file, output_directory)
