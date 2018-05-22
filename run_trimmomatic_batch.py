#!/usr/bin/evn/python

"""A Python script that runs Trimmomatic.

python run_trimmomatic_batch.py <input_directory> <output_directory>

Keyword arguments:
- input_directory
- output_directory
Returns:
- filled output_directory containing Trimmomatic results
"""

import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "18th of April 2018"
__version__ = "1.0"

def run_trimmomatic(file, output_directory):
    file_output = str(output_directory) + str(file.split(".")[0]) + "_trimmed.fastq"
    cmd = "java -jar /metagenomics/lottewitjes/programs/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 4 -phred33 {} {} ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36".format(file, file_output)
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

    file_list = os.listdir(input_directory)
    for file in file_list:
        if not file.startswith("NG-5593_4") and not file == "TruSeq3-SE" and not file.startswith("fastqc"):
            run_trimmomatic(file, output_directory)


