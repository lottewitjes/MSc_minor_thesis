#! /usr/bin/evn python

"""A Python script to parse and filter EnzDP results, making a CSV with gene/protein and EC number.

python enzdp_result_parser.py <enzdp_results> <filtered_enzdp_results>

Keyword arguments:
enzdp_results --> a CSV or TSV with columns for gene/protein, EC number, likelihood score and max bitscore as given by EnzDP
filtered_enzdp_results --> the name of the output file containing filtered EnzDP results

Returns:
filtered_enzdp_results --> a CSV with columns for gene/protein and EC number

"""

import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "29th of June 2018"

def enzdp_parser(enzdp_res):
    alist = []
    with open(enzdp_res, "r") as thefile:
        for line in thefile:
            if not line.startswith(">") or line.startswith("gene"): #skip headers
                gene_protein, ec_number, likelihood, bitscore = line.strip.split(",") #CSV or TSV
                likelihood = float(likelihood)
                bitscore = float(bitscore)
                if ec_number.split(".")[-1] != "-" and likelihood >= 0.3 and bitscore >= 74:
                    approved = [gene_protein, ec_number]
                    alist.append(approved)
    return alist

if __name__ == "__main__":
    #Get arguments from command line
    enzdp_results = sys.argv[1]
    filtered_enzdp_results_filename = sys.argv[2]

    #Parse and filter EnzDP results

    #Write filtered EnzDP results
