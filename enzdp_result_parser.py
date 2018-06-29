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
    dic = {}
    with open(enzdp_res, "r") as thefile:
        for line in thefile:
            if not line.startswith(">") and not line.startswith("gene"): #skip headers
                gene_protein, ec_number, likelihood, bitscore = line.strip().split(",") #CSV or TSV
                likelihood = float(likelihood)
                bitscore = float(bitscore)
                if ec_number.split(".")[-1] != "-" and likelihood >= 0.1 and bitscore >= 74: #specify thresholds
                    if gene_protein in dic:
                        if likelihood >= dic[gene_protein][1] and bitscore >= dic[gene_protein][2]:
                            dic[gene_protein] = [ec_number, likelihood, bitscore]
                    else:
                        dic[gene_protein] = [ec_number, likelihood, bitscore]
    return dic

def write_filtered_results(filtered_enzdp_res, filename):
    with open(filename, "w") as thefile:
        for element in filtered_enzdp_res:
            line = ",".join([element, filtered_enzdp_res[element][0]])
            thefile.write(line + "\n")

if __name__ == "__main__":
    #Get arguments from command line
    enzdp_results = sys.argv[1]
    filtered_enzdp_results_filename = sys.argv[2]

    #Parse and filter EnzDP results
    filtered_enzdp_results = enzdp_parser(enzdp_results)
    print len(filtered_enzdp_results)

    #Write filtered EnzDP results
    write_filtered_results(filtered_enzdp_results, filtered_enzdp_results_filename)
