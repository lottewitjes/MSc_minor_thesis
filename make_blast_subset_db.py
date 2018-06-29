#! /usr/bin/evn python

"""A Python script to parse BLASTN/X results in tabular (6) format to create a subset sequence database.

python make_blast_subset_db.py <blast_result_dir> <gene_protein_db> <subset_db>

Keyword arguments:
blast_result_dir --> PATH to directory containing BLASTN/X results in format 6
gene_protein_db --> gene/protein sequence db in FASTA format
subset_db --> name of the output FASTA file containing the subset

Returns:
subset_db --> subset gene/protein sequence db in FASTA format based on BLASTN/X results

"""

import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "27th of June 2018"

def blast_output_parser(blast_output):
    blast_set = set()
    with open(blast_output, "r") as thefile:
        for line in thefile:
            elements = line.strip().split("\t")
            gene_protein = elements[1]
            bit_score = float(elements[11])
            e_value = float(elements[10])
            if e_value <= 0.0000001 and bit_score >= 74: #set bitscore and e-value cutoff
                if not gene_protein in blast_set:
                    blast_set.add(gene_protein)
    return blast_set

def write_subset_db(blast_total, gene_protein_db, subset_db):
    with open(gene_protein_db, "r") as f1, open(subset_db, "w") as f2:
        for line in f1:
            gene_protein, sequence = line.strip().split(",")[0:2]
            if gene_protein in blast_total:
                header = ">{}".format(gene_protein)
                f2.write(header + "\n" + sequence + "\n")

if __name__ == "__main__":
    #Get variables from command line
    blast_result_dir = sys.argv[1]
    gene_protein_db = sys.argv[2]
    subset_db = sys.argv[3]

    #Make a set of genes/proteins from the BLASTN/X results
    blast_total = set()
    file_list = os.listdir(blast_result_dir)
    for file in file_list:
        if file.startswith("NG-5"):
            file_name = "".join([blast_result_dir, file])
            print file_name
            blast_subset = blast_output_parser(file_name)
            blast_total = blast_total.union(blast_subset)

    #Make a subset db containing only those gene/protein sequences present in the BLASTN/X results
    subset_db = "".join([blast_result_dir, subset_db])
    write_subset_db(blast_total, gene_protein_db, subset_db)
