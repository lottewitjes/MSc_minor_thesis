#! /usr/bin/evn python

"""A Python script to parse BLASTN/X results in tabular (6) format returning Pfam/EC counts per sample.

python blast_pfam_ec_counts.py <sample_id> <blast_output> <pfam_ec_db> <output_file>

Keyword arguments:
sample_id --> ID of the sample being analyzed
blast_output --> path to the BLASTN/X output file in tabular (6) format
pfam_ec_db --> path to the Pfam/EC database file in CSV with gene-Pfam, gene-EC, Sha384key-Pfam, Sha384key-EC matches 
output_file --> desired name of the output file

Returns:
output_file --> TSV file (3 columns: sample, Pfam/EC, count)
"""

import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "28th of May 2018"
__version__ = "1.0"

def pfam_ec_db_parser(pfam_ec_db):
    dic = {}
    with open(pfam_ec_db, "r") as thefile:
        for line in thefile:
            gene_protein, pfam_ec = line.strip().split(",")
            if gene_protein in dic:
                dic[gene_protein].append(pfam_ec)
            else:
                dic[gene_protein] = [pfam_ec]
    return dic

def blast_output_parser(blast_output):
    alist = []
    with open(blast_output, "r") as thefile:
        for line in thefile:
            elements = line.strip().split("\t")
            gene_protein = elements[1]
            bit_score = float(elements[11])
            e_value = float(elements[10])
            if e_value <= 0.0000001 and bit_score >= 74:
                alist.append(gene_protein)
    return alist

def count_pfam_ec(pfam_ec_dic, blast_list):
    dic = {}
    for element in blast_list:
        if element in pfam_ec_dic:
            for pfam_ec in pfam_ec_dic[element]:
                if pfam_ec in dic:
                    dic[pfam_ec] += 1
                else:
                    dic[pfam_ec] = 1
    return dic

def write_pfam_ec_counts(pfam_ec_count_dic, sample_id, output_file):
    with open(output_file, "w") as thefile:
        for pfam_ec in pfam_ec_count_dic:
            line = "\t".join([sample_id, pfam_ec, str(pfam_ec_count_dic[pfam_ec])])
            thefile.write(line + "\n")

if __name__ == "__main__":
    #Get variables from command line
    print "Loading arguments from command line..."
    #sample_id = sys.argv[1]
    output_dir = sys.argv[1]
    pfam_ec_db = sys.argv[2]
    #output_file = sys.argv[3]
    print "Arguments loaded"

    #Parse the Pfam/EC database
    print "Parsing the Pfam/EC database..."
    pfam_ec_dic = pfam_ec_db_parser(pfam_ec_db)
    print "Pfam/EC database parsed"

    file_list = os.listdir(output_dir)
    for file in file_list:
        if file.startswith("NG-5"):
            sample_id = file.split(".")[0]
            output_file = "{}_ec.tsv".format(sample_id)
            file_name = "".join([output_dir, file])
            print sample_id, output_file

            #Parse the BLASTN/X output
            print "Parsing the BLASTN/X output..."
            blast_list = blast_output_parser(file_name)
            print "BLASTN/X output parsed"

            #Make a Pfam/EC count dictionary
            print "Counting Pfam/EC occurrences..."
            pfam_ec_count_dic = count_pfam_ec(pfam_ec_dic, blast_list)
            print "Pfam/EC occurrences counted"

            #Write the count dictionary to a TSV
            print "Writing Pfam/EC occurrences to file..."
            write_pfam_ec_counts(pfam_ec_count_dic, sample_id, output_file)
            print "Pfam/EC occurrences written"
