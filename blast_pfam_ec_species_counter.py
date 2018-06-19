#! /usr/bin/evn python

"""A Python script to parse BLASTN/X results in tabular (6) format returning Pfam/EC/species counts per sample.

python blast_pfam_ec_species_counts.py <sample_id> <blast_output> <pfam_ec_species_db> <output_file>

Keyword arguments:
sample_id --> ID of the sample being analyzed
blast_output --> path to the BLASTN/X output file in tabular (6) format
pfam_ec_species_db --> path to the Pfam/EC/species database file in CSV with gene-Pfam, gene-EC, Sha384key-Pfam, Sha384key-EC, gene-species matches 
output_file --> desired name of the output file

Returns:
output_file --> TSV file (3 columns: sample, Pfam/EC/species, count)
"""

import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "28th of May 2018"

def blast_output_parser(blast_output):
    dic = {}
    with open(blast_output, "r") as thefile:
        for line in thefile:
            elements = line.strip().split("\t")
            gene_protein = elements[1]
            bit_score = float(elements[11])
            e_value = float(elements[10])
            if e_value <= 0.0000001 and bit_score >= 110: #set bitscore and e-value cutoff
                if gene_protein in dic:
                    dic[gene_protein] += 1
                else:
                    dic[gene_protein] = 1
    return dic

def count_pfam_ec_species(pfam_ec_species_db, blast_dic):
    dic = {}
    with open(pfam_ec_species_db, "r") as thefile:
        for line in thefile:
            gene_protein, pfam_ec_species = line.strip().split(",")
            if gene_protein in blast_dic:
                if pfam_ec_species in dic:
                    dic[pfam_ec_species] += (1*blast_dic[gene_protein])
                else:
                    dic[pfam_ec_species] = (1*blast_dic[gene_protein])
    return dic

def write_pfam_ec_species_counts(pfam_ec_species_count_dic, sample_id, output_file):
    with open(output_file, "w") as thefile:
        for pfam_ec_species in pfam_ec_species_count_dic:
            line = "\t".join([sample_id, pfam_ec_species, str(pfam_ec_species_count_dic[pfam_ec_species])])
            thefile.write(line + "\n")

if __name__ == "__main__":
    #Get variables from command line
    print "Loading arguments from command line..."
    #sample_id = sys.argv[1]
    output_dir = sys.argv[1]
    pfam_ec_species_db = sys.argv[2]
    #output_file = sys.argv[3]
    print "Arguments loaded"

    file_list = os.listdir(output_dir)
    for file in file_list:
        if file.startswith("NG-5"):
            sample_id = file.split(".")[0]
            output_file = "{}_species.tsv".format(sample_id) #set the extension
            file_name = "".join([output_dir, file])
            print sample_id, output_file

            #Parse the BLASTN/X output
            print "Parsing the BLASTN/X output..."
            blast_dic = blast_output_parser(file_name)
            print "BLASTN/X output parsed"

            #Make a Pfam/EC count dictionary
            print "Counting Pfam/EC/species occurrences..."
            pfam_ec_species_count_dic = count_pfam_ec_species(pfam_ec_species_db, blast_dic)
            print "Pfam/EC/species occurrences counted"

            #Write the count dictionary to a TSV
            print "Writing Pfam/EC/species occurrences to file..."
            write_pfam_ec_species_counts(pfam_ec_species_count_dic, sample_id, output_file)
            print "Pfam/EC/species occurrences written"
