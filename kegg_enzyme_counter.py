#! /usr/bin/evn python

"""A Python script to count how many enzymes in a given KEGG pathway are present in your favourite dataset.

python kegg_enzyme_counter.py <target_dir> <enzyme_count_file> <output_file>

Keywords:
target_dir --> the path to the directory containing the KEGG pathways in XML format
enzyme_count_file --> a TSV file with the following columns: sample_id, EC_number, color_code
output_file --> desired name of output file

Returns:
output_file --> a TSV with the following columns: pathway name, mapped_enzymes, enzymes_in_pathway

"""
import sys
import subprocess
import os.path
import xml.etree.ElementTree as ET

def parse_pathway_xml(dir_xml):
    dic = {}
    for file in os.listdir(dir_xml):
	if file.endswith(".xml"):
            tree = ET.parse("{}/{}".format(dir_xml, file))
            root = tree.getroot()
            dic[root.attrib["title"]] = []
            for child in root:
	        if child.tag == "entry":
                    if child.attrib["name"].startswith("ec"):
                        enzyme_name = child.attrib["name"].strip("ec:")
			if "ec:" in enzyme_name:
			    enzymes = enzyme_name.split(" ")
			    enzymes = [enzyme.strip("ec:") for enzyme in enzymes]
			    dic[root.attrib["title"]].extend(enzymes)
			else:
	                    dic[root.attrib["title"]].append(enzyme_name)
    return dic

def parse_ec_count(ec_file):
    alist = []
    with open(ec_file, "r") as thefile:
        for line in thefile:
            ec = line.split("\t")[1]
            alist.append(ec)
    return alist

def count_mapped_total_enzymes(pathway_dic, ec_list):
    dic = {}
    for pathway in pathway_dic:
        set1 = set(pathway_dic[pathway])
        set2 = set(ec_list)
        intersection = set1.intersection(set2)
        dic[pathway] = [len(intersection), len(set1)] #len(pathway_dic[pathway])]
    return dic

def write_tsv(count_dic, output_file):
    with open(output_file, "w") as thefile:
        line = "\t".join(["pathway", "present_enzymes", "total_enzymes"])
        thefile.write(line + "\n")
        for pathway in count_dic:
            line = "\t".join([pathway, str(count_dic[pathway][0]), str(count_dic[pathway][1])])
            thefile.write(line + "\n")

if __name__ == "__main__":
    target_dir = sys.argv[1]
    enzyme_count_file = sys.argv[2]
    output_file = sys.argv[3]

    pathway_dic = parse_pathway_xml(target_dir)

    ec_list= parse_ec_count(enzyme_count_file)

    pathway_count_dic = count_mapped_total_enzymes(pathway_dic, ec_list)

    write_tsv(pathway_count_dic, output_file)

    kegg_enzymes = set()
    for pathway in pathway_dic:
	unique_enzymes = set(pathway_dic[pathway])
	kegg_enzymes = kegg_enzymes.union(unique_enzymes)
    print kegg_enzymes
    print len(kegg_enzymes)

