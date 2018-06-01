#! /usr/bin/evn python

import sys
import subprocess
import os.path
import xml.etree.ElementTree as ET

def parse_pathway_xml(dir_xml):
    dic = {}
    for file in os.listdir(dir_xml):
	if file.endswith(".xml"):
            tree = ET.parse(file)
            root = tree.getroot()
            dic[root.attrib["title"]] = []
            for child in root:
	        if child.tag == "entry":
                    if child.attrib["name"].startswith("ec"):
                        enzyme = child.attrib["name"].strip("ec:")
	                dic[root.attrib["title"]].append(enzyme)
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
        dic[pathway] = [len(intersection), len(pathway_dic[pathway])]
    return dic

def write_tsv(count_dic):
    with open("pathway_counts.tsv", "w") as thefile:
        line = "\t".join(["pathway", "present_enzymes", "total_enzymes"])
        thefile.write(line + "\n")
        for pathway in count_dic:
            line = "\t".join([pathway, str(count_dic[pathway][0]), str(count_dic[pathway][1])])
            thefile.write(line + "\n")

if __name__ == "__main__":
    target_dir = sys.argv[1]
    enzyme_count_file = sys.argv[2]

    pathway_dic = parse_pathway_xml(target_dir)
    ec_list= parse_ec_count(enzyme_count_file)

    pathway_count_dic = count_mapped_total_enzymes(pathway_dic, ec_list)

    write_tsv(pathway_count_dic)

