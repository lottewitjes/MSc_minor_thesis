#! /usr/bin/evn python

"""A Python script to parse a SILVA db to extract IDs with accompanying genus.

python silva_id_genus_parser.py <silva_db> <output_name>

Keyword arguments:
silva_db -->
output_name -->

Returns:
output_name -->
"""

import sys
import subprocess
import os.path

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "4th of July 2018"

def silva_parser(silva_db):
    dic = {}
    with open(silva_db, "r") as thefile:
        for line in thefile:
            if line.startswith(">"):
                elements = line.strip().split(" ")
                if elements[1].startswith("Bacteria") or elements[1].startswith("Archaea"):
                    silva_id = elements[0].strip(">")
                    genus = elements[1].split(";")[-2].strip("[").strip("]")
                    dic[silva_id] = genus
    return dic

def write_genus_db(genus_dic, output_name):
    with open(output_name, "w") as thefile:
        for element in genus_dic:
            line = ",".join([element, genus_dic[element]])
            thefile.write(line + "\n")

if __name__ == "__main__":
    #Get arguments from command line
    silva_db = sys.argv[1]
    output_file = sys.argv[2]

    #Parse SILVA db
    silva_db_dic = silva_parser(silva_db)

    #Write SILVA ID, genus name to CSV
    write_genus_db(silva_db_dic, output_file)
