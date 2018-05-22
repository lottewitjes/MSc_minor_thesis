#! /usr/bin/evn python

"""A Python script to query a RDF graph returning enzyme/domain counts per sample.

python rdf_query_enzyme_domain_counts.py

"""

import sys
import subprocess
import os.path
from SPARQLWrapper import SPARQLWrapper, JSON

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "22th of May 2018"
__version__ = "1.0"

if __name__ == "__main__":

