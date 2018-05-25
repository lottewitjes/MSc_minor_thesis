#! /usr/bin/evn python

"""A Python script to query a RDF graph returning a CSV with gene ID and sequence.

python rdf_query_gene_sequence.py <sample> <output_file>

Keyword arguments:
sample --> name of the sample of which sequences have to be extracted
output_file --> name of the CSV output file

Returns:
filled output_file --> CSV file with gene ID and gene sequence
"""

import sys
import subprocess
import os.path
from SPARQLWrapper import SPARQLWrapper, JSON

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "25th of May 2018"
__version__ = "1.0"

def run_query(graph, query):
    sparql = SPARQLWrapper(graph)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    return results

def write_sequence_file(query_output, output_file):
    with open(output_file, "w") as thefile:
        for result in query_output["results"]["bindings"]:
            gene = result["gene"]["value"]
            sequence = result["lcsequence"]["value"]
            line = ",".join([gene, sequence])
            thefile.write(line + "\n")

if __name__ == "__main__":
    graph_url = "http://10.117.11.77:7200/repositories/metagenomics_ileum"
    sample_id = sys.argv[1]
    output_file = sys.argv[2]

    sparql_query = "PREFIX gbol: <http://gbol.life/0.1/> PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> SELECT DISTINCT ?gene ?lcsequence WHERE { VALUES ?sample {<http://gbol.life/0.1/%s>} . ?sample a gbol:Sample . ?dnaobject gbol:sample ?sample . ?dnaobject gbol:feature ?gene . ?gene a gbol:Gene . ?gene gbol:transcript ?transcript . ?transcript gbol:sequence ?sequence . BIND (lcase(?sequence) as ?lcsequence) .}" %(str(sample_id))

    results = run_query(graph_url, sparql_query)

    write_sequence_file(results, output_file)
