#! /usr/bin/evn python

"""A Python script to query a RDF graph returning enzyme/domain counts per sample.

python rdf_query_enzyme_domain_counts.py <query_type> <output_file>

Keyword arguments:
query_type --> pfam/ec for domain or enzyme counts
output_file --> filename for output

Returns:
filled output_file --> tab-separated file with sample, pfam/ec accession, count

"""

import sys
import subprocess
import os.path
from SPARQLWrapper import SPARQLWrapper, JSON

__author__ = "Lotte Witjes"
__email__ = "lottewitjes@outlook.com"
__date__ = "22th of May 2018"
__version__ = "1.0"

def run_query(graph, query):
    sparql = SPARQLWrapper(graph)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    return results

def write_count_file(query_output, output_file):
    with open(output_file, "w") as thefile:
        for result in query_output["results"]["bindings"]:
            sample = result["sample"]["value"].split("/")[-1]
            acc = result["acc"]["value"]
            acc_count = result["acc_count"]["value"]
            line = "\t".join([sample, acc, acc_count])
            thefile.write(line + "\n")

if __name__ == "__main__":
    graph_url = "http://10.117.11.77:7200/repositories/metagenomics_ileum"
    query_type = sys.argv[1] #pfam/ec
    output_file = sys.argv[2]

    sparql_query = "PREFIX gbol: <http://gbol.life/0.1/> PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> SELECT ?sample ?acc (COUNT(?acc) as ?acc_count) WHERE { VALUES ?db {<http://gbol.life/0.1/db/%s>} ?sample a gbol:Sample . ?dnaobject gbol:sample ?sample . ?dnaobject gbol:feature ?gene . ?gene a gbol:Gene . ?gene gbol:transcript/gbol:feature ?cds . ?cds gbol:protein ?protein . ?protein gbol:xref ?xref . ?xref gbol:db ?db . ?xref gbol:accession ?acc . } GROUP BY ?sample ?acc ORDER BY ?sample ?acc" % (str(query_type))

    results = run_query(graph_url, sparql_query)

    write_count_file(results, output_file)
