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

domain_query = """PREFIX gbol: <http://gbol.life/0.1/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
SELECT ?sample ?acc (COUNT(?acc) as ?acc_count) WHERE {
    VALUES ?db {<http://gbol.life/0.1/db/pfam>}
    ?sample a gbol:Sample .
    ?dnaobject gbol:sample ?sample .
    ?dnaobject gbol:feature ?gene .
    ?gene a gbol:Gene .
    ?gene gbol:transcript/gbol:feature ?cds .
    ?cds gbol:protein ?protein .
    ?protein gbol:xref ?xref .
    ?xref gbol:db ?db .
    ?xref gbol:accession ?acc .
}
GROUP BY ?sample ?acc
ORDER BY ?sample ?acc
"""

enzyme_query = """PREFIX gbol: <http://gbol.life/0.1/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
SELECT ?sample ?acc (COUNT(?acc) as ?signature_count) WHERE {
    VALUES ?db {<http://gbol.life/0.1/db/ec>}
    ?sample a gbol:Sample .
    ?dnaobject gbol:sample ?sample .
    ?dnaobject gbol:feature ?gene .
    ?gene a gbol:Gene .
    ?gene gbol:transcript/gbol:feature ?cds .
    ?cds gbol:protein ?protein .
    ?protein gbol:xref ?xref .
    ?xref gbol:db ?db .
    ?xref gbol:accession ?acc .
}
GROUP BY ?sample ?acc
ORDER BY ?sample ?acc
"""

def sparql_query(graph, query):
    sparql = SPARQLWrapper(graph)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    return results

def write_count_matrix(query_output, output_file):
    results_dic = {}
    for result in query_output["results"]["bindings"]:
        sample = result["sample"]["value"].split("/")[-1]
        if sample not in results_dic:
            results_dic[sample] = []
        acc = result["acc"]["value"]
        acc_count = result["acc_count"]["value"]
        results_dic[sample].append([acc, acc_count])
    print results_dic["NG-5450_A"]


if __name__ == "__main__":
    graph_url = "http://10.117.11.77:7200/repositories/metagenomics_ileum"
    query_type = sys.argv[1] #domain or enzyme
    output_file = sys.argv[2]

    if query_type == "domain":
        results = sparql_query(graph_url, domain_query)
    elif query_type == "enzyme":
        results = sparql_query(graph_url, enzyme_query)

    write_count_matrix(results, output_file)
