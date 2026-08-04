"""
EpiMatch GEO Discovery

Searches NCBI GEO for datasets matching a keyword and organism,
using the public E-utilities API. Returns candidate accessions for
a human to (or a QC pipeline to, later) review - this module only
discovers, it never downloads or ingests anything itself.

NCBI asks that automated tools identify themselves and stay under
3 requests/second without an API key - both respected here.
"""

import time
import xml.etree.ElementTree as ET

import requests

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

TOOL_NAME = "EpiMatch"
CONTACT_EMAIL = "ghanshyamyadav1107@gmail.com"


def search_geo(keyword, organism="Homo sapiens", max_results=20):
    """
    Returns a list of dicts: {accession, title, gds_id}, for GEO
    Series (GSE) entries matching the keyword and organism.
    """

    query = f'{keyword} AND "{organism}"[Organism] AND gse[Entry Type]'

    search_params = {
        "db": "gds",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
        "tool": TOOL_NAME,
        "email": CONTACT_EMAIL,
    }

    response = requests.get(ESEARCH_URL, params=search_params, timeout=30)
    response.raise_for_status()

    id_list = response.json()["esearchresult"]["idlist"]

    if not id_list:
        return []

    time.sleep(0.4)  # stay under NCBI's rate limit

    summary_params = {
        "db": "gds",
        "id": ",".join(id_list),
        "tool": TOOL_NAME,
        "email": CONTACT_EMAIL,
    }

    summary_response = requests.get(ESUMMARY_URL, params=summary_params, timeout=30)
    summary_response.raise_for_status()

    root = ET.fromstring(summary_response.content)

    results = []

    for doc_sum in root.findall(".//DocSum"):

        gds_id = doc_sum.findtext("Id")
        accession = None
        title = None

        for item in doc_sum.findall("Item"):
            if item.get("Name") == "Accession":
                accession = item.text
            if item.get("Name") == "title":
                title = item.text

        if accession:
            results.append({
                "accession": accession,
                "title": title,
                "gds_id": gds_id,
            })

    return results