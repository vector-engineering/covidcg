#!/usr/bin/env python3
# coding: utf-8

"""
(2020-12-21) Derived from:
https://github.com/nextstrain/ncov-ingest: fetch-from-genbank script

License is in LICENSE_NEXTSTRAIN, same folder

------------------

Download all SARS-CoV-2 sequences and their curated metadata from GenBank via
NCBI Virus.

Outputs newline-delimited JSON records.

The request this program makes is based on observing the network activity that

    https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Alphainfluenzavirus,%20taxid:197911&VirusLineage_ss=Betainfluenzavirus,%20taxid:197912&VirusLineage_ss=Gammainfluenzavirus,%20taxid:197913

performs after clicking through the download interface.  Some tweaks were made
by comparing different download requests and guessing, which allows us to
download the metadata + sequence in the same request instead of two.
"""

import requests

endpoint = "https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/"
params = {
    # Search criteria
    "fq": [
        '{!tag=SeqType_s}SeqType_s:("Nucleotide")',  # Nucleotide sequences (as opposed to protein)
        # NCBI Taxon id for:
        #   Alphainfluenzavirus: 197911
        #   Betainfluenzavirus: 197912
        #   Gammainfluenzavirus: 197913
        # "VirusLineageId_ss:(197911 OR 197912 OR 197913)",  
        "VirusLineageId_ss:(197911 OR 197912)",  
    ],
    # Unclear, but seems necessary.
    "q": "*:*",
    # "taxresolve": "true", # Maybe necessary for taxonomy/lineage info?
    # Response format
    "cmd": "download",
    "dlfmt": "csv",
    "fl": ",".join(
        ":".join(names)
        for names in [
            # Pairs of (output column name, source data field).  These are pulled
            # from watching requests from the UI.
            #
            # XXX TODO: Is the full set source data fields documented
            # somewhere?  Is there more info we could be pulling that'd be
            # useful?
            #   -trs, 13 May 2020
            ("genbank_accession", "id"),
            ("database", "SourceDB_s"),
            
            # Serotype info is here as well... maybe get these if the serotype field
            # is missing for some of these
            # ("species", "VirusSpecies_s"),
            # ("species_id", "VirusSpeciesId_i"),
            # ("lineage", "VirusLineage_ss"),
            # ("lineage_id", "VirusLineageId_ss"),
            ("genus", "VirusGenus_s"),
            ("serotype", "Serotype_s"),
            ("strain", "Strain_s"),

            # Sequencing/Assembly
            ("length", "SLen_i"),
            ("is_segmented", "Segmented_s"),
            ("complete", "GenomeCompleteness_s"),
            ("segments", "Segments_ss"),
            # ("protein_names", "ProtNames_ss"),
            ("genome_coverage", "Genome_js"),

            # Geographical info
            ("region", "Region_s"),
            ("country", "Country_s"),
            ("location", "CountryFull_s"),

            # Date info
            ("collected", "CollectionDate_s"),
            ("submitted", "CreateDate_dt"),
            ("updated", "UpdateDate_dt"),

            # Additional metadata
            ("host", "Host_s"),
            ("isolation_source", "Isolation_csv"),
            ("biosample_accession", "BioSample_s"),
            ("title", "Definition_s"),
            ("authors", "Authors_csv"),
            ("publications", "PubMed_csv"),
            ("sequence", "Nucleotide_seq"),
        ]
    ),
    # Stable sort with newest last so diffs work nicely.  Columns are source
    # data fields, not our output columns.
    "sort": "SourceDB_s desc, CollectionDate_s asc, id asc",
    # This isn't Entrez, but include the same email parameter it requires just
    # to be nice.
    "email": "covidcg@broadinstitute.org",
}

headers = {
    "User-Agent": "https://github.com/vector-engineering/covid-cg (covidcg@broadinstitute.org)",
}

response = requests.get(endpoint, params=params, headers=headers, stream=True)
response.raise_for_status()

response_content = response.iter_content(chunk_size=1024, decode_unicode=True)

for chunk in response_content:
    print(chunk, end="")
