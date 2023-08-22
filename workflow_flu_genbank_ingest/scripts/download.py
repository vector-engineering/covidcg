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

------------------

From inspecting network data on NCBI virus site
Example record:

{
    "ExportDate_dt":"2021-11-29T07:01:48.204Z",
    "QualNum_i":0,
    "IncompleteCdsCnt_i":0,
    "Host_s":"Homo sapiens",
    "HostSpecies_s":"Homo sapiens (human), taxid:9606|",
    "HostLineage_ss":["cellular organisms, taxid:131567| biota", ...],
    "HostLineageId_ss":["131567", ...],
    "Locus_s":"NC_026431",
    "OrgId_i":641809,
    "VirusFamily_s":"Orthomyxoviridae",
    "VirusGenus_s":"Alphainfluenzavirus",
    "VirusSpecies_s":"Influenza A virus",
    "VirusSpeciesId_i":11320,
    "VirusLineage_ss":["Viruses, taxid:10239| Vira Viridae viruses", ...],
    "VirusLineageId_ss":["10239", ...],
    "VirusL0_s":"RNA viruses",
    "VirusL1_s":"Orthornavirae, taxid:2732396",
    "VirusL2_s":"Negarnaviricota (Negative-strand RNA viruses), taxid:2497569",
    "VirusL3_s":"Polyploviricotina, taxid:2497571",
    "VirusL4_s":"Insthoviricetes, taxid:2497577",
    "VirusL5_s":"Articulavirales, taxid:2499411",
    "VirusL6_s":"Orthomyxoviridae, taxid:11308",
    "VirusL7_s":"Alphainfluenzavirus, taxid:197911",
    "VirusL8_s":"Influenza A virus, taxid:11320",
    "ViralHost_ss":["human",
        "vertebrates"],
    "GenomicMoltype_s":"ssRNA(-)",
    "SLen_i":982,
    "Flags_ss":["refseq",
        "complete"],
    "Flags_csv":"refseq, complete",
    "FlagsCount_i":2,
    "SetAcc_s":"GCF_001343785.1",
    "Strain_s":"A/California/07/2009",
    "Authors_ss":[ ... ],
    "Authors_csv":"...",
    "AuthorsCount_i":68,
    "Country_s":"USA",
    "Segment_s":"7",
    "Division_s":"VRL",
    "Keywords_ss":["RefSeq"],
    "KeywordsCount_i":1,
    "TaxName_s":"Influenza A virus (A/California/07/2009(H1N1))",
    "Segments_ss":["7"],
    "SegmentsCount_i":1,
    "Serotype_s":"H1N1",
    "Region_s":"North America",
    "ParentAcc_s":"set:2d58cc4fd08e5bef66b3d72d7545a9531dad1057",
    "Segmented_s":"true",
    "SetPosition_i":6,
    "SourceDB_s":"RefSeq",
    "USAState_s":"CA",
    "Definition_s":"Influenza A virus (A/California/07/2009(H1N1)) segment 7 matrix protein 2 (M2) and matrix protein 1 (M1) genes, complete cds",
    "HostId_i":9606,
    "CreateDate_dt":"2015-02-23T00:00:00Z",
    "MolType_s":"cRNA",
    "ProtAcc_ss":["YP_009118622",
        "YP_009118623"],
    "ProtAccCount_i":2,
    "UpdateDate_dt":"2018-08-13T00:00:00Z",
    "PubMed_ss":["19465683",
        "19423869"],
    "PubMed_csv":"19465683, 19423869",
    "PubMedCount_i":2,
    "Completeness_s":"complete",
    "CountryFull_s":"USA: California state",
    "ProtNames_ss":["matrix protein 2",
        "matrix protein 1"],
    "ProtNamesCount_i":2,
    "NuclAcc_ss":["NC_026431"],
    "NuclAccCount_i":1,
    "CollectionDate_dr":"2009-04-09",
    "BioProject_ss":["PRJNA485481"],
    "BioProject_csv":"PRJNA485481",
    "BioProjectCount_i":1,
    "AccVer_s":"NC_026431.1",
    "CollectionDate_s":"2009-04-09",
    "GenomeCompleteness_s":"complete",
    "BioProject_s":"PRJNA485481",
    "AccNV_s":"NC_026431",
    "id":"NC_026431",
    "SeqType_s":"Nucleotide",
    "FastaMD5_s":"6c90e1a0ccdb70640f8df30b331aa239",
    "ids_ss":["GCF_001343785",
        "GCF_001343785.1",
        "NC_026431",
        "NC_026431.1",
        "PRJNA485481",
        "YP_009118622",
        "YP_009118623",
        "set:2d58cc4fd08e5bef66b3d72d7545a9531dad1057"],
    "gi_i":758899349,
    "Genome_js":[{"id": "NC_026431.1", "segment": "7", "proteins": [{"id": "YP_009118622.1", "name": "matrix protein 2", "location": "join(1..26,715..982)"}, {"id": "YP_009118623.1", "name": "matrix protein 1", "location": "1..759"}]}]},
}

"""

import argparse
import requests


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--start-time",
        type=str,
        required=True,
        help="Start time in format YYYY-MM-DDTHH:MM:SS.00Z",
    )
    parser.add_argument(
        "--end-time",
        type=str,
        required=True,
        help="Start time in format YYYY-MM-DDTHH:MM:SS.00Z",
    )
    args = parser.parse_args()

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
            f"{{!tag=CreateDate_dt}}CreateDate_dt:([{args.start_time} TO {args.end_time}])",
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
                # AccNV_s - same as ID ("OR450050")
                # AccVer_s - ID + version ("OR450050.1")
                ("database", "SourceDB_s"),
                ("set_id", "SetAcc_s"),
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
                # Completeness_s - same as GenomeCompleteness_s?
                ("segments", "Segments_ss"),
                # ("protein_names", "ProtNames_ss"),
                ("genome_coverage", "Genome_js"),
                # Geographical info
                ("region", "Region_s"),
                ("country", "Country_s"),
                ("location", "CountryFull_s"),
                # Date info
                ("collected", "CollectionDate_s"),
                # CollectionDate_dr - date range
                # CollectionDate_dt - date as datetime
                # CollectionYear_i - year
                ("submitted", "CreateDate_dt"),
                # CreateYear_i - year
                ("updated", "UpdateDate_dt"),
                # Additional metadata
                ("host", "Host_s"),
                ("isolation_source", "Isolation_csv"),
                ("biosample_accession", "BioSample_s"),
                ("title", "Definition_s"),
                # AuthorsCount_i - number of authors
                ("authors", "Authors_csv"),
                # Authors_ss - array of authors
                ("publications", "PubMed_csv"),
                ("sequence", "Nucleotide_seq"),
                # FastaMD5_s - MD5 of sequence
                # BioProjectCount_i - number of BioProjects
                # BioProject_csv - list of BioProjects - csv
                # BioProject_s - list of BioProjects - string
                # BioProject_ss - list of BioProjects - array
                # Division_s - ???
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


if __name__ == "__main__":
    main()
