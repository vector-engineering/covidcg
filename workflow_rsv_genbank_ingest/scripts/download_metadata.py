#!/usr/bin/env python3
# coding: utf-8

"""
(2020-12-21) Derived from:
https://github.com/nextstrain/ncov-ingest: fetch-from-genbank script

License is in LICENSE_NEXTSTRAIN, same folder

------------------

Download all RSV sequences and their curated metadata from GenBank via
NCBI Virus.

Outputs newline-delimited JSON records.

The request this program makes is based on observing the network activity that

    https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Severe%20acute%20respiratory%20syndrome%20coronavirus%202,%20taxid:2697049

performs after clicking through the download interface.  Some tweaks were made
by comparing different download requests and guessing, which allows us to
download the metadata + sequence in the same request instead of two.
"""


import argparse
import datetime
import requests

RETRY_ATTEMPTS = 5


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

    max_date = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S.00Z")
    if args.end_time > max_date:
        args.end_time = max_date

    endpoint = "https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/"
    params = {
        # Search criteria
        "fq": [
            '{!tag=SeqType_s}SeqType_s:("Nucleotide")',  # Nucleotide sequences (as opposed to protein)
            "VirusLineageId_ss:(11250)",  # NCBI Taxon id for RSV
            f"{{!tag=CreateDate_dt}}CreateDate_dt:([{args.start_time} TO {args.end_time}])",
        ],
        # Unclear, but seems necessary.
        "q": "*:*",
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
                ("strain", "Isolate_s"),
                ("region", "Region_s"),
                ("location", "CountryFull_s"),
                ("collected", "CollectionDate_s"),
                ("submitted", "CreateDate_dt"),
                ("length", "SLen_i"),
                ("host", "Host_s"),
                ("isolation_source", "Isolation_csv"),
                ("biosample_accession", "BioSample_s"),
                ("title", "Definition_s"),
                ("authors", "Authors_csv"),
                ("publications", "PubMed_csv"),
                # ("sequence", "Nucleotide_seq"),
                ("protein_names", "ProtNames_ss"),
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

    for _ in range(RETRY_ATTEMPTS):
        try:
            response = requests.get(
                endpoint, params=params, headers=headers, stream=True
            )
            response.raise_for_status()
            break
        except requests.exceptions.RequestException as e:
            print(f"Error downloading metadata: {e}")

        if _ == RETRY_ATTEMPTS - 1:
            raise Exception(
                f"Failed to download metadata after {RETRY_ATTEMPTS} attempts"
            )
        else:
            print(f"Retrying download...")

    response_content = response.iter_content(chunk_size=1024, decode_unicode=True)

    for chunk in response_content:
        print(chunk, end="")


if __name__ == "__main__":
    main()
