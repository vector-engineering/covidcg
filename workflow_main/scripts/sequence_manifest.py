#!/usr/bin/env python3
# coding: utf-8

"""Sequence manifest (all sequence-reference pairs)

Used as the left side for joining mutation/coverage data

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json
import gzip
from pathlib import Path

import pandas as pd


def extract_ids(fasta_file):
    """Extract Accession IDs (entry names) from a fasta file

    Parameters
    ----------
    fasta_file: str

    Returns
    -------
    out: list of tuples, (Accession ID, date)

    """

    out = []

    # Read sequences
    cur_entry = ""
    cur_seq = ""

    # Get the date from the fasta file name, as a string
    file_name = Path(fasta_file).name.replace(".fa.gz", "")

    with gzip.open(fasta_file, "rt") as fp:
        lines = fp.readlines()
        for i, line in enumerate(lines):
            # Strip whitespace
            line = line.strip()

            # Ignore empty lines that aren't the last line
            if not line and i < (len(lines) - 1):
                continue

            # If not the name of an entry, add this line to the current sequence
            # (some FASTA files will have multiple lines per sequence)
            if line[0] != ">":
                cur_seq = cur_seq + line

            # Start of another entry = end of the previous entry
            if line[0] == ">" or i == (len(lines) - 1):
                # Avoid capturing the first one and pushing an empty sequence
                if cur_entry:
                    out.append((cur_entry, file_name,))

                # Clear the entry and sequence
                cur_entry = line[1:]
                # Ignore anything past the first whitespace
                if cur_entry:
                    cur_entry = cur_entry.split()[0]
                cur_seq = ""

    # print("Read {} entries for file {}".format(len(out), fasta_file))

    return out


def main():
    """Get all sequence-reference pairs
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fasta",
        type=str,
        required=True,
        help="Path to processed FASTA file directory",
    )
    parser.add_argument(
        "--reference", type=str, required=True, help="Path to reference JSON file",
    )
    parser.add_argument(
        "--out", type=str, required=True, help="Output manifest CSV file"
    )

    args = parser.parse_args()

    # Get reference sequence names
    with open(args.reference, "rt") as fp:
        references = json.load(fp)

    subtypes = list(sorted(references.keys()))

    # Get references for each subtype
    subtype_refs = {}
    for subtype in subtypes:
        subtype_refs[subtype] = list(sorted(references[subtype].keys()))

    manifest = []
    for fasta_file in sorted(Path(args.fasta).glob("*.fa.gz")):
        manifest.extend(extract_ids(fasta_file))
    manifest = pd.DataFrame.from_records(
        manifest, columns=["Accession ID", "file_name"]
    )
    pruned_manifest = manifest.drop_duplicates(["Accession ID"], keep="last")

    pruned_manifest = (
        pruned_manifest.reset_index(drop=True)
        .reset_index()
        .rename(columns={"index": "sequence_id"})
    )

    # Add references to manifest
    # Extract subtype and segment from file name
    # Then map subtype references
    pruned_manifest = pruned_manifest.assign(
        segment=lambda x: x["file_name"].str.split("_").str.get(0),
        subtype=lambda x: x["file_name"].str.split("_").str.get(1),
        date=lambda x: x["file_name"].str.split("_").str.get(2),
        reference=lambda x: x["subtype"].map(subtype_refs),
    )

    pruned_manifest = pruned_manifest.explode("reference")

    pruned_manifest.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
