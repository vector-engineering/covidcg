#!/usr/bin/env python3
# coding: utf-8

"""Write reference JSON files

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json
import pandas as pd

from fasta import read_fasta_file


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reference-fasta", type=str, required=True, help="Reference FASTA file"
    )
    parser.add_argument(
        "--primers-csv", type=str, required=True, help="Primers CSV file"
    )
    parser.add_argument(
        "--reference-json", type=str, required=True, help="Reference JSON file"
    )
    parser.add_argument(
        "--primers-json", type=str, required=True, help="Primers JSON file"
    )
    args = parser.parse_args()

    # Write the reference fasta file to json
    # Load the reference sequence
    with open(args.reference_fasta, "r") as fp:
        lines = fp.readlines()
        ref_seqs = read_fasta_file(lines)

    with open(args.reference_json, "w") as fp:
        fp.write(json.dumps(ref_seqs))

    # Load primers, write to JSON
    primers_df = pd.read_csv(args.primers_csv, comment="#")
    # Only take a subset of the data to kee file sizes down
    primers_df[["Institution", "Name", "Sequence", "Reverse", "Start", "End"]].to_json(
        args.primers_json, orient="records"
    )


if __name__ == "__main__":
    main()
