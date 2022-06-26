#!/usr/bin/env python3
# coding: utf-8

"""Write reference JSON files

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json
import pandas as pd
from pathlib import Path

from fasta import read_fasta_file


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reference-path", type=str, required=True, help="Reference FASTA file"
    )
    parser.add_argument(
        "--subtypes", type=str, nargs="+", required=True, help="All subtypes"
    )
    parser.add_argument(
        "--segments", type=str, nargs="+", required=True, help="All segments"
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

    ref_obj = {}
    for subtype in args.subtypes:
        ref_obj[subtype] = {}

        # Find individual references within each subtype
        reference_paths = [
            d
            for d in sorted((Path(args.reference_path) / subtype).iterdir())
            if d.is_dir()
        ]
        for reference_path in reference_paths:
            reference = reference_path.name
            ref_obj[subtype][reference] = {}
            ref_obj[subtype][reference]["name"] = reference
            ref_obj[subtype][reference]["segments"] = {}
            for segment in args.segments:
                with (reference_path / f"{segment}.fa").open("r") as fp:
                    lines = fp.readlines()
                    records = read_fasta_file(lines)
                    ref_obj[subtype][reference]["segments"][segment] = list(
                        records.values()
                    )[0]["sequence"]

    with open(args.reference_json, "w") as fp:
        fp.write(json.dumps(ref_obj))

    # Load primers, write to JSON
    primers_df = pd.read_csv(args.primers_csv, comment="#")
    # Only take a subset of the data to kee file sizes down
    primers_df[["Institution", "Name", "Sequence", "Reverse", "Start", "End"]].to_json(
        args.primers_json, orient="records"
    )


if __name__ == "__main__":
    main()
