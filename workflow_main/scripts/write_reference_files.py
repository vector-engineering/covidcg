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

        # Find individual references within each subtype
        reference_paths = [
            d
            for d in sorted((Path(args.reference_path) / subtype).iterdir())
            if d.is_dir()
        ]
        for reference_path in reference_paths:
            reference = reference_path.name
            ref_obj[reference] = {}
            ref_obj[reference]["name"] = reference
            ref_obj[reference]["subtype"] = subtype
            ref_obj[reference]["segments"] = {}
            for segment in args.segments:
                with (reference_path / f"{segment}.fa").open("r") as fp:
                    lines = fp.readlines()
                    records = read_fasta_file(lines)
                    ref_obj[reference]["segments"][segment] = list(records.values())[0]

            # Get description
            description = ""
            if (reference_path / "DESCRIPTION").exists():
                with (reference_path / "DESCRIPTION").open("r") as fp:
                    description = fp.read()
            ref_obj[reference]["description"] = description

    with open(args.reference_json, "w") as fp:
        fp.write(json.dumps(ref_obj, indent=2))

    # Load primers, write to JSON
    primers_df = pd.read_csv(args.primers_csv, comment="#")
    # Only take a subset of the data to kee file sizes down
    primers_df[["Institution", "Name", "Sequence", "Reverse", "Start", "End"]].to_json(
        args.primers_json, orient="records"
    )


if __name__ == "__main__":
    main()
