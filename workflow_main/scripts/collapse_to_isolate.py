#!/usr/bin/env python3
# coding: utf-8

"""Collapse sequence data by isolate

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument(
        "--sequence-data",
        type=str,
        required=True,
        help="Path to sequence data JSON file",
    )
    # Params
    parser.add_argument(
        "--metadata-cols", type=str, nargs="+", default=[], help="Metadata columns"
    )
    parser.add_argument(
        "--group-cols", type=str, nargs="+", default=[], help="Grouping columns"
    )
    # Output
    parser.add_argument(
        "--isolate-data",
        type=str,
        required=True,
        help="Path to output isolate data JSON file",
    )
    parser.add_argument(
        "--isolate-data-csv",
        type=str,
        required=True,
        help="Path to output isolate data CSV file",
    )

    args = parser.parse_args()

    # Load sequence data
    sequence_df = pd.read_json(args.sequence_data)

    column_aggs = {
        # Sequence metadata lists
        "sequence_ids": ("sequence_id", list),
        "accession_ids": ("Accession ID", list),
        "segments": ("segment", list),
        # Isolate data
        "subtype": ("subtype", "first"),
        "virus_name": ("virus_name", "first"),
        "collection_date": ("collection_date", "first"),
        "submission_date": ("submission_date", "first"),
        # Location data
        "region": ("region", "first"),
        "country": ("country", "first"),
        "division": ("division", "first"),
        "location": ("location", "first"),
        # Mutation data
        "dna_mutation": ("dna_mutation", sum),
        "gene_aa_mutation": ("gene_aa_mutation", sum),
        "protein_aa_mutation": ("protein_aa_mutation", sum),
        # Coverage data
        "dna_range": ("dna_range", list),
        "gene_aa_range": ("gene_aa_range", sum),
        "protein_aa_range": ("protein_aa_range", sum),
    }

    # Metadata/groupings should be the same across all sequences
    # of the same isolate
    for col in args.group_cols + args.metadata_cols:
        column_aggs[col] = (col, "first")

    isolate_df = sequence_df.groupby(["isolate_id", "reference"], as_index=False).agg(
        **column_aggs
    )

    # Create isolate IDs, rename isolate_id to isolate_name
    isolate_df.rename(columns={"isolate_id": "isolate_name"}, inplace=True)
    isolate_ids, _ = pd.factorize(isolate_df["isolate_name"])
    isolate_df.insert(0, "isolate_id", isolate_ids)

    isolate_df.to_json(args.isolate_data, orient="records")
    isolate_df.to_csv(args.isolate_data_csv, index=False)


if __name__ == "__main__":
    main()
