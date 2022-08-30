#!/usr/bin/env python3
# coding: utf-8

import argparse
import csv
import gzip
import pandas as pd
import sys

csv.field_size_limit(sys.maxsize)


def main():

    """
    python3 scripts/split_sequences_by_subtype.py \
            --feed {input.feed} \
            --metadata {input.metadata} \
            --out-A {output.seq_A} \
            --out-B {output.seq_B} \
            --out-metadata-A {output.metadata_A} \
            --out-metadata-B {output.metadata_B}
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--feed", type=str, required=True, help="Input feed CSV")
    parser.add_argument(
        "--metadata", type=str, required=True, help="Input metadata + subtype CSV"
    )
    parser.add_argument("--out-A", type=str, required=True, help="Output A fasta file")
    parser.add_argument("--out-B", type=str, required=True, help="Output B fasta file")
    parser.add_argument(
        "--out-metadata-A", type=str, required=True, help="Output A metadata CSV"
    )
    parser.add_argument(
        "--out-metadata-B", type=str, required=True, help="Output B metadata CSV"
    )

    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, index_col="Accession ID")

    seq_a = gzip.open(args.out_A, "wt")
    seq_b = gzip.open(args.out_B, "wt")

    with open(args.feed, "r", newline="") as fp_in:
        feed_reader = csv.DictReader(fp_in, delimiter=",", quotechar='"')
        for row in feed_reader:

            # If this entry isn't present in the cleaned metadata, then skip
            accession_id = row["genbank_accession"]
            if accession_id not in metadata.index:
                continue

            subtype = metadata.at[accession_id, "subtype"]
            if not subtype:
                continue

            if subtype == "A":
                seq_a.write(">{}\n{}\n".format(accession_id, row["sequence"]))
            else:
                seq_b.write(">{}\n{}\n".format(accession_id, row["sequence"]))

    seq_a.close()
    seq_b.close()

    metadata.loc[metadata["subtype"] == "A"].to_csv(args.out_metadata_A)
    metadata.loc[metadata["subtype"] == "B"].to_csv(args.out_metadata_B)


if __name__ == "__main__":
    main()
