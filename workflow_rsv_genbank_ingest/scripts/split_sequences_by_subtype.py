#!/usr/bin/env python3
# coding: utf-8

import argparse
import gzip
import pandas as pd


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

    parser.add_argument("--feed", type=str, required=True, help="Input feed FASTA")
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

    # Read sequences
    cur_entry = ""
    cur_seq = ""

    with open(args.feed, "rt", newline="") as fp_in:
        lines = fp_in.readlines()

        for i, line in enumerate(lines):
            # Strip whitespace
            line = line.strip()

            # Ignore empty lines that aren't the last line
            if not line and i <= (len(lines) - 1):
                continue

            # If not the name of an entry, add this line to the current sequence
            # (some FASTA files will have multiple lines per sequence)
            if line[0] != ">":
                cur_seq = cur_seq + line

            # Start of another entry = end of the previous entry
            if line[0] == ">" or i == (len(lines) - 1):
                # Avoid capturing the first one and pushing an empty sequence
                if cur_entry:
                    accession_id = cur_entry.split("|")[0].strip()
                    if accession_id in metadata.index:
                        subtype = metadata.at[accession_id, "subtype"]
                        if subtype:
                            if subtype == "A":
                                seq_a.write(">{}\n{}\n".format(accession_id, cur_seq))
                            else:
                                seq_b.write(">{}\n{}\n".format(accession_id, cur_seq))

                # Clear the entry and sequence
                cur_entry = line[1:]
                # Ignore anything past the first whitespace
                if cur_entry:
                    cur_entry = cur_entry.split()[0]
                cur_seq = ""

    seq_a.close()
    seq_b.close()

    metadata.loc[metadata["subtype"] == "A"].to_csv(args.out_metadata_A)
    metadata.loc[metadata["subtype"] == "B"].to_csv(args.out_metadata_B)


if __name__ == "__main__":
    main()
