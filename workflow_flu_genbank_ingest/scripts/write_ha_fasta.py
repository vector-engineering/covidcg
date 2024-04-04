#!/usr/bin/env python3
# coding: utf-8

import argparse
import gzip
import pandas as pd

def check_validity(metadata, cur_entry):
    # If this entry isn't present in the cleaned metadata, then skip
    accession_id = cur_entry.split("|")[0].strip()

    if accession_id not in metadata.index:
        return False

    segment = metadata.at[accession_id, "segment"]
    if segment != 4:
        return False
    
    return True

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--feed", type=str, required=True, help="Feed fasta file")
    parser.add_argument(
        "--metadata-in", type=str, required=True, help="Path to metadata csv file"
    )
    parser.add_argument(
        "--out", type=str, required=True, help="Output fasta file suffix"
    )

    args = parser.parse_args()

    # Load metadata
    metadata = pd.read_csv(args.metadata_in, index_col="Accession ID")

    write_count = 0

    with open(args.feed, "r", newline="") as fp_in, gzip.open(args.out, "at") as fp_out:
        # Read sequences
        cur_entry = ""
        cur_seq = ""

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
                if cur_entry and check_validity(metadata, cur_entry):
                    accession_id = cur_entry.split("|")[0].strip()
                    fp_out.write(">{}\n{}\n".format(accession_id, cur_seq))
                    write_count += 1
                
                # Clear the entry and sequence
                cur_entry = line[1:]
                # Ignore anything past the first whitespace
                if cur_entry:
                    cur_entry = cur_entry.split()[0]
                cur_seq = ""

    print(f"Wrote {write_count} sequences")


if __name__ == "__main__":
    main()
