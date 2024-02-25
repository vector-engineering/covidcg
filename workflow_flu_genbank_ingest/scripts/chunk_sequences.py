#!/usr/bin/env python3
# coding: utf-8

"""Chunk sequences in data feed

Split up the data feed's individual objects into metadata and fasta files. 
Chunk the fasta files so that every day we only reprocess the subset of fasta files that have changed. 
The smaller the chunk size, the more efficient the updates, but the more files on the filesystem.

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import gzip
import pandas as pd

from collections import defaultdict
from pathlib import Path


def flush_chunk(output_path, fasta_by_month_serotype_segment):
    for (month, serotype, segment), seqs in fasta_by_month_serotype_segment.items():
        # print(month, serotype, segment)
        # Open the output fasta file for this month chunk
        fasta_out_path = str(
            output_path / (str(segment) + "_" + serotype + "_" + month + ".fa.gz")
        )
        # Mode 'at' is append, in text mode
        with gzip.open(fasta_out_path, "at") as fp_out:
            for seq in seqs:
                fp_out.write(">{}\n{}\n".format(seq[0], seq[1]))


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--feed", type=str, required=True, help="Path to feed FASTA file"
    )
    parser.add_argument(
        "--metadata-in", type=str, required=True, help="Path to metadata csv file"
    )
    parser.add_argument(
        "--metadata-virus-in",
        type=str,
        required=True,
        help="Path to metadata-virus csv file",
    )
    parser.add_argument(
        "--out-fasta", type=str, required=True, help="Path to fasta output directory"
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10_000,
        help="Number of records to hold in RAM before flushing to disk (default: 10,000",
    )

    args = parser.parse_args()

    # Load metadata
    metadata = pd.read_csv(args.metadata_in, index_col="Accession ID")
    metadata_virus = pd.read_csv(args.metadata_virus_in, index_col="virus_name")

    output_path = Path(args.out_fasta)

    # Make the output directory, if it hasn't been made yet
    # Snakemake won't make the directory itself, since it's a special
    # directory output
    if not output_path.exists():
        output_path.mkdir(exist_ok=True)
    else:
        # Erase all files in the output directory
        for fasta_file in output_path.iterdir():
            if fasta_file.is_file():
                fasta_file.unlink()

    # Keep track of how far we're along the current chunk
    chunk_i = 0

    with open(args.feed, "r", newline="") as fp_in:
        # Open up the initial fasta file for the first chunk
        fasta_by_month_serotype_segment = defaultdict(list)

        # Read sequences
        cur_entry = ""
        cur_seq = ""
        entries = []

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
                    entries.append((cur_entry, cur_seq))

                # Clear the entry and sequence
                cur_entry = line[1:]
                # Ignore anything past the first whitespace
                if cur_entry:
                    cur_entry = cur_entry.split()[0]
                cur_seq = ""

        line_counter = 0
        skip_counter = 0

        for cur_entry, cur_seq in entries:
            if chunk_i == args.chunk_size:
                print("Writing {} sequences".format(chunk_i))
                flush_chunk(output_path, fasta_by_month_serotype_segment)
                # Reset chunk counter
                chunk_i = 0
                # Reset sequence dictionary
                fasta_by_month_serotype_segment = defaultdict(list)

            # If this entry isn't present in the cleaned metadata, then skip
            accession_id = cur_entry.split("|")[0].strip()
            if accession_id not in metadata.index:
                skip_counter += 1
                continue

            # Get the serotype
            virus_name = metadata.at[accession_id, "virus_name"]
            if not virus_name or virus_name not in metadata_virus.index:
                skip_counter += 1
                continue

            serotype = metadata_virus.at[virus_name, "serotype"]
            # Get segment
            segment = metadata.at[accession_id, "segment"]

            if not serotype or not segment:
                skip_counter += 1
                continue

            month = metadata.at[accession_id, "submission_date"][0:7]

            # Store sequence in dictionary (By month, not day)
            fasta_by_month_serotype_segment[(month, serotype, segment)].append(
                (accession_id, cur_seq)
            )

            # Iterate the intra-chunk counter
            chunk_i += 1

            line_counter += 1

        # Flush the last chunk
        print("Writing {} sequences".format(chunk_i))
        flush_chunk(output_path, fasta_by_month_serotype_segment)

        print("Skipped {} sequences".format(skip_counter))


if __name__ == "__main__":
    main()
