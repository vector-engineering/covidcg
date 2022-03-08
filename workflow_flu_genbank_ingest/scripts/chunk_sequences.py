#!/usr/bin/env python3
# coding: utf-8

"""Chunk sequences in data feed

Split up the data feed's individual objects into metadata and fasta files. 
Chunk the fasta files so that every day we only reprocess the subset of fasta files that have changed. 
The smaller the chunk size, the more efficient the updates, but the more files on the filesystem.

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import csv
import gzip
import pandas as pd
import sys

from collections import defaultdict
from pathlib import Path


csv.field_size_limit(sys.maxsize)

def flush_chunk(output_path, fasta_by_month_and_segment):
    for (month, segment), seqs in fasta_by_month_and_segment.items():
        # print(month, segment)
        # Open the output fasta file for this month chunk
        fasta_out_path = str(output_path / (str(segment) + '_' + month + ".fa.gz"))
        # Mode 'at' is append, in text mode
        with gzip.open(fasta_out_path, "at") as fp_out:
            for seq in seqs:
                fp_out.write('>{}|{}\n{}\n'.format(seq[0], seq[1], seq[2]))


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--feed', type=str, required=True, help='Path to data feed csv file')
    parser.add_argument('--metadata-in', type=str, required=True, help='Path to metadata csv file')
    parser.add_argument('--out-fasta', type=str, required=True, help='Path to fasta output directory')
    parser.add_argument('--chunk-size', type=int, default=10_000, 
        help='Number of records to hold in RAM before flushing to disk (default: 10,000')

    args = parser.parse_args()

    # Load metadata
    metadata = pd.read_csv(args.metadata_in, index_col='Accession ID')

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
        fasta_by_month_and_segment = defaultdict(list)

        line_counter = 0
        skip_counter = 0

        feed_reader = csv.DictReader(fp_in, delimiter=",", quotechar='"')
        for row in feed_reader:

            # Flush results if chunk is full
            if chunk_i == args.chunk_size:
                print('Writing {} sequences'.format(chunk_i))
                flush_chunk(output_path, fasta_by_month_and_segment)
                # Reset chunk counter
                chunk_i = 0
                # Reset sequence dictionary
                fasta_by_month_and_segment = defaultdict(list)

            # If this entry isn't present in the cleaned metadata, then skip
            accession_id = row['genbank_accession']
            if accession_id not in metadata.index:
                skip_counter += 1
                continue

            # Store sequence in dictionary
            fasta_by_month_and_segment[
                (
                    # By month, not day
                    metadata.at[accession_id, 'submission_date'][0:7], 
                    metadata.at[accession_id, 'segment'])
            ].append(
                (accession_id, metadata.at[accession_id, 'virus_name'], row["sequence"])
            )

            # Iterate the intra-chunk counter
            chunk_i += 1

            line_counter += 1

        # Flush the last chunk
        print('Writing {} sequences'.format(chunk_i))
        flush_chunk(output_path, fasta_by_month_and_segment)

        print('Skipped {} sequences'.format(skip_counter))


if __name__ == '__main__':
    main()