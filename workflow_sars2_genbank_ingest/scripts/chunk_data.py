# coding: utf-8

"""Chunk data feed, collect metadata

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import csv
import datetime
import gzip
import multiprocessing as mp
import pandas as pd
import sys

from collections import defaultdict
from functools import partial
from pathlib import Path


csv.field_size_limit(sys.maxsize)


def write_sequences_day(fasta_out_path, seqs):
    # Mode 'at' is append, in text mode
    with gzip.open(fasta_out_path, "at") as fp_out:
        for seq in seqs:
            fp_out.write(">" + seq[0] + "\n" + seq[1] + "\n")


def main():
    """Split up the data feed's individual objects into metadata and fasta files. Chunk the fasta files so that every day we only reprocess the subset of fasta files that have changed. The smaller the chunk size, the more efficient the updates, but the more files on the filesystem.
    On a 48-core workstation with 128 GB RAM, aligning 200 sequences takes about 10 minutes, and this is more acceptable than having to align 1000 sequences, which takes ~1 hour. We end up with hundreds of files, but the filesystem seems to be handling it well.

    Parameters
    ----------
    data_feed: str
        - Path to data feed csv file
    out_fasta: str
        - Path to fasta output directory
    out_metadata: str
        - Path to metadata.csv output file
    chunk_size: int
        - Number of records to hold in RAM before flushing to disk
    processes: int
        - Number of processes to spawn when writing to disk

    Returns
    -------
    None
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-feed", type=str, required=True, help="Path to data feed CSV file"
    )
    parser.add_argument(
        "--out-fasta", type=str, required=True, help="Path to output fasta directory"
    )
    parser.add_argument(
        "--out-metadata",
        type=str,
        required=True,
        help="Path to output metadata CSV file",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=100000,
        help="Number of records to hold in RAM before flushing to disk",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=1,
        help="Number of processes to spawn when writing to disk",
    )
    args = parser.parse_args()

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

    # Store metadata entries as a list of dictionaries, for now
    # We'll wrap it in a pandas DataFrame later for easier serialization
    metadata_df = []

    def flush_chunk(fasta_by_subm_date):
        with mp.get_context("spawn").Pool(processes=args.processes) as pool:
            for date, seqs in fasta_by_subm_date.items():
                # Open the output fasta file for this date chunk
                fasta_out_path = str(output_path / f"1_SARS-CoV-2_{date}.fa.gz")
                pool.apply_async(partial(write_sequences_day, fasta_out_path, seqs))

            pool.close()
            pool.join()

    with open(args.data_feed, "r", newline="") as fp_in:
        # Open up the initial fasta file for the first chunk
        fasta_by_subm_date = defaultdict(list)

        line_counter = 0

        feed_reader = csv.DictReader(fp_in, delimiter=",", quotechar='"')
        for row in feed_reader:
            # Flush results if chunk is full
            if chunk_i == args.chunk_size:
                flush_chunk(fasta_by_subm_date)
                # Reset chunk counter
                chunk_i = 0
                # Reset sequence dictionary
                fasta_by_subm_date = defaultdict(list)

            # Add to metadata list
            metadata_df.append({k: row[k] for k in row.keys() if k != "sequence"})

            # Store sequence in dictionary
            # Chop off the "Z" at the end of the submission time string, then parse
            # as an ISO datetime format, then return just the year-month-day
            subm_date = datetime.datetime.fromisoformat(row["submitted"][:-1]).strftime(
                "%Y-%m-%d"
            )
            fasta_by_subm_date[subm_date].append(
                (row["genbank_accession"], row["sequence"])
            )

            # Iterate the intra-chunk counter
            chunk_i += 1

            line_counter += 1

        # Flush the last chunk
        flush_chunk(fasta_by_subm_date)

    # Cast the list of dictionaries (list of metadata entries) into a pandas
    # DataFrame, and then serialize it to disk
    # Do this step since pandas can handle some special serialization options
    # that I didn't want to implement manually (such as wrapping certain strings
    # in double quotes)
    metadata_df = pd.DataFrame(metadata_df)
    metadata_df.to_csv(args.out_metadata, index=False)


if __name__ == "__main__":
    main()
