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
import multiprocessing as mp
import os
import pandas as pd

from collections import defaultdict
from functools import partial
from pathlib import Path


def flush_chunk(output_path, fasta_by_year_subtype_segment):
    for (year, subtype, segment), seqs in fasta_by_year_subtype_segment.items():
        # print(year, subtype, segment)
        # Open the output fasta file for this year chunk
        fasta_out_path = str(
            output_path / (str(segment) + "_" + subtype + "_" + year + ".fa.gz")
        )
        # Mode 'at' is append, in text mode
        with gzip.open(fasta_out_path, "at", compresslevel=6) as fp_out:
            for seq in seqs:
                fp_out.write(">{}\n{}\n".format(seq[0], seq[1]))


def process_fasta_file(fasta_file, metadata, output_path, chunk_size=10_000):
    """Process raw fasta file, chunk sequences by year, subtype, and segment

    Parameters
    ----------
    fasta_file: str
    metadata: pandas DataFrame
    output_path: str
    chunk_size: int

    Returns
    -------
    dict:
        lines: int
        skips: int

    """

    entries = []

    # Read sequences
    cur_entry = ""
    cur_seq = ""

    # Get the date from the fasta file name, as a string
    # file_date = Path(fasta_file).name.replace(".fa.gz", "")

    with open(fasta_file, "rt", encoding="latin-1") as fp:
        lines = fp.readlines()
        for i, line in enumerate(lines):
            # Strip whitespace
            line = line.strip()

            # Ignore empty lines that aren't the last line
            if not line and i < (len(lines) - 1):
                continue

            # If not the name of an entry, add this line to the current sequence
            # (some FASTA files will have multiple lines per sequence)
            if line[0] != ">":
                cur_seq = cur_seq + line

            # Start of another entry = end of the previous entry
            if line[0] == ">" or i == (len(lines) - 1):
                # Avoid capturing the first one and pushing an empty sequence
                if cur_entry:
                    entries.append(
                        (
                            cur_entry,
                            cur_seq,
                            # file_date,
                        )
                    )

                # Clear the entry and sequence
                cur_entry = line[1:]
                # Ignore anything past the first whitespace
                if cur_entry:
                    cur_entry = cur_entry.split()[0]
                cur_seq = ""

    print("Read {} entries for file {}".format(len(entries), fasta_file))

    # return entries

    # Keep track of how far we're along the current chunk
    chunk_i = 0

    # Open up the initial fasta file for the first chunk
    fasta_by_year_subtype_segment = defaultdict(list)

    line_counter = 0
    skip_counter = 0

    # Create lookup dictionaries
    subtype_lookup = dict(zip(metadata.index, metadata["serotype"]))
    segment_lookup = dict(zip(metadata.index, metadata["segment"]))
    date_lookup = dict(zip(metadata.index, metadata["submission_date"]))

    for name, seq in entries:
        # Flush results if chunk is full
        if chunk_i == chunk_size:
            print("Writing {} sequences".format(chunk_i))
            flush_chunk(output_path, fasta_by_year_subtype_segment)
            # Reset chunk counter
            chunk_i = 0
            # Reset sequence dictionary
            fasta_by_year_subtype_segment = defaultdict(list)

        accession_id = "EPI" + name.split("|")[0]

        # If this entry isn't present in the cleaned metadata, then skip
        if accession_id not in metadata.index:
            skip_counter += 1
            continue

        # subtype = metadata.at[accession_id, "serotype"]
        # segment = metadata.at[accession_id, "segment"]
        # year = metadata.at[accession_id, "submission_date"][0:4]
        subtype = subtype_lookup[accession_id]
        segment = segment_lookup[accession_id]
        year = date_lookup[accession_id][0:4]

        # Store sequence in dictionary (By year, not day)
        fasta_by_year_subtype_segment[(year, subtype, segment)].append(
            (accession_id, seq)
        )

        # Iterate the intra-chunk counter
        chunk_i += 1

        line_counter += 1

    # Flush the last chunk
    print("Writing {} sequences".format(chunk_i))
    flush_chunk(output_path, fasta_by_year_subtype_segment)
    print("Skipped {} sequences".format(skip_counter))

    return {"lines": line_counter, "skips": skip_counter}


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sequences", type=str, nargs="+", required=True, help="Path to sequence files"
    )
    parser.add_argument(
        "--metadata-in", type=str, required=True, help="Path to metadata csv file"
    )
    parser.add_argument(
        "--out-fasta", type=str, required=True, help="Path to fasta output directory"
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10_000,
        help="Number of records to hold in RAM before flushing to disk (default: 10,000)",
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=None,
        help="Number of parallel processes to run. Default: number of CPU cores available",
    )

    args = parser.parse_args()

    # Default to the number of available CPUs in the system
    if args.processes is None or args.processes <= 0:
        max_num_processes = os.cpu_count()
    else:
        max_num_processes = args.processes

    # Load metadata
    metadata = pd.read_csv(args.metadata_in, index_col="Accession ID")

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

    line_counter = []
    skip_counter = []

    def callback(ret):
        line_counter.append(ret["lines"])
        skip_counter.append(ret["skips"])

    def error_callback(e):
        raise Exception(e)

    with mp.get_context("spawn").Pool(processes=max_num_processes) as pool:
        for f in args.sequences:
            # Open the output fasta file for this date chunk
            pool.apply_async(
                partial(
                    process_fasta_file,
                    str(f),
                    metadata,
                    output_path,
                    chunk_size=args.chunk_size,
                ),
                callback=callback,
                error_callback=error_callback,
            )

        pool.close()
        pool.join()

    line_counter = sum(line_counter)
    skip_counter = sum(skip_counter)

    print("Wrote {} sequences".format(line_counter))
    print("Skipped {} sequences".format(skip_counter))


if __name__ == "__main__":
    main()
