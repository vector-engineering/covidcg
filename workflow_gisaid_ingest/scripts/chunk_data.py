# coding: utf-8

"""Chunk data feed, collect metadata

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import gzip
import json
import multiprocessing as mp
import pandas as pd

from collections import defaultdict
from functools import partial
from pathlib import Path


def write_sequences_day(fasta_out_path, seqs):
    # Mode 'at' is append, in text mode
    with gzip.open(fasta_out_path, "at") as fp_out:
        for seq in seqs:
            fp_out.write(">" + seq[0] + "\n" + seq[1] + "\n")


def chunk_data(data_feed, out_fasta, out_metadata, chunk_size=100000, processes=1):
    """Split up the data feed's individual JSON objects into metadata and fasta files. Chunk the fasta files so that every day we only reprocess the subset of fasta files that have changed. The smaller the chunk size, the more efficient the updates, but the more files on the filesystem.
    On a 48-core workstation with 128 GB RAM, aligning 200 sequences takes about 10 minutes, and this is more acceptable than having to align 1000 sequences, which takes ~1 hour. We end up with hundreds of files, but the filesystem seems to be handling it well.

    Parameters
    ----------
    data_feed: str
        - Path to data feed "JSON" file
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
    output_path = Path(out_fasta)

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

    # Get fields for each isolate
    fields = []
    with open(data_feed, "r") as fp_in:
        isolate = json.loads(fp_in.readline().strip())
        for i, key in enumerate(isolate.keys()):
            # Skip the special sequence column
            if key == "sequence":
                continue
            fields.append(key)

    # Store metadata entries as a list of dictionaries, for now
    # We'll wrap it in a pandas DataFrame later for easier serialization
    metadata_df = []

    def flush_chunk(fasta_by_subm_date):
        with mp.get_context("spawn").Pool(processes=processes) as pool:
            for date, seqs in fasta_by_subm_date.items():
                # Open the output fasta file for this date chunk
                fasta_out_path = str(output_path / (date + ".fa.gz"))
                pool.apply_async(partial(write_sequences_day, fasta_out_path, seqs))

            pool.close()
            pool.join()

    with open(data_feed, "r") as fp_in:

        # Open up the initial fasta file for the first chunk
        fasta_by_subm_date = defaultdict(list)

        line_counter = 0
        for line in fp_in:

            # Flush results if chunk is full
            if chunk_i == chunk_size:
                flush_chunk(fasta_by_subm_date)
                # Reset chunk counter
                chunk_i = 0
                # Reset sequence dictionary
                fasta_by_subm_date = defaultdict(list)

            try:
                isolate = json.loads(line.strip())
            except json.JSONDecodeError as err:
                print("ERROR PARSING LINE", line_counter)
                print(line)
                continue

            # Add to metadata list
            metadata_df.append({k: isolate[k] for k in fields})

            # Store sequence in dictionary
            fasta_by_subm_date[isolate["covv_subm_date"]].append(
                (isolate["covv_accession_id"], isolate["sequence"])
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
    metadata_df = pd.DataFrame(metadata_df, columns=fields)
    metadata_df.to_csv(out_metadata, index=False)
