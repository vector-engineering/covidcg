#!/usr/bin/env python3
# coding: utf-8

"""Get a list of existing sequences, collect metadata

Do this diff with the newly downloaded sequences
so we only do the costly alignment process on new sequences

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import datetime
import gzip
import hashlib
import json
import lzma
import multiprocessing as mp
import os
import pandas as pd

from functools import partial
from pathlib import Path


def extract_ids_and_hashed_sequences(fasta_file):
    """Extract Accession IDs (entry names) and MD5 hashed
    sequences from a fasta file

    Parameters
    ----------
    fasta_file: str

    Returns
    -------
    out: list of tuples, (Accession ID, MD5 sequence, date)

    """

    out = []

    # Read sequences
    cur_entry = ""
    cur_seq = ""

    # Get the date from the fasta file name, as a string
    file_date = Path(fasta_file).name.replace(".fa.gz", "")

    with gzip.open(fasta_file, "rt") as fp:
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
                    out.append(
                        (
                            cur_entry,
                            hashlib.md5(cur_seq.encode("utf-8")).hexdigest(),
                            file_date,
                        )
                    )

                # Clear the entry and sequence
                cur_entry = line[1:]
                # Ignore anything past the first whitespace
                if cur_entry:
                    cur_entry = cur_entry.split()[0]
                cur_seq = ""

    # print("Read {} entries for file {}".format(len(out), fasta_file))

    return out


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-d", "--data-feed", type=str, required=True, help="Path to feed.json file"
    )

    parser.add_argument(
        "-f",
        "--fasta-file-dir",
        type=str,
        required=True,
        help="Path to fasta file directory",
    )

    parser.add_argument(
        "-m",
        "--metadata-out",
        type=str,
        required=True,
        help="Path to metadata CSV output file",
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

    # Make the fasta file directory, if it hasn't been made yet
    # Snakemake won't make the directory itself, since it's a special
    # directory output
    # Path(args.fasta_file_dir).mkdir(exist_ok=True)

    fasta_files = sorted(Path(args.fasta_file_dir).iterdir())
    fasta_files = [f for f in fasta_files if f.is_file() and ".fa.gz" in f.name]

    # Exclude today's file from the manifest
    # Get today's date in ISO format (YYYY-MM-DD)
    today_str = "1_SARS-CoV-2_" + datetime.date.today().isoformat()
    fasta_files = [f for f in fasta_files if f.name != today_str + ".fa.gz"]

    # Store sequence manifest as a list of tuples
    # (Accession ID, MD5 sequence, date)
    entries = []

    def callback(ret):
        entries.extend(ret)

    def error_callback(e):
        raise Exception(e)

    with mp.get_context("spawn").Pool(processes=max_num_processes) as pool:
        for f in fasta_files:
            # Open the output fasta file for this date chunk
            pool.apply_async(
                partial(extract_ids_and_hashed_sequences, str(f)),
                callback=callback,
                error_callback=error_callback,
            )

        pool.close()
        pool.join()

    manifest = pd.DataFrame(entries, columns=["Accession ID", "sequence_hash", "date"])
    #manifest["date"] = pd.to_datetime(manifest["date"])
    # Sort by date, and drop duplicate Accession IDs, by keeping the last copy
    # (i.e., the latest copy)
    manifest = manifest.sort_values("date", ascending=True).drop_duplicates(
        "Accession ID", keep="last"
    )
    # Convert dataframe to dictionary
    manifest = dict(
        zip(manifest["Accession ID"].values, manifest["sequence_hash"].values)
    )

    # Get fields for each isolate
    fields = []
    with lzma.open(args.data_feed, "xt") as fp_in:
        isolate = json.loads(fp_in.readline().strip())
        for i, key in enumerate(isolate.keys()):
            # Skip the special sequence column
            if key == "sequence":
                continue
            fields.append(key)

    # Store metadata entries as a list of dictionaries, for now
    # We'll wrap it in a pandas DataFrame later for easier serialization
    metadata_df = []

    # Store new entries in this list
    # to write into today's fasta file
    # as list of tuples, (Accession ID, sequence)
    new_entries = []

    with open(args.data_feed, "r") as fp_in:
        line_counter = 0

        for line in fp_in:

            if line_counter > 0 and line_counter % 100000 == 0:
                print("Read {} entries so far".format(line_counter))

            try:
                isolate = json.loads(line.strip())
            except json.JSONDecodeError as err:
                print("ERROR PARSING LINE", line_counter)
                print(line)
                continue

            # Add to metadata list
            metadata_df.append({k: isolate[k] for k in fields})

            accession_id = isolate["covv_accession_id"]
            sequence_hash = hashlib.md5(
                isolate["sequence"].strip().replace("\n", "").encode("utf-8")
            ).hexdigest()

            # Does the sequence not exist in our manifest?
            # Or has the sequence changed since the last download of that sequence?
            # If so, then add to the new_entries list
            if accession_id not in manifest:
                new_entries.append((accession_id, isolate["sequence"]))
            elif manifest[accession_id] != sequence_hash:
                print("SEQUENCE HASH MISMATCH: {}".format(accession_id))
                new_entries.append((accession_id, isolate["sequence"]))

            line_counter += 1

    print("{} new/modified sequences".format(len(new_entries)))

    # Write new sequences
    fasta_out = Path(args.fasta_file_dir) / (today_str + ".fa.gz")
    with gzip.open(fasta_out, "wt") as fp_out:
        for accession_id, seq in new_entries:
            fp_out.write(">" + accession_id + "\n" + seq + "\n")

    # Cast the list of dictionaries (list of metadata entries) into a pandas
    # DataFrame, and then serialize it to disk
    # Do this step since pandas can handle some special serialization options
    # that I didn't want to implement manually (such as wrapping certain strings
    # in double quotes)
    metadata_df = pd.DataFrame(metadata_df, columns=fields)
    metadata_df.to_csv(args.metadata_out, index=False)


if __name__ == "__main__":
    main()
