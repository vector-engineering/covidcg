#!/usr/bin/env python3
# coding: utf-8

"""Using the `get_changed_chunks` function, only copy fasta files which have changed
from the purgatory `fasta_temp` folder to the `fasta_raw` folder. By copying over the files,
it will flag to snakemake that they (and only they - not the others) will need to be
reprocessed and realigned.

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import re
import shutil

from pathlib import Path


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--files", type=str, nargs="*", required=True, help="Files to copy"
    )
    parser.add_argument(
        "--data-folder", type=str, required=True, help="Path to data folder"
    )

    args = parser.parse_args()

    # Make the fasta_raw folder if it doesn't already exist yet
    # snakemake won't do this automatically since no fasta_raw files
    # are explicitly defined in any output
    (Path(args.data_folder) / "fasta_raw").mkdir(exist_ok=True)

    # For each changed chunk (as defined by `get_changed_chunks`),
    # copy over the fasta file from the `fasta_temp` folder to the `fasta_raw` folder
    for chunk in args.files:
        chunk = re.sub(r"\.fa\.gz$", "", Path(chunk).name)
        fasta_temp_path = Path(args.data_folder) / "fasta_temp" / (chunk + ".fa.gz")
        fasta_raw_path = Path(args.data_folder) / "fasta_raw" / (chunk + ".fa.gz")

        shutil.copyfile(fasta_temp_path, fasta_raw_path)
        shutil.copystat(fasta_temp_path, fasta_raw_path)


if __name__ == '__main__':
    main()
