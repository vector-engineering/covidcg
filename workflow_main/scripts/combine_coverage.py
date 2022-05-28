#!/usr/bin/env python3
# coding: utf-8

"""Combine coverage data

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import io
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--manifest", type=str, required=True, help="Path to manifest CSV file"
    )
    parser.add_argument(
        "--coverage-files",
        type=str,
        nargs="+",
        required=True,
        help="DNA coverage files",
    )
    parser.add_argument(
        "--mode", type=str, required=True, help="Mode (dna, gene_aa, protein_aa)"
    )

    # Output
    parser.add_argument(
        "--out", type=str, required=True, help="Combined coverage file CSV output"
    )

    args = parser.parse_args()

    # Get sequence manifest
    manifest = pd.read_csv(args.manifest)

    # Dump all mutation chunks into a text buffer
    coverage_df_io = io.StringIO()
    for i, chunk in enumerate(args.coverage_files):
        file_date = Path(chunk).name.replace("_coverage_" + args.mode + ".csv", "")
        with open(chunk, "r") as fp_in:
            # Write dates, so we can remove duplicate sequences
            # and default to the mutations of the latest sequence, by date
            for j, line in enumerate(fp_in):
                # Write the header of the first file
                if i == 0 and j == 0:
                    coverage_df_io.write(line.strip() + ",date\n")
                # Or write any line that's not the header
                # (to avoid writing the header more than once)
                elif j > 0:
                    coverage_df_io.write(line.strip() + "," + file_date + "\n")

    # Read the buffer into a dataframe, then discard the buffer
    coverage_df_io.seek(0)
    coverage_df = pd.read_csv(coverage_df_io)
    coverage_df_io.close()

    # --------------------------
    # Remove duplicate sequences
    # --------------------------

    # Also has the effect of adding rows for sequences without mutations
    # (pos, ref, alt, etc filled with NaNs)
    coverage_df = manifest.set_index(["Accession ID", "reference", "date"]).join(
        coverage_df.set_index(["Accession ID", "reference", "date"]), how="left"
    )
    coverage_df.reset_index(inplace=True)

    # Remove rows without coverage data
    coverage_df.drop(coverage_df.index[pd.isna(coverage_df["start"])], inplace=True)
    coverage_df.reset_index(drop=True, inplace=True)
    # Convert start/end to ints
    coverage_df[["start", "end"]] = coverage_df[["start", "end"]].astype(int)

    # Drop date column
    coverage_df.drop(columns=["date"], inplace=True)

    coverage_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()

