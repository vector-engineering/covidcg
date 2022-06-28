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
        "--coverage-dir", type=str, required=True, help="Coverage files",
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

    # If not in DNA mode, the segment/subtype cols are provided by the
    # coverage files themselves
    # Remove now to prevent duplicate columns during the join
    if args.mode != "dna":
        manifest.drop(columns=["segment", "subtype"], inplace=True)

    # Dump all mutation chunks into a text buffer
    coverage_df_io = io.StringIO()
    for i, chunk in enumerate(sorted(Path(args.coverage_dir).glob("*.csv"))):
        chunk_name = Path(chunk).name.replace("_coverage_" + args.mode + ".csv", "")
        with open(chunk, "r") as fp_in:
            # Write file names, so we can remove duplicate sequences
            # and default to the mutations of the latest sequence, by file name
            for j, line in enumerate(fp_in):
                # Write the header of the first file
                if i == 0 and j == 0:
                    coverage_df_io.write(line.strip() + ",file_name\n")
                # Or write any line that's not the header
                # (to avoid writing the header more than once)
                elif j > 0:
                    coverage_df_io.write(line.strip() + "," + chunk_name + "\n")

    # Read the buffer into a dataframe, then discard the buffer
    coverage_df_io.seek(0)

    # For gene/protein mode, read segment as string
    dtype_opts = {}
    if args.mode == "gene_aa" or args.mode == "protein_aa":
        dtype_opts["segment"] = str
        dtype_opts["feature"] = str

    # Specify only '' for null values
    # Otherwise gene name of 'NA' will be interpreted as missing
    coverage_df = pd.read_csv(
        coverage_df_io, dtype=dtype_opts, keep_default_na=False, na_values=[""],
    )
    coverage_df_io.close()

    # --------------------------
    # Remove duplicate sequences
    # --------------------------

    # Also has the effect of adding rows for sequences without mutations
    # (pos, ref, alt, etc filled with NaNs)

    coverage_df = manifest.merge(
        coverage_df,
        how="left",
        left_on=["Accession ID", "reference", "file_name"],
        right_on=["Accession ID", "reference", "file_name"],
    )

    # Remove rows without coverage data
    coverage_df.drop(coverage_df.index[pd.isna(coverage_df["start"])], inplace=True)
    coverage_df.reset_index(drop=True, inplace=True)
    # Convert start/end to ints
    coverage_df[["start", "end"]] = coverage_df[["start", "end"]].astype(int)

    # Drop unnecessary metadata columns
    coverage_df.drop(columns=["file_name", "date"], inplace=True)

    coverage_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()

