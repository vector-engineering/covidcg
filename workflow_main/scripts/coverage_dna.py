#!/usr/bin/env python3
# coding: utf-8

import argparse

import pandas as pd
import pysam


def coverage_dna(bam_file):
    """Extract coverage of each sequence on their respective references.
    This assumes the sequence (query) is one contiguous sequence and is
    not split up into segments (e.g., paired-end reads that do not overlap).
    Not sure how to adapt to that case without having to crawl through the
    CIGAR/MD string...

    Parameters
    ----------
    bam_file: str
        Path to BAM file

    Returns
    -------
    coverage_df: pandas.DataFrame
    """

    bamfile = pysam.AlignmentFile(bam_file, "r")  # pylint: disable=no-member

    coverage_df = []
    for read in bamfile.fetch(until_eof=True):
        # Skip if unmapped
        if read.is_unmapped:
            continue

        coverage_df.append(
            (
                read.query_name,
                read.reference_name,
                read.reference_start,
                read.reference_end,
            )
        )

    bamfile.close()

    coverage_df = pd.DataFrame.from_records(
        coverage_df, columns=["Accession ID", "reference", "start", "end"]
    )

    return coverage_df


def main():
    """Entry-point"""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bam", type=str, required=True, help="Path to BAM file",
    )
    parser.add_argument("--out", type=str, required=True, help="Path to output file")
    args = parser.parse_args()

    coverage_df = coverage_dna(args.bam)
    coverage_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
