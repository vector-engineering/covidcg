#!/usr/bin/env python3
# coding: utf-8

"""Filter out SARS2 sequences

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import gzip
import math
import re

import pandas as pd


def main():
    """Filter out sequences (adapted from van Dorp et al, 2020)
    1. Filter against nextstrain exclusion list
    2. Can't be less than 29700NT
    3. Can't have more than 5% ambiguous NT
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", type=str, required=True, help="Input FASTA file")
    parser.add_argument(
        "--nextstrain-exclusion-file",
        type=str,
        required=True,
        help="Path to nextstrain exclusion file",
    )
    parser.add_argument("--exclude-list", type=str, required=True, help="Debug")
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output FASTA file",
    )

    args = parser.parse_args()

    # Load lines, ignoring comments and empty lines
    df = pd.read_csv(
        "https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/exclude.txt",
        comment="#",
        header=None,
        skip_blank_lines=True,
    )

    exclude_taxons = list(set(df[0].tolist()))

    num_excluded = 0
    fp_in = gzip.open(args.input, "rt")
    fp_out = gzip.open(args.output, "wt")
    fp_exclude = open(args.exclude_list, "a")

    cur_entry = ""
    cur_seq = ""
    while True:
        line = fp_in.readline()

        # Beginning of a new entry, or EOF = end of current entry
        if not line or line[0] == ">":
            if cur_entry:
                num_ambiguous = 0
                for char in cur_seq:
                    if char == "N":
                        num_ambiguous += 1

                exclude_reasons = []
                # 1: Check against nextstrain exclusion list
                if cur_entry in exclude_taxons:
                    num_excluded += 1
                    exclude_reasons.append("in_exclusion_list")

                # 2: Can't be less than 29700 NT
                if len(cur_seq) < 29700:
                    num_excluded += 1
                    exclude_reasons.append(f"too_short:{str(len(cur_seq))}")

                # 3: Can't have more than 5% ambiguous (N) NT
                if num_ambiguous > math.floor(len(cur_seq) * 0.05):
                    num_excluded += 1
                    exclude_reasons.append(f"too_many_ambiguous:{str(num_ambiguous)}")

                if len(exclude_reasons) > 1:
                    fp_exclude.write(f"{cur_entry},{';'.join(exclude_reasons)}\n")
                else:
                    # It passed, write to output
                    fp_out.write(">" + cur_entry + "\n")
                    fp_out.write(cur_seq + "\n")

            # If it's the end, then break out
            if not line:
                break

            # Reset sequence and name
            cur_seq = ""
            # Extract the name (up to the first whitespace)
            # [1:] excludes the first '>'
            # .split() breaks up the line into chunks separated by whitespace
            # [0] gets the first chunk
            # cur_entry = line[1:].split()[0]
            # Nevermind, the fasta entries sometimes have spaces.....
            # Just rstrip to remove the newline, that should work good enough
            cur_entry = line[1:].rstrip()

        # Otherwise add sequence to the current entry
        elif cur_entry:
            cur_seq += re.sub(r"\s+", "", line).strip()

    fp_in.close()
    fp_out.close()
    fp_exclude.close()

    print("Removed {:,} sequences".format(num_excluded), flush=True)


if __name__ == "__main__":
    main()
