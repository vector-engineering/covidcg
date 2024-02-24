#!/usr/bin/env python3
# coding: utf-8

"""Collect statistics on sequence quality

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import re


def main():
    """Collect statistics on sequence quality"""

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", type=str, required=True, help="Input FASTA file")
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output CSV file",
    )

    args = parser.parse_args()

    fp_in = open(args.input, "rt")
    fp_out = open(args.output, "wt")

    fp_out.write("Accession ID,length,num_ambiguous\n")

    num_sequences = 0
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

                # Clean up entry ID
                cur_entry = cur_entry.split("|")[0].strip()

                fp_out.write(
                    f'"{cur_entry}",{str(len(cur_seq))},{str(num_ambiguous)}\n'
                )
                num_sequences += 1

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

    print("Processed {:,} sequences".format(num_sequences), flush=True)


if __name__ == "__main__":
    main()
