# coding: utf-8

"""Filter out sequences

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import gzip
import math
import re

<<<<<<< HEAD
=======
from pathlib import Path

>>>>>>> e6dd8312 (Rsvg workflow main (#420))

def preprocess_sequences(input_file, output_file):
    """Filter out sequences (adapted from van Dorp et al, 2020)
    1. Can't have more than 5% ambiguous NT
    """

    num_excluded = 0
    fp_in = gzip.open(input_file, "rt")
    fp_out = gzip.open(output_file, "wt")

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

                if (
                    # 1: Can't have more than 5% ambiguous (N) NT
                    num_ambiguous > math.floor(len(cur_seq) * 0.05)
                ):
                    num_excluded += 1
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

    print("Removed {:,} sequences".format(num_excluded), flush=True)
