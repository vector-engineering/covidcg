# coding: utf-8

"""Write reference files

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import json
import pandas as pd

from scripts.fasta import read_fasta_file


def write_reference_files(reference_fasta, reference_json):

    # Write the reference fasta file to json
    # Load the reference sequence
    with open(reference_fasta, "r") as fp:
        lines = fp.readlines()
        ref = read_fasta_file(lines)
        ref_seq = list(ref.values())[0]

    with open(reference_json, "w") as fp:
        fp.write(json.dumps({"ref_seq": ref_seq}))
