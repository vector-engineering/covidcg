# coding: utf-8

"""Write reference files

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import json
import pandas as pd


def write_reference_files(reference_fasta, primers_csv, reference_json, primers_json):

    # Write the reference fasta file to json
    # Load the reference sequence
    with open(reference_fasta, "r") as fp:
        lines = fp.readlines()
        ref = read_fasta_file(lines)
        ref_seq = list(ref.values())[0]

    with open(reference_json, "w") as fp:
        fp.write(json.dumps({"ref_seq": ref_seq}))

    # Load primers, write to JSON
    primers_df = pd.read_csv(primers_csv, comment="#")
    # Only take a subset of the data to kee file sizes down
    primers_df[["Institution", "Name", "Sequence", "Reverse", "Start", "End"]].to_json(
        primers_json, orient="records"
    )
