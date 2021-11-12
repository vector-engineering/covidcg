#!/usr/bin/env python3
# coding: utf-8

"""Get ARTIC Network primer sequences, map them onto the reference sequence,
and get them into the same format as the other primers

This should only be run once - and not part of the pipeline generally

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import os
import numpy as np
import pandas as pd

from scripts.fasta import read_fasta_file
from scripts.util import reverse_complement

reference_sequence_path = os.path.join("../static_data", "reference.fasta")

artic_files = [
    "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V1/nCoV-2019.tsv",
    "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V2/nCoV-2019.tsv",
    "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V3/nCoV-2019.tsv",
]
artic_prefixes = ["V1", "V2", "V3"]
artic_dates = ["2020-01-24", "2020-03-13", "2020-03-20"]

with open(reference_sequence_path, "r") as fp:
    lines = fp.readlines()
    ref = read_fasta_file(lines)
    ref_seq = list(ref.values())[0]

ref_seq_rev = reverse_complement(ref_seq)

artic_df = pd.DataFrame()

for i, f in enumerate(artic_files):
    temp_df = pd.read_csv(f, sep="\t")
    temp_df["version"] = artic_prefixes[i]
    temp_df["Institution"] = "ARTIC Network " + artic_prefixes[i]
    temp_df["Name"] = "ARTIC-" + artic_prefixes[i] + "_" + temp_df["name"]
    temp_df["Date"] = artic_dates[i]
    temp_df["Source"] = f
    artic_df = pd.concat([artic_df, temp_df], ignore_index=True)

artic_df["Start"] = -1
artic_df["End"] = -1
artic_df["Reverse"] = ""

forward_seqs = artic_df["Name"].str.contains("LEFT")
artic_df.loc[forward_seqs, "Start"] = artic_df.loc[forward_seqs, "seq"].apply(
    lambda x: ref_seq.index(x) + 1
)
artic_df.loc[forward_seqs, "End"] = artic_df.loc[forward_seqs, "seq"].apply(
    lambda x: ref_seq.index(x) + len(x)
)

rev_seqs = artic_df["Name"].str.contains("RIGHT")
artic_df.loc[rev_seqs, "Start"] = artic_df.loc[rev_seqs, "seq"].apply(
    lambda x: len(ref_seq_rev) - ref_seq_rev.index(x) - len(x) + 1
)
artic_df.loc[rev_seqs, "End"] = artic_df.loc[rev_seqs, "seq"].apply(
    lambda x: len(ref_seq_rev) - ref_seq_rev.index(x)
)
artic_df.loc[rev_seqs, "Reverse"] = "+"

# print(artic_df)
# print((artic_df["Start"] == -1).sum())
# print((artic_df["End"] == -1).sum())

# Get it into the shape: Institution, Date, Name, Description, Sequence, Reverse, Start, End, Label, Final Conc, Comments, Source

artic_df["Description"] = (
    "version: "
    + artic_df["version"]
    + ", pool: "
    + artic_df["pool"]
    + ", %gc: "
    + artic_df["%gc"].astype(str)
    + ", tm (use 65): "
    + artic_df["tm (use 65)"].astype(str)
    + ", length: "
    + artic_df["length"].astype(str)
)
artic_df = artic_df.rename(columns={"seq": "Sequence"})
artic_df["Label"] = None
artic_df["Final Conc"] = None
artic_df["Comments"] = None

artic_df = artic_df[
    [
        "Institution",
        "Date",
        "Name",
        "Description",
        "Sequence",
        "Reverse",
        "Start",
        "End",
        "Label",
        "Final Conc",
        "Comments",
        "Source",
    ]
]

artic_df.to_csv(os.path.join("../static_data", "artic_primers.csv"))

