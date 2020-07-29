#!/usr/bin/env python3
# coding: utf-8

"""Get ARTIC Network primer sequences, map them onto the reference sequence,
and get them into the same format as the other primers
"""

import pandas as pd
import numpy as np

from cg_scripts.fasta import read_fasta_file
from cg_scripts.util import reverse_complement


def process_artic_primers(artic_files, artic_prefixes, artic_dates, reference_file):
    # Load the reference sequence
    with open(reference_file, "r") as fp:
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

    # print(artic_df)
    return artic_df

