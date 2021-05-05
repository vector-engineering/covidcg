#!/usr/bin/env python3
# coding: utf-8

import argparse
import pandas as pd
import numpy as np


def load_genes_or_proteins(file):

    df = pd.read_json(file)
    df = df.set_index("name")
    df["protein_coding"] = df["protein_coding"].astype(bool)
    df["segments"] = df["segments"]
    df["len_nt"] = df["segments"].apply(
        lambda x: sum([rng[1] - rng[0] + 1 for rng in x], 0)
    )
    df["len_aa"] = df["len_nt"] // 3
    df.loc[~df["protein_coding"], "len_aa"] = -1
    df.loc[:, "len_aa"] = df["len_aa"].astype(int)

    df["aa_ranges"] = None
    df["nt_ranges"] = None
    for name, row in df.iterrows():
        cur_residue_index = 1

        if not row["protein_coding"]:
            continue

        aa_segments = []
        nt_segments = []
        for rng in row["segments"]:
            aa_range = [
                cur_residue_index,
                cur_residue_index - 1 + (rng[1] - rng[0] + 1) // 3,
            ]
            nt_range = [
                (aa_range[0] * 3) - 2,
                aa_range[1] * 3
            ]
            cur_residue_index = aa_range[1] + 1
            aa_segments.append(aa_range)
            nt_segments.append(nt_range)

        df.at[name, "aa_ranges"] = aa_segments
        df.at[name, "nt_ranges"] = nt_segments

    return df


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=str,
        help="Path to genes.json or proteins.json input file",
    )
    parser.add_argument(
        "-o", "--output", required=True, type=str, help="Path to output file"
    )

    args = parser.parse_args()

    df = load_genes_or_proteins(args.input)
    df.reset_index().to_json(args.output, orient="records")


if __name__ == "__main__":
    main()
