#!/usr/bin/env python3
# coding: utf-8

import argparse
import pandas as pd
import json

# import numpy as np


def bars_overlap(df, bar_row, segments):
    # Determines if a given gene/protein would overlap others in the
    # entropy domain plot
    buffer = 5
    items_in_row = df[df["row"] == bar_row]

    for item, row in items_in_row.iterrows():
        item_segments = []
        if "segments" in df.columns:
            item_segments = df.at[item, "segments"][0]
        else:
            item_segments = df.at[item, "ranges"][0]

        # Flatten ORF1ab segments
        if isinstance(item_segments[0], (list, tuple)):
            item_segments = [item_segments[0][0], item_segments[-1][1]]

        # Segments is the ranges of the bar we are attempting to place.
        # item_segments is the ranges of the bar already placed.
        if (
            (
                segments[0] >= item_segments[0] - buffer
                and segments[0] <= item_segments[1] + buffer
            )
            or (
                segments[1] >= item_segments[0] - buffer
                and segments[1] <= item_segments[1] + buffer
            )
            or (segments[0] > item_segments[0] and segments[1] < item_segments[1])
            or (item_segments[0] > segments[0] and item_segments[1] < segments[1])
        ):
            # Bars overlap
            return True
    # No bars overlap in current row
    return False


def load_genes_or_proteins(file):

    feature_df = pd.read_csv(file, na_values=[""], keep_default_na=False)
    feature_df["domains"] = feature_df["domains"].apply(json.loads)
    feature_df["segments"] = feature_df["segments"].apply(
        lambda x: [list([int(y) for y in z.split("..")]) for z in x.split(";")]
    )
    feature_df["segment"] = feature_df["segment"].apply(str)

    if "notes" in feature_df.columns:
        feature_df["notes"] = feature_df["notes"].fillna("").apply(str)

    file = {}
    subtypes = feature_df["subtype"].unique()
    for subtype in subtypes:
        file[subtype] = {}
        references = feature_df.loc[
            feature_df["subtype"] == subtype, "reference"
        ].unique()
        for reference in references:
            file[subtype][reference] = (
                feature_df.loc[
                    (feature_df["subtype"] == subtype)
                    & (feature_df["reference"] == reference)
                ]
                .drop(columns=["subtype", "reference"])
                .to_dict(orient="records")
            )

    out = dict()

    for subtype in file.keys():
        for reference in file[subtype].keys():

            df = pd.DataFrame.from_records(file[subtype][reference])
            df = df.set_index("name")
            df["protein_coding"] = df["protein_coding"].astype(bool)
            df["segments"] = df["segments"]
            df["len_nt"] = df["segments"].apply(
                lambda x: sum([rng[1] - rng[0] + 1 for rng in x], 0)
            )
            df["len_aa"] = df["len_nt"] // 3
            df["residue_offset_range"] = None
            df["nt_range"] = None
            df.loc[~df["protein_coding"], "len_aa"] = -1
            df.loc[:, "len_aa"] = df["len_aa"].astype(int)

            df["aa_ranges"] = None
            df["nt_ranges"] = None
            df["row"] = None
            for name, row in df.iterrows():
                cur_residue_index = 1

                # Determine rows for entropy plot for All Genes/All Proteins
                bar_row = 0
                segments = df.at[name, "segments"][0]
                # Flatten ORF1ab segments
                if isinstance(segments[0], (list, tuple)):
                    segments = [segments[0][0], segments[-1][1]]

                # If bars would overlap, move down a row and check again
                while bars_overlap(df, bar_row, segments):
                    bar_row += 1

                df.at[name, "row"] = bar_row

                if not row["protein_coding"]:
                    continue

                aa_segments = []
                nt_segments = []
                for rng in row["segments"]:
                    aa_range = [
                        cur_residue_index,
                        cur_residue_index - 1 + (rng[1] - rng[0] + 1) // 3,
                    ]
                    nt_range = [(aa_range[0] * 3) - 2, aa_range[1] * 3]
                    cur_residue_index = aa_range[1] + 1
                    aa_segments.append(aa_range)
                    nt_segments.append(nt_range)

                df.at[name, "aa_ranges"] = aa_segments
                df.at[name, "nt_ranges"] = nt_segments

                # Get the maximal NT extent of this feature on the linear segment
                min_nt = min([x[0] for x in df.at[name, "segments"]])
                max_nt = max([x[1] for x in df.at[name, "segments"]])
                df.at[name, "nt_range"] = [int(min_nt), int(max_nt)]

                # Get the AA extent of this feature, modified by the residue offset
                residue_offset = df.at[name, "residue_offset"]
                df.at[name, "residue_offset_range"] = [
                    int(1 - residue_offset),
                    int(df.at[name, "len_aa"] - residue_offset),
                ]

                # Determine rows for entropy plots for domains
                # All genes/proteins that get here will be protein coding
                domain_df = pd.DataFrame.from_records(row["domains"])
                domain_df["row"] = None
                domain_df["nt_ranges"] = None
                for [index, domain] in domain_df.iterrows():

                    # Domain ranges, in the nucleotide space
                    ranges = domain_df.at[index, "ranges"]
                    nt_ranges = []
                    for rng in ranges:

                        for i, segment in enumerate(aa_segments):
                            # Range is wholly contained within this segment
                            if rng[0] >= segment[0] and rng[1] <= segment[1]:
                                # Adjust with the corresponding segment
                                nt_ranges.append(
                                    [
                                        row["segments"][i][0]
                                        + ((rng[0] - segment[0]) * 3),
                                        row["segments"][i][0]
                                        + ((rng[1] - segment[0] + 1) * 3)
                                        - 1,
                                    ]
                                )
                            # Segment is wholly contained within this range
                            elif rng[0] <= segment[0] and rng[1] >= segment[1]:
                                nt_ranges.append(row["segments"][i])
                            # Middle --> end of segment
                            elif (
                                rng[0] >= segment[0] and rng[0] <= segment[1]
                            ) and rng[1] >= segment[1]:
                                nt_ranges.append(
                                    [
                                        row["segments"][i][0]
                                        + ((rng[0] - segment[0]) * 3),
                                        row["segments"][i][1],
                                    ]
                                )
                            # Beginning --> middle of segment
                            elif rng[0] <= segment[0] and (
                                rng[1] >= segment[0] and rng[1] <= segment[1]
                            ):
                                nt_ranges.append(
                                    [
                                        row["segments"][i][0],
                                        row["segments"][i][0]
                                        + ((rng[1] - segment[0] + 1) * 3)
                                        - 1,
                                    ]
                                )

                    domain_df.at[index, "nt_ranges"] = nt_ranges

                    # Adjust domain ranges with feature residue offset
                    ranges = [
                        [start - residue_offset, end - residue_offset]
                        for start, end in ranges
                    ]
                    domain_df.at[index, "ranges"] = ranges

                    domain_row = 0
                    if (
                        "all" in domain_df.at[index, "name"]
                        or "All" in domain_df.at[index, "name"]
                    ):
                        continue
                    ranges = domain_df.at[index, "ranges"][0]
                    while bars_overlap(domain_df, domain_row, ranges):
                        domain_row += 1
                    domain_df.at[index, "row"] = domain_row

                domain_df.rename(columns={"index": "name"})
                df.at[name, "domains"] = json.loads(domain_df.to_json(orient="records"))

            out[reference] = df.reset_index().to_dict(orient="records")

    return out


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

    out = load_genes_or_proteins(args.input)
    with open(args.output, "w") as fp:
        fp.write(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
