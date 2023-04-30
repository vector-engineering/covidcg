#!/usr/bin/env python3
# coding: utf-8

import argparse
import json

import pandas as pd


def coverage_aa(coverage_dna_df, feature_dfs, active_segment):
    """Extract coverage of sequences on the AA level

    Parameters
    ----------
    coverage_dna_df: pandas.DataFrame
        - Accession ID,reference,start,end
    feature_df: dict
        - key: str
            Reference name
        - value: pandas.DataFrame
    active_segment: str
    references: dict

    Returns
    -------
    coverage_aa_df: pandas.DataFrame
        - One row per sequence-reference-feature combo
        - Each row has start and end index in AAs
        - If the entire feature is not covered, then start and end both = null
    """

    coverage_aa_df = []

    for _, row in coverage_dna_df.iterrows():
        reference_name = row["reference"]
        nt_start = row["start"]
        nt_end = row["end"]

        # print(row["Accession ID"], reference_name, nt_start, nt_end)

        # Loop through all features of the sequence's reference
        for feature_name, feature_row in feature_dfs[reference_name].iterrows():

            # Only process genes in this segment/chromosome
            if feature_row["segment"] != active_segment:
                # print(
                #     f'Feature segment {feature_row["segment"]} outside of active segment {active_segment}'
                # )
                continue

            # print(feature_name)

            segments = feature_row["segments"]
            aa_ranges = feature_row["aa_ranges"]

            # Return nulls for both start and end if the NT range
            # is not within this feature's segments
            # This assumes segments are in the order of low -> high index,
            # so only check the first segment for the lower bound and
            # the last segment for the upper bound

            # Check lower bound + one codon of padding OR
            # Check upper bound - one codon of padding
            if (
                nt_end < segments[0][0] + 3
                or nt_start > segments[len(segments) - 1][1] - 3
            ):
                # print(
                #     f"Out of range. sequence: {nt_start}--{nt_end}. "
                #     f"{feature_name}: {segments[0][0]}--{segments[len(segments) - 1][1]} "
                # )
                coverage_aa_df.append(
                    (row["Accession ID"], reference_name, feature_name, None, None)
                )
                continue

            # Loop through segments and collect coverage
            aa_start = None
            aa_end = None
            for segment_i, segment in enumerate(segments):
                # Check to skip segment
                if nt_end < segment[0] + 3 or nt_start > segment[1] - 3:
                    # print(
                    #     "Segment out of range "
                    #     f"sequence: {nt_start}--{nt_end}. "
                    #     f"segment {segment_i}: {segment[0]}--{segment[1]} "
                    # )
                    continue

                # Start is before this segment?
                # Skip if we already have a start from a previous segment
                if nt_start <= segment[0] and aa_start is None:
                    # print(
                    #     f"Start before segment. sequence start: {nt_start}. segment start: {segment[0]}"
                    # )
                    aa_start = aa_ranges[segment_i][0]
                # Start is within this segment?
                # Skip if we already have a start from a previous segment
                elif nt_start > segment[0] and aa_start is None:
                    # print(
                    #     f"Start inside segment. sequence start: {nt_start}. segment start: {segment[0]}"
                    # )
                    # Get the codon index of the NT start within this segment
                    # e.g., nt_start = 4, segment = [1, 9] --> codon #2
                    rel_codon_ind = (nt_start - segment[0]) // 3
                    # Get absolute codon index via. the aa_ranges field
                    # and set it to the AA start
                    aa_start = aa_ranges[segment_i][0] + rel_codon_ind

                # End is after this segment?
                if nt_end > segment[1]:
                    # print(
                    #     f"End after segment. sequence end: {nt_end}. segment end: {segment[1]}"
                    # )
                    aa_end = aa_ranges[segment_i][1]
                # End is within this segment?
                # Don't skip if end is already defined - overwrite previous one instead
                # Again, this assumes segments are defined in order
                # TODO: just sort segments beforehand to ensure this...
                elif nt_end <= segment[1]:
                    # print(
                    #     f"End inside segment. sequence end: {nt_end}. segment end: {segment[1]}"
                    # )
                    # Get the codon index of the NT end within this segment
                    # e.g., nt_end = 8, segment = [1, 9] --> codon #2
                    #       (codon #3 only partially covered)
                    rel_codon_ind = (nt_end - segment[0]) // 3
                    # Get absolute codon index via. the aa_ranges field
                    # and set it to the AA start
                    aa_end = aa_ranges[segment_i][0] + rel_codon_ind
            # END FOR SEGMENT

            # print(aa_start, aa_end)
            # print("")

            # If both start ane end are undefined, then no coverage for this feature
            if aa_start is None and aa_start is None:
                continue
            # If only one of start or end are undefined, something went wrong...
            elif aa_start is None or aa_end is None:
                raise Exception("AA start/end not defined")

            # Push the coverage tuple
            coverage_aa_df.append(
                (
                    row["Accession ID"],
                    reference_name,
                    feature_name,
                    # Apply feature residue offset here
                    # for renumbering residues off of, e.g., signal peptide
                    aa_start - feature_row["residue_offset"],
                    aa_end - feature_row["residue_offset"],
                )
            )
        # END FOR FEATURE
    # END FOR NT COVERAGE ROW

    coverage_aa_df = pd.DataFrame.from_records(
        coverage_aa_df, columns=["Accession ID", "reference", "feature", "start", "end"]
    )

    # print(coverage_aa_df)

    return coverage_aa_df


def main():
    """Entry-point"""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--coverage-dna", type=str, required=True, help="Path to coverage_dna.csv file",
    )
    parser.add_argument(
        "--reference", type=str, required=True, help="Path to reference file"
    )
    parser.add_argument(
        "--gene-protein-def",
        type=str,
        required=True,
        help="Path to gene/protein definition JSON",
    )
    parser.add_argument("--subtype", type=str, required=True, help="Subtype")
    parser.add_argument("--segment", type=str, required=True, help="Segment/chromosome")
    parser.add_argument(
        "--mode", type=str, required=True, help="Mode ('gene' or 'protein')"
    )
    parser.add_argument("--out", type=str, required=True, help="Path to output file")
    args = parser.parse_args()

    # Load NT coverage dataframe
    coverage_dna_df = pd.read_csv(args.coverage_dna)

    # Load the reference sequences
    with open(args.reference, "r") as fp:
        references = json.loads(fp.read())

    # Get references for this subtype
    references = {k: v for k, v in references.items() if v["subtype"] == args.subtype}

    # Load gene/protein defs
    # JSON to dataframe
    with open(args.gene_protein_def, "r") as fp:
        feature_dicts = json.loads(fp.read())

    # Get only features for the above references
    feature_dicts = {k: v for k, v in feature_dicts.items() if k in references.keys()}

    feature_dfs = {}
    for k, v in feature_dicts.items():
        v = pd.DataFrame.from_records(v)

        if args.mode == "gene":
            # Only take protein-coding genes
            v = (
                v.loc[v["protein_coding"] == 1, :]
                # set the gene as the index
                .set_index("name")
            )
        else:
            v = v.set_index("name")

        feature_dfs[k] = v

    coverage_aa_df = coverage_aa(coverage_dna_df, feature_dfs, args.segment)

    coverage_aa_df.insert(1, "subtype", args.subtype)
    coverage_aa_df.insert(2, "segment", args.segment)

    coverage_aa_df.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
