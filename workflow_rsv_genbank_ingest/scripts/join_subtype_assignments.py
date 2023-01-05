#!/usr/bin/env python3
# coding: utf-8

"""Join RSV subtype assignments to metadata

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import pandas as pd
import pysam


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--in-bam", type=str, required=True, help="Input BAM file")
    parser.add_argument("--metadata", type=str, required=True, help="Metadata CSV file")
    parser.add_argument(
        "--out-metadata", type=str, required=True, help="Output metadata CSV file"
    )
    parser.add_argument(
        "--out-failed",
        type=str,
        required=True,
        help="Output metadata CSV file with only records that did not get assigned a subtype. For debugging purposes only",
    )

    args = parser.parse_args()

    all_subtypes = []

    bamfile = pysam.AlignmentFile(
        args.in_bam, "r", check_sq=False
    )  # pylint: disable=no-member

    for read in bamfile.fetch(until_eof=True):
        # print(read.query_name, read.reference_name)
        template = read.reference_name.split("/")[0]
        subtype = template.split("_")[0]
        genotype = template.split("_")[1]

        accession_id = read.query_name.split("|")[0]

        all_subtypes.append((accession_id, subtype, genotype))

    bamfile.close()

    serotype_assignments = (
        pd.DataFrame.from_records(
            all_subtypes, columns=["Accession ID", "subtype", "genotype"]
        )
        .drop_duplicates("Accession ID")
        .set_index("Accession ID")
    )

    metadata = pd.read_csv(args.metadata, index_col="Accession ID")

    # Join assignments onto virus DF
    metadata = metadata.join(serotype_assignments, how="left")

    # metadata.loc[:, 'serotype'] = metadata['serotype'].combine_first(metadata['assign_serotype'])
    # metadata.drop(columns=['assign_genus', 'assign_serotype'], inplace=True)

    # Special cases for sequencing of non-G genes but subtype is in the virus/isolate name
    # Only making these exceptions for where the name is very clearly deliniating A vs. B
    # TODO: in the GenBank entries for some of these viruses that are don't have A/B in the name,
    #       the GenBank entry has, in the "source" field under "FEATURES", clearly marked "genotype"
    #       or "subgroup" fields. Maybe we can mine and extract these from the NCBI virus API?
    subtype_A_isolate_matches = [
        "RSVA",
        r"^A89",
        r"^A9[0-9]-[0-9]{2,3}-[0-9]{2}",
        # C0607-1015-A
        r"^[CH][0-9]{4}-[0-9]{3,4}-A",
        # RSV/177-A/Nairobi/2008
        r"^RSV\/[0-9]{2,3}-A",
        r"RSVs?\/Tehran\.IRN\/A\/",
        # Ghom/A/12/1363
        r"Ghom\/A\/",
    ]
    subtype_B_isolate_matches = [
        "RSVB",
        r"^B89",
        r"^B9[0-9]-[0-9]{2,3}-[0-9]{2}",
        # B941083-05
        r"^B94[0-9]{4}-[0-9]{2}",
        r"^[CH][0-9]{4}-[0-9]{3,4}-B",
    ]

    for subtype_A_isolate_match in subtype_A_isolate_matches:
        A_match = (
            metadata["virus_name"].str.contains(subtype_A_isolate_match)
            & metadata["subtype"].isna()
        )
        metadata.loc[A_match, "subtype"] = "A"
        metadata.loc[A_match, "genotype"] = "Unknown"

    for subtype_B_isolate_match in subtype_B_isolate_matches:
        B_match = (
            metadata["virus_name"].str.contains(subtype_B_isolate_match)
            & metadata["subtype"].isna()
        )
        metadata.loc[B_match, "subtype"] = "B"
        metadata.loc[B_match, "genotype"] = "Unknown"

    # Save failed subtype assignment to separate file for debugging
    metadata.loc[metadata["subtype"].isna() | metadata["genotype"].isna()].to_csv(
        args.out_failed
    )

    # Remove viruses without serotype assignment
    metadata.drop(metadata.index[metadata["subtype"].isna()], inplace=True)
    metadata.drop(metadata.index[metadata["genotype"].isna()], inplace=True)

    # Rename some serotypes
    subtype_rename_map = {"RSVA": "A", "RSVB": "B"}
    metadata["subtype"] = metadata["subtype"].replace(subtype_rename_map)

    # Save to disk
    metadata.to_csv(args.out_metadata)


if __name__ == "__main__":
    main()
