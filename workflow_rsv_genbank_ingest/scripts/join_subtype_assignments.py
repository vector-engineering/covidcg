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
    metadata = metadata.join(serotype_assignments, how="inner")

    # metadata.loc[:, 'serotype'] = metadata['serotype'].combine_first(metadata['assign_serotype'])
    # metadata.drop(columns=['assign_genus', 'assign_serotype'], inplace=True)

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
