#!/usr/bin/env python3
# coding: utf-8

"""Join serotype assignments to metadata

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json
import pandas as pd
import pysam


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--in-bam", type=str, required=True, help="Input BAM file")
    parser.add_argument("--metadata", type=str, required=True, help="Metadata CSV file")
    parser.add_argument(
        "--metadata-virus", type=str, required=True, help="Metadata virus CSV file"
    )
    parser.add_argument(
        "--out-metadata", type=str, required=True, help="Output metadata CSV file"
    )
    parser.add_argument(
        "--out-metadata-virus",
        type=str,
        required=True,
        help="Output metadata virus CSV file",
    )

    args = parser.parse_args()

    all_serotypes = []

    bamfile = pysam.AlignmentFile(
        args.in_bam, "r", check_sq=False
    )  # pylint: disable=no-member

    for read in bamfile.fetch(until_eof=True):
        # print(read.query_name, read.reference_name)
        template = read.reference_name.split("/")[0]
        serotype = template.split("_")[0]
        genus = template.split("_")[1]

        accession_id = read.query_name.split("|")[0]

        all_serotypes.append((accession_id, genus, serotype))

    bamfile.close()

    serotype_assignments = pd.DataFrame.from_records(
        all_serotypes, columns=["Accession ID", "genus", "serotype"]
    ).set_index("Accession ID")

    metadata = pd.read_csv(args.metadata, index_col="Accession ID")
    metadata_virus = pd.read_csv(args.metadata_virus, index_col="virus_name")

    # Read Accession ID JSONs
    metadata_virus.loc[:, "accession_ids"] = metadata_virus["accession_ids"].apply(
        json.loads
    )

    # Join virus names
    serotype_assignments = serotype_assignments.merge(
        metadata_virus[["accession_ids"]].explode("accession_ids").reset_index(),
        how="inner",
        left_index=True,
        right_on="accession_ids",
        copy=False,
        sort=False,
    )
    serotype_assignments.drop_duplicates("virus_name", keep="first", inplace=True)

    # Join assignments onto virus DF
    metadata_virus = metadata_virus.join(
        (
            serotype_assignments.drop(columns=["accession_ids"])
            .rename(columns={"genus": "assign_genus", "serotype": "assign_serotype"})
            .set_index("virus_name")
        ),
        how="inner",
    )
    # Join assignments onto sequence DF
    metadata["assign_serotype"] = metadata["isolate_id"].map(
        metadata_virus["assign_serotype"]
    )
    metadata["assign_genus"] = metadata["isolate_id"].map(
        metadata_virus["assign_genus"]
    )
    metadata["serotype"] = metadata["serotype"].fillna(metadata["assign_serotype"])
    metadata.drop(columns=["assign_serotype", "assign_genus"], inplace=True)

    metadata_virus.loc[:, "serotype"] = metadata_virus["serotype"].combine_first(
        metadata_virus["assign_serotype"]
    )
    metadata_virus.drop(columns=["assign_genus", "assign_serotype"], inplace=True)

    # Remove viruses without serotype assignment
    metadata_virus.drop(
        metadata_virus.index[metadata_virus["serotype"].isna()], inplace=True
    )
    metadata.drop(metadata.index[metadata["serotype"].isna()], inplace=True)

    # Rename some serotypes
    # serotype_rename_map = {
    #     'Yamagata': 'B-yam',
    #     'Yamagata-like': 'B-yam',
    #     'Victoria': 'B-vic'
    # }
    # metadata_virus['serotype'] = metadata_virus['serotype'].replace(serotype_rename_map)

    # Drop serotypes not in approved list
    valid_serotypes = ["H1N1", "H3N2", "H5NX", "H7NX", "H9NX", "B-yam", "B-vic"]
    metadata_virus.drop(
        metadata_virus.index[~metadata_virus["serotype"].isin(valid_serotypes)],
        inplace=True,
    )
    metadata.drop(
        metadata.index[~metadata["serotype"].isin(valid_serotypes)], inplace=True
    )

    # Save to disk
    # First reserialize accession IDs
    metadata_virus.loc[:, "accession_ids"] = metadata_virus["accession_ids"].apply(
        json.dumps
    )
    metadata_virus.to_csv(args.out_metadata_virus)
    metadata.to_csv(args.out_metadata)


if __name__ == "__main__":
    main()
