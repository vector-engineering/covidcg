# coding: utf-8

"""Combine metadata and SNP data, create metadata maps

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import hashlib
import json
import numpy as np
import pandas as pd

from scripts.process_snps import process_snps


def hash_accession_id(accession_id):
    m = hashlib.sha256()
    m.update(str(accession_id).encode("utf-8"))
    return m.hexdigest()


def combine_all_data(
    # Input
    metadata,
    dna_snp_files,
    gene_aa_snp_files,
    protein_aa_snp_files,
    # Output
    accession_hashmap,
    metadata_map,
    location_map,
    case_data,
    case_data_csv,
    count_threshold=3,
):

    # Count SNPs
    dna_snp_group_df, dna_snp_map = process_snps(
        dna_snp_files, mode="dna", count_threshold=count_threshold
    )
    gene_aa_snp_group_df, gene_aa_snp_map = process_snps(
        gene_aa_snp_files, mode="gene_aa", count_threshold=count_threshold
    )
    protein_aa_snp_group_df, protein_aa_snp_map = process_snps(
        protein_aa_snp_files, mode="protein_aa", count_threshold=count_threshold
    )

    # Load metadata
    df = pd.read_csv(metadata).set_index("Accession ID")

    # Filter out "None" lineages
    df = df.loc[df["lineage"] != "None", :]
    df = df.loc[df["lineage"] != "nan", :]
    # Exclude sequences without a lineage/clade assignment
    df = df.loc[~pd.isnull(df["lineage"]), :]
    df = df.loc[~pd.isnull(df["clade"]), :]

    # Join SNPs to main dataframe
    # inner join to exclude filtered out sequences
    df = df.join(
        dna_snp_group_df[["snp_str"]], on="Accession ID", how="inner", sort=False,
    ).rename(columns={"snp_str": "dna_snp_str"})
    df = df.join(
        gene_aa_snp_group_df[["snp_str"]], on="Accession ID", how="inner", sort=False,
    ).rename(columns={"snp_str": "gene_aa_snp_str"})
    df = df.join(
        protein_aa_snp_group_df[["snp_str"]],
        on="Accession ID",
        how="inner",
        sort=False,
    ).rename(columns={"snp_str": "protein_aa_snp_str"})

    # Semicolon-delimited string to array of SNPs
    df[["dna_snp_str", "gene_aa_snp_str", "protein_aa_snp_str"]] = (
        df[["dna_snp_str", "gene_aa_snp_str", "protein_aa_snp_str"]]
        .astype(str)
        .applymap(lambda x: [int(_x) for _x in x.split(";")])
    )

    # Process location metadata
    # Create complete location column from the separate parts
    # This time with no padding
    loc_str = df["region"].str.cat(
        [
            df["country"].astype(str),
            df["division"].astype(str),
            df["location"].astype(str),
        ],
        sep="/",
    )
    loc_str.name = "loc_str"

    unique_locations = loc_str.drop_duplicates().sort_values(ignore_index=True)

    # Location data is stored in one column, "region/country/division/location"
    location_map_df = pd.concat(
        [
            unique_locations,
            (
                unique_locations.str.split("/", expand=True).iloc[
                    :, :4
                ]  # Only take 4 columns
                # Rename columns
                .rename(
                    columns={0: "region", 1: "country", 2: "division", 3: "location"}
                )
            ),
        ],
        axis=1,
    )

    # Save location map
    location_map_df.drop(columns=["loc_str"]).to_json(location_map, orient="records")

    # Map location IDs back to taxon_locations dataframe
    df["location_id"] = loc_str.map(
        pd.Series(location_map_df.index.values, index=unique_locations.values,)
    )
    # Drop original location columns
    df = df.drop(columns=["region", "country", "division", "location"])

    # Hash Accession IDs. Only take the first 8 chars, that's good enough
    df["hashed_id"] = np.random.rand(len(df))
    df["hashed_id"] = (
        df["hashed_id"].astype(str).apply(hash_accession_id).str.slice(stop=8)
    )
    # Create map of hash -> Accession ID
    accession_hash_df = df[["hashed_id"]]
    accession_hash_df.to_csv(accession_hashmap, index_label="Accession ID")

    # Delete old accession ID column, reassign to hashed ID
    df = (
        df.reset_index()
        .drop(columns=["Accession ID"])
        .rename(columns={"hashed_id": "Accession ID"})
        .set_index("Accession ID")
    )

    # Factorize some more metadata columns
    map_cols = [
        "gender",
        "patient_status",
        "passage",
        "specimen",
        "sequencing_tech",
        "assembly_method",
        "comment_type",
        "authors",
        "originating_lab",
        "submitting_lab",
    ]
    metadata_maps = {}

    for i, col in enumerate(map_cols):
        factor = pd.factorize(df[col])

        id_col = col + "_id"
        df[id_col] = factor[0]

        metadata_maps[col] = pd.Series(factor[1]).to_dict()

    # Drop the original metadata columns
    df = df.drop(columns=map_cols)

    # Add SNP maps into the metadata map
    metadata_maps["dna_snp"] = dna_snp_map.to_dict()
    metadata_maps["gene_aa_snp"] = gene_aa_snp_map.to_dict()
    metadata_maps["protein_aa_snp"] = protein_aa_snp_map.to_dict()

    # Write the metadata map to a JSON file
    with open(metadata_map, "w") as fp:
        fp.write(json.dumps(metadata_maps))

    # Write final dataframe
    df.to_csv(case_data_csv, index_label="Accession ID")
    df.reset_index().to_json(case_data, orient="records")
