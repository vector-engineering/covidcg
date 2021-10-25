# coding: utf-8

"""Combine metadata and SNP data, create metadata maps

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import json
import numpy as np
import pandas as pd

from pathlib import Path

from scripts.process_snps import process_snps


def combine_all_data(
    # Input
    processed_fasta_files,
    metadata,
    dna_snp_files,
    gene_aa_snp_files,
    protein_aa_snp_files,
    # Output
    metadata_map,
    case_data,
    case_data_csv,
    # Parameters
    count_threshold=3,
    group_cols=[],
    metadata_cols=[],
):

    # Count SNPs
    dna_snp_group_df, dna_snp_map = process_snps(
        processed_fasta_files,
        dna_snp_files,
        mode="dna",
        count_threshold=count_threshold,
    )
    gene_aa_snp_group_df, gene_aa_snp_map = process_snps(
        processed_fasta_files,
        gene_aa_snp_files,
        mode="gene_aa",
        count_threshold=count_threshold,
    )
    protein_aa_snp_group_df, protein_aa_snp_map = process_snps(
        processed_fasta_files,
        protein_aa_snp_files,
        mode="protein_aa",
        count_threshold=count_threshold,
    )

    # Load metadata
    df = pd.read_csv(metadata).set_index("Accession ID")

    # Exclude sequences without a group assignment
    # (i.e., lineage or clade assignment)
    # "group_cols" is defined in the "group_cols" field in the
    # config.yaml file
    # for col in group_cols:
    #     df = df.loc[~pd.isnull(df[col]), :]

    # Join SNPs to main dataframe
    # inner join to exclude filtered out sequences
    df = df.join(
        dna_snp_group_df[["snp_id"]], on="Accession ID", how="inner", sort=False,
    ).rename(columns={"snp_id": "dna_snp_str"})
    df = df.join(
        gene_aa_snp_group_df[["snp_id"]], on="Accession ID", how="inner", sort=False,
    ).rename(columns={"snp_id": "gene_aa_snp_str"})
    df = df.join(
        protein_aa_snp_group_df[["snp_id"]], on="Accession ID", how="inner", sort=False,
    ).rename(columns={"snp_id": "protein_aa_snp_str"})

    # Factorize some more metadata columns
    metadata_maps = {}

    # Metadata cols passed in as kwarg, and defined
    # in config.yaml as "metadata_cols"
    for col in metadata_cols:
        factor, labels = pd.factorize(df[col])
        df.loc[:, col] = factor
        metadata_maps[col] = pd.Series(labels).to_dict()

    # Special processing for locations - leave missing data as -1
    for col in ["region", "country", "division", "location"]:
        missing_inds = df[col] == "-1"
        factor, labels = pd.factorize(df.loc[~missing_inds, col])
        df.loc[~missing_inds, col] = factor
        metadata_maps[col] = pd.Series(labels).to_dict()
        df.loc[:, col] = df[col].astype(int)

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
