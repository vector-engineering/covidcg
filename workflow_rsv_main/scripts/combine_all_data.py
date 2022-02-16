# coding: utf-8

"""Combine metadata and mutation data, create metadata maps
Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import json
import pandas as pd

from scripts.process_mutations import process_mutations


def combine_all_data(
    # Input
    processed_fasta_files,
    metadata,
    dna_mutation_files,
    gene_aa_mutation_files,
    protein_aa_mutation_files,
    genotypes,
    # Output
    metadata_map,
    case_data,
    case_data_csv,
    # Parameters
    count_threshold=3,
    group_cols=[],
    metadata_cols=[],
):

    # Count mutations
    dna_mutation_group_df, dna_mutation_map = process_mutations(
        processed_fasta_files,
        dna_mutation_files,
        mode="dna",
        count_threshold=count_threshold,
    )
    gene_aa_mutation_group_df, gene_aa_mutation_map = process_mutations(
        processed_fasta_files,
        gene_aa_mutation_files,
        mode="gene_aa",
        count_threshold=count_threshold,
    )
    protein_aa_mutation_group_df, protein_aa_mutation_map = process_mutations(
        processed_fasta_files,
        protein_aa_mutation_files,
        mode="protein_aa",
        count_threshold=count_threshold,
    )

    # Load metadata
    df = pd.read_csv(metadata).set_index("Accession ID")

    # Load genotypes
    genotype_df_list = []
    for file in genotypes:
        genotype_df_list.append(pd.read_csv(file).set_index("Accession ID"))
    genotype_df = pd.concat(genotype_df_list)

    # Exclude sequences without a group assignment
    # (i.e., lineage or clade assignment)
    # "group_cols" is defined in the "group_cols" field in the
    # config.yaml file
    # for col in group_cols:
    #     df = df.loc[~pd.isnull(df[col]), :]

    # Add known genotypes to df
    df = df.join(genotype_df[["genotype", "subtype"]], on="Accession ID", how="inner", sort=False,)

    # Join mutations to main dataframe
    # inner join to exclude filtered out sequences
    df = df.join(
        dna_mutation_group_df[["mutation_id"]], on="Accession ID", how="inner", sort=False,
    ).rename(columns={"mutation_id": "dna_mutation_str"})
    df = df.join(
        gene_aa_mutation_group_df[["mutation_id"]], on="Accession ID", how="inner", sort=False,
    ).rename(columns={"mutation_id": "gene_aa_mutation_str"})
    df = df.join(
        protein_aa_mutation_group_df[["mutation_id"]], on="Accession ID", how="inner", sort=False,
    ).rename(columns={"mutation_id": "protein_aa_mutation_str"})

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

    # Add mutation maps into the metadata map
    metadata_maps["dna_mutation"] = dna_mutation_map.to_dict()
    metadata_maps["gene_aa_mutation"] = gene_aa_mutation_map.to_dict()
    metadata_maps["protein_aa_mutation"] = protein_aa_mutation_map.to_dict()

    # Write the metadata map to a JSON file
    with open(metadata_map, "w") as fp:
        fp.write(json.dumps(metadata_maps))

    # Write final dataframe
    df.to_csv(case_data_csv, index_label="Accession ID")
    df.reset_index().to_json(case_data, orient="records")
