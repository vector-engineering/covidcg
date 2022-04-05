#!/usr/bin/env python3
# coding: utf-8

"""Combine metadata and mutation data, create metadata maps

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json
import pandas as pd

from process_mutations import process_mutations


def main():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument("--fasta", type=str, required=True, help="FASTA file directory")
    parser.add_argument("--metadata", type=str, required=True, help="Metadata file")

    parser.add_argument(
        "--dna-mutation-files",
        type=str,
        nargs="+",
        required=True,
        help="DNA mutation files",
    )
    parser.add_argument(
        "--gene-aa-mutation-files",
        type=str,
        nargs="+",
        required=True,
        help="Gene AA mutation files",
    )
    parser.add_argument(
        "--protein-aa-mutation-files",
        type=str,
        nargs="+",
        required=True,
        help="Protein AA mutation files",
    )

    # Output
    parser.add_argument(
        "--metadata-map", type=str, required=True, help="Metadata map file output"
    )
    parser.add_argument(
        "--case-data", type=str, required=True, help="Case data JSON file output"
    )
    parser.add_argument(
        "--case-data-csv", type=str, required=True, help="Case data CSV file output"
    )

    # Parameters
    parser.add_argument(
        "--count-threshold", type=int, default=3, help="Global mutation count threshold"
    )
    parser.add_argument(
        "--metadata-cols", type=str, nargs="+", default=[], help="Metadata columns"
    )

    args = parser.parse_args()

    # Count mutations
    dna_mutation_group_df, dna_mutation_map = process_mutations(
        args.fasta,
        args.dna_mutation_files,
        mode="dna",
        count_threshold=args.count_threshold,
    )
    gene_aa_mutation_group_df, gene_aa_mutation_map = process_mutations(
        args.fasta,
        args.gene_aa_mutation_files,
        mode="gene_aa",
        count_threshold=args.count_threshold,
    )
    protein_aa_mutation_group_df, protein_aa_mutation_map = process_mutations(
        args.fasta,
        args.protein_aa_mutation_files,
        mode="protein_aa",
        count_threshold=args.count_threshold,
    )

    # Load metadata
    df = pd.read_csv(args.metadata).set_index("Accession ID")

    # Exclude sequences without a group assignment
    # (i.e., lineage or clade assignment)
    # "group_cols" is defined in the "group_cols" field in the
    # config.yaml file
    # for col in group_cols:
    #     df = df.loc[~pd.isnull(df[col]), :]

    # Join mutations to main dataframe
    # inner join to exclude filtered out sequences
    df = df.join(
        dna_mutation_group_df[["mutation_id"]],
        on="Accession ID",
        how="inner",
        sort=False,
    ).rename(columns={"mutation_id": "dna_mutation_str"})
    df = df.join(
        gene_aa_mutation_group_df[["mutation_id"]],
        on="Accession ID",
        how="inner",
        sort=False,
    ).rename(columns={"mutation_id": "gene_aa_mutation_str"})
    df = df.join(
        protein_aa_mutation_group_df[["mutation_id"]],
        on="Accession ID",
        how="inner",
        sort=False,
    ).rename(columns={"mutation_id": "protein_aa_mutation_str"})

    # Factorize some more metadata columns
    metadata_maps = {}

    # Metadata cols passed in as kwarg, and defined
    # in config.yaml as "metadata_cols"
    for col in args.metadata_cols:
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
    with open(args.metadata_map, "w") as fp:
        fp.write(json.dumps(metadata_maps))

    # Write final dataframe
    df.to_csv(args.case_data_csv, index_label="Accession ID")
    df.reset_index().to_json(args.case_data, orient="records")


if __name__ == "__main__":
    main()
