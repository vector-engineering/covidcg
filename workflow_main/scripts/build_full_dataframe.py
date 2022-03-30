#!/usr/bin/env python3
# coding: utf-8

import argparse
import pandas as pd
import json


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--case-data", type=str, required=True, help="Case data JSON file"
    )
    parser.add_argument(
        "--metadata-map", type=str, required=True, help="Metadata map JSON file"
    )
    parser.add_argument("--df-out", type=str, required=True, help="Data frame output")

    args = parser.parse_args()

    # Load data
    df = pd.read_json(args.case_data).set_index("Accession ID")
    with open(args.metadata_map, "r") as fp:
        metadata_map = json.loads(fp.read())

    # Join mutation information
    dna_mutation_map = {v: k for k, v in metadata_map["dna_mutation"].items()}
    gene_aa_mutation_map = {v: k for k, v in metadata_map["gene_aa_mutation"].items()}
    protein_aa_mutation_map = {
        v: k for k, v in metadata_map["protein_aa_mutation"].items()
    }

    df.loc[:, "dna_mutation_str"] = df["dna_mutation_str"].apply(
        lambda x: ";".join([dna_mutation_map[i] for i in x])
    )
    df.loc[:, "gene_aa_mutation_str"] = df["gene_aa_mutation_str"].apply(
        lambda x: ";".join([gene_aa_mutation_map[i] for i in x])
    )
    df.loc[:, "protein_aa_mutation_str"] = df["protein_aa_mutation_str"].apply(
        lambda x: ";".join([protein_aa_mutation_map[i] for i in x])
    )

    df = df.rename(
        columns={
            "dna_mutation_str": "dna_mutation",
            "gene_aa_mutation_str": "gene_aa_mutation",
            "protein_aa_mutation_str": "protein_aa_mutation",
        }
    )

    # Metadata
    for col in metadata_map.keys():
        # We already did these
        if col in ["dna_mutation", "gene_aa_mutation", "protein_aa_mutation"]:
            continue

        mmap = {int(k): v for k, v in metadata_map[col].items()}
        df[col] = df[col].map(mmap)

    df.to_csv(args.df_out)


if __name__ == "__main__":
    main()
