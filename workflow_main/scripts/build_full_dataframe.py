#!/usr/bin/env python3
# coding: utf-8

import argparse
import pandas as pd
import json


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--isolate-data", type=str, required=True, help="Isolate data JSON file"
    )
    parser.add_argument(
        "--metadata-map", type=str, required=True, help="Metadata map JSON file"
    )
    parser.add_argument("--df-out", type=str, required=True, help="Data frame output")

    args = parser.parse_args()

    # Load data
    df = pd.read_json(args.isolate_data).set_index("isolate_id")
    with open(args.metadata_map, "r") as fp:
        metadata_map = json.loads(fp.read())

    # Serialize sequences per isolate
    df.loc[:, "sequence_ids"] = df["sequence_ids"].apply(
        lambda x: ";".join([str(_x) for _x in x])
    )
    df.loc[:, "accession_ids"] = df["accession_ids"].apply(lambda x: ";".join(x))
    df.loc[:, "segments"] = df["segments"].apply(
        lambda x: ";".join([str(_x) for _x in x])
    )

    # Join mutation information
    dna_mutation_map = {v: k for k, v in metadata_map["dna_mutation"].items()}
    gene_aa_mutation_map = {v: k for k, v in metadata_map["gene_aa_mutation"].items()}
    protein_aa_mutation_map = {
        v: k for k, v in metadata_map["protein_aa_mutation"].items()
    }

    df.loc[:, "dna_mutation"] = df["dna_mutation"].apply(
        lambda x: ";".join([dna_mutation_map[i] for i in x])
    )
    df.loc[:, "gene_aa_mutation"] = df["gene_aa_mutation"].apply(
        lambda x: ";".join([gene_aa_mutation_map[i] for i in x])
    )
    df.loc[:, "protein_aa_mutation"] = df["protein_aa_mutation"].apply(
        lambda x: ";".join([protein_aa_mutation_map[i] for i in x])
    )

    # Serialize coverage
    df.loc[:, "dna_range"] = df["dna_range"].apply(
        lambda rngs: ";".join([f"{rng[0]}-{rng[1]}" for rng in rngs])
    )
    df.loc[:, "gene_aa_range"] = df["gene_aa_range"].apply(
        lambda rngs: ";".join([f"{rng[0]}:{rng[1]}-{rng[2]}" for rng in rngs])
    )
    df.loc[:, "protein_aa_range"] = df["protein_aa_range"].apply(
        lambda rngs: ";".join([f"{rng[0]}:{rng[1]}-{rng[2]}" for rng in rngs])
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
